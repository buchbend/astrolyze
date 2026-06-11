"""The :class:`Stack` container — aligned cutouts gathered from a corpus (#64, PRD #56 stage 1).

A :class:`Stack` is the **stage-1** product of the two-stage stacking design: a container of
sky-coordinate postage-stamp cutouts gathered from every cube that covers a position (or a sample
of sources), built by :meth:`~astrolyze.collection.Collection.stack`. The defining property is
that it is **always safe to construct and browse**, no matter how heterogeneous its members are —
a CO(2-1) single-dish stamp and an HI interferometer stamp sit side by side in the same Stack, and
you can iterate, filter, and :meth:`~Stack.plot_grid` them without any homogeneity precondition
(PRD #56 user stories 10/11/12). Co-addition is a *separate* explicit stage (#65) that gates on
physical homogeneity and is **not** implemented here.

Each member is a :class:`StackMember`: a cutout :class:`~astrolyze.core.Cube` paired with its
**identity** (object / survey / species / transition / beam). The identity is sourced from the
:class:`~astrolyze.collection.Record` it came from — ``survey`` and ``transition`` live on the
catalog row, not on the cube's :class:`~astrolyze.io.Metadata`, so the record is the authority for
them — and falls back to the cube's own metadata for a member without a record. Every member cube
also still carries its own origin provenance in :class:`~astrolyze.io.Metadata` (origin store URI +
catalog version), stamped by :meth:`Record.open`, so a member always traces back to its exact
corpus snapshot.

**Stage 2 (#65) — explicit alignment + homogeneity-gated co-addition.** On top of the stage-1
container, this module adds the *physics* of stacking as **separate, explicit, auditable** steps,
each a per-member broadcast over :meth:`Stack.map`:

- :meth:`Stack.to_common_beam` — convolve every member to a common (largest) beam, reusing the
  existing convolution machinery (:meth:`Cube.convolve_to_beam`) and inheriting its
  no-super-resolution guard (asking for a finer beam raises :class:`~astrolyze.core.cube.LossyDirectionError`).
- :meth:`Stack.to_velocity_grid` — resample every member onto one shared velocity grid via the
  existing regrid machinery (:meth:`Cube.to_velocity_grid`).
- :meth:`Stack.shift_to_rest` — shift each member to rest velocity; ``v_sys`` is **required** and
  resolved per member (catalog → cube metadata → user value), raising if none is available — a
  systemic velocity is never guessed (PRD #56 user story 17).
- :meth:`Stack.coadd` — combine members into one :class:`~astrolyze.core.Cube`. It **raises** unless
  the members are homogeneous (same species / transition / bunit *and* compatible spatial+spectral
  grids), so a physically meaningless average can never be produced silently (user story 16).
  ``weights="noise"`` inverse-variance weights from each member's noise companion; ``"uniform"``
  weights equally. The only path to a valid coadd is the explicit one: ``filter()`` to a
  homogeneous species/transition, then ``to_common_beam()`` / ``to_velocity_grid()`` (and
  ``shift_to_rest()`` when stacking across systemic velocities), then ``coadd()``.

The Stack carries **selection provenance** (:class:`Selection`): what was asked of
:meth:`~astrolyze.collection.Collection.stack` (the position or the source list, the stamp size,
the catalog version), accessible as :attr:`Stack.selection`. :meth:`~Stack.filter` narrows a stack
to a homogeneous subset (``species=`` / ``survey=`` / ``transition=`` / ``object=``) — the path
from a browse-everything stack to a co-addable one (user story 13) — and preserves that provenance.

Design seams these stage-2 methods build on:

- **Homogeneity report.** :attr:`Stack.is_homogeneous` / :meth:`Stack.homogeneity_report` describe
  whether the members share a species / transition / bunit; :meth:`Stack.coadd` turns the same
  report into a **hard precondition** (plus a spatial+spectral grid-compatibility check the
  categorical report does not cover).
- **Per-member broadcast.** :meth:`Stack.map` applies a per-member ``Cube -> Cube`` operation
  across the stack and returns a new Stack; the alignment methods (``to_common_beam`` /
  ``to_velocity_grid`` / ``shift_to_rest``) are exactly such per-member operations and are
  expressed over it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Callable

if TYPE_CHECKING:
    from astropy.coordinates import SkyCoord

    from astrolyze.core import Cube

    from ._facade import Record

#: Fractional tolerance for the coadd grid-compatibility check (spatial WCS world coordinates and
#: the velocity axis). The members are resampled by the alignment methods onto one grid, so a tiny
#: floating-point drift is expected and benign; a real grid mismatch (un-aligned members) is orders
#: of magnitude larger and is what the gate must catch. A pre-1.0 DECISION (issue #65): loose enough
#: to absorb regrid round-off, tight enough that two genuinely different grids never pass.
_GRID_RTOL = 1e-6


# The identity axes a stack can be filtered on / grouped by. These are exactly the catalog-row
# categorical fields a member carries; keeping the list in one place keeps filter() and the
# homogeneity report agreeing on what "the same kind of data" means (DRY).
_IDENTITY_AXES = ("object", "survey", "species", "transition")


@dataclass(frozen=True)
class StackMember:
    """One member of a :class:`Stack`: a cutout cube bound to its identity.

    Pairs the cutout :class:`~astrolyze.core.Cube` (full spectral axis, sub-WCS, origin
    provenance carried on its Metadata) with the identity fields a stack is browsed and filtered
    by. The identity is resolved at gather time from the member's :class:`Record` where it has one
    (``survey`` / ``transition`` come off the catalog row, not the cube Metadata) and from the
    cube's own :class:`~astrolyze.io.Metadata` otherwise (a bare-SkyCoord member with no catalog
    row). :attr:`record` is kept for traceability (``None`` for a record-less member)."""

    cube: "Cube"
    object: str | None
    survey: str | None
    species: str | None
    transition: str | None
    record: "Record | None" = field(default=None, repr=False)

    @property
    def beam(self):
        """The member cube's beam (``radio_beam.Beam`` or ``None``) — read off its Metadata."""
        return self.cube.metadata.beam

    @property
    def bunit(self):
        """The member cube's brightness unit — read off its Metadata (the coadd-gate axis)."""
        return self.cube.metadata.bunit

    @property
    def origin_store_uri(self) -> str | None:
        """The source store URI stamped on the cutout's Metadata (its corpus provenance)."""
        return self.cube.metadata.origin_store_uri

    def resolve_v_sys(self, supplied=None):
        """The systemic velocity to use for this member, by a fixed precedence — or ``None``.

        The resolution order (PRD #56 user story 17; documented on :meth:`Stack.shift_to_rest`):

        1. **the catalog's curated v_sys** (the authority): the member's :class:`Record` row's
           ``v_sys_kms`` where the catalog schema carries it (curated once at publish time);
        2. **the cube's own Metadata** ``systemic_velocity`` (the store schema's optional v_sys);
        3. **a user-supplied** ``supplied`` value (the explicit fallback the user passes to
           :meth:`Stack.shift_to_rest`).

        The catalog/store curated values win over the user value because they are *per source* —
        a single user-supplied scalar cannot be right for a multi-source sample, so a curated value
        must not be silently overridden by it. Returns a velocity ``Quantity`` or ``None`` when no
        v_sys exists at any tier; :meth:`Stack.shift_to_rest` turns the ``None`` into a descriptive
        raise (a systemic velocity is never guessed — ADR-0003)."""
        import astropy.units as u

        # 1) Curated catalog v_sys, when the catalog schema carries it (schema 1.1: v_sys_kms).
        # getattr-guarded: the v1.0/1.1 CatalogRow this build reads may not model the column yet,
        # so a missing field is simply "no curated value", not an error.
        if self.record is not None:
            v_kms = getattr(self.record.row, "v_sys_kms", None)
            if v_kms is not None:
                return float(v_kms) * (u.km / u.s)
        # 2) The cube's own store-schema systemic velocity (an optional Quantity).
        store_v = self.cube.metadata.systemic_velocity
        if store_v is not None:
            return store_v
        # 3) The user-supplied fallback (may itself be None -> the caller raises).
        return supplied

    def _label(self) -> str:
        """A compact human label for a plot panel: object / survey / species (transition)."""
        bits = [b for b in (self.object, self.survey, self.species) if b]
        label = " ".join(bits) if bits else "(member)"
        if self.transition:
            label = f"{label} ({self.transition})"
        return label

    @classmethod
    def from_record(cls, record: "Record", cube: "Cube") -> "StackMember":
        """Build a member from a covering :class:`Record` and the cutout it produced.

        Identity comes from the record's catalog row (the authority for ``survey`` and
        ``transition``); the cube already carries its origin provenance from :meth:`Record.open`."""
        return cls(
            cube=cube,
            object=record.object,
            survey=record.survey,
            species=record.species,
            transition=record.transition,
            record=record,
        )


@dataclass(frozen=True)
class Selection:
    """The selection provenance of a :class:`Stack`: what produced it (PRD #56 user story 19).

    Records the *request* a stack answered — the sky position or the source list it was built for,
    the requested stamp ``size``, and the corpus ``catalog_version`` it ran against — so a stack
    (and any product derived from it) traces back to exactly what was asked of the corpus. It is
    carried unchanged through :meth:`Stack.filter` / :meth:`Stack.map` so a narrowed or
    broadcast stack still names its origin request."""

    #: The kind of request: ``"position"`` (a single SkyCoord) or ``"sources"`` (a list of names
    #: and/or coordinates). Lets a consumer interpret :attr:`targets` without guessing.
    kind: str
    #: The requested target(s), verbatim as passed to ``stack()`` — a single SkyCoord (``kind ==
    #: "position"``) or the list of names / SkyCoords (``kind == "sources"``). Kept for provenance,
    #: not re-resolved.
    targets: object
    #: The requested angular stamp size (the ``size`` argument), verbatim.
    size: object
    #: The ``catalog_schema_version`` of the corpus the stack was gathered from.
    catalog_version: str | None = None


class Stack:
    """A browse-everything container of aligned cutout members (#64; co-addition is #65).

    Built by :meth:`~astrolyze.collection.Collection.stack`; holds an ordered list of
    :class:`StackMember`\\ s (each a cutout :class:`~astrolyze.core.Cube` plus its identity) and the
    :class:`Selection` provenance of the request that produced it. It is a **plain container** at
    this stage — always constructible and browsable regardless of member heterogeneity, with **no**
    homogeneity precondition on construction, iteration, :meth:`filter`, :meth:`map`, or
    :meth:`plot_grid`.

    Sequence protocol: ``len(stack)`` is the member count, ``stack[i]`` is the *i*-th
    :class:`StackMember`, and iterating yields the members in gather order.

    Homogeneity (:attr:`is_homogeneous` / :meth:`homogeneity_report`) is **read-only** here — it
    describes the members for the consumer and as the seam #65's ``coadd`` will gate on, but never
    constrains this stage."""

    def __init__(self, members, selection: Selection):
        # Materialise to a tuple so the container is immutable in member order (a Stack is a
        # snapshot of a gather; filter()/map() return new Stacks rather than mutating in place).
        self._members: tuple[StackMember, ...] = tuple(members)
        self._selection = selection

    # -- selection provenance -----------------------------------------------------------
    @property
    def selection(self) -> Selection:
        """The :class:`Selection` provenance: what request produced this stack (and its subsets)."""
        return self._selection

    # -- the members --------------------------------------------------------------------
    @property
    def members(self) -> tuple[StackMember, ...]:
        """The stack's :class:`StackMember`\\ s in gather order (the axis #65's coadd reduces over)."""
        return self._members

    @property
    def cubes(self) -> list["Cube"]:
        """Just the member cutout :class:`~astrolyze.core.Cube`\\ s, in order (a browsing convenience)."""
        return [member.cube for member in self._members]

    # -- sequence protocol --------------------------------------------------------------
    def __len__(self) -> int:
        return len(self._members)

    def __iter__(self):
        return iter(self._members)

    def __getitem__(self, index):
        # An int indexes a single member; a slice returns a sub-Stack carrying the same selection
        # provenance (slicing is a browse op, not a new gather).
        result = self._members[index]
        if isinstance(index, slice):
            return Stack(result, self._selection)
        return result

    def __repr__(self) -> str:
        return f"<Stack {len(self._members)} member(s); {self._selection.kind}>"

    # -- filtering to a homogeneous subset (user story 13) ------------------------------
    def filter(self, **criteria) -> "Stack":
        """Select the members matching *criteria*; return a new :class:`Stack` (the homogenising step).

        The path from a browse-everything stack to a physically co-addable subset (PRD #56 user
        story 13): ``stack.filter(species="CO", transition="2-1")`` keeps only the CO(2-1) members.
        The filter axes are the identity fields a member carries — ``object`` / ``survey`` /
        ``species`` / ``transition`` — and multiple criteria are **conjunctive** (a member must
        match all). The selection provenance is carried through unchanged, so a filtered stack
        still names the request it descended from.

        An **unknown criterion** raises :class:`ValueError` naming the offending key and the valid
        axes — a silent no-op would hide a typo (``specie=`` vs ``species=``) as a wrongly-empty
        subset (ADR-0003: refuse, never silently mis-match). A criterion that simply matches
        nothing returns an **empty** stack (an honest "no such members", not an error)."""
        unknown = set(criteria) - set(_IDENTITY_AXES)
        if unknown:
            raise ValueError(
                f"unknown stack filter criterion(a): {', '.join(sorted(unknown))}; "
                f"valid axes are {', '.join(_IDENTITY_AXES)}"
            )
        kept = [
            member
            for member in self._members
            if all(getattr(member, key) == value for key, value in criteria.items())
        ]
        return Stack(kept, self._selection)

    # -- per-member broadcast (the #65 alignment seam) ----------------------------------
    def map(self, fn: Callable[["Cube"], "Cube"]) -> "Stack":
        """Apply *fn* to every member cube; return a new :class:`Stack` of the results.

        The broadcast primitive (PRD #56 stage-1 "broadcasting per-member operations"): *fn* is a
        ``Cube -> Cube`` callable applied to each member's cutout in turn, and the results form a
        new Stack with each member's **identity preserved** (only the cube is replaced). The
        selection provenance is carried through unchanged.

        This is the seam #65's alignment methods (``to_common_beam`` / ``to_velocity_grid`` /
        ``shift_to_rest``) hook onto — each is a per-member ``Cube -> Cube`` operation, so they are
        expressed as ``stack.map(...)`` (or thin sugar over it) without this class needing to know
        the physics. *fn* must return a :class:`~astrolyze.core.Cube`; anything else raises, so a
        broken broadcast fails loudly rather than producing a malformed stack (ADR-0003)."""
        from astrolyze.core import Cube

        mapped = []
        for member in self._members:
            result = fn(member.cube)
            if not isinstance(result, Cube):
                raise TypeError(
                    "Stack.map(fn) requires fn to return a Cube (a per-member Cube -> Cube "
                    f"operation); got {type(result).__name__} for member {member._label()!r}"
                )
            mapped.append(
                StackMember(
                    cube=result,
                    object=member.object,
                    survey=member.survey,
                    species=member.species,
                    transition=member.transition,
                    record=member.record,
                )
            )
        return Stack(mapped, self._selection)

    # -- homogeneity: read-only here, the #65 coadd gate --------------------------------
    @property
    def is_homogeneous(self) -> bool:
        """Whether every member shares one species, transition, and brightness unit (the coadd gate).

        ``True`` when the members are physically the *same kind* of data — one species, one
        transition, one bunit — the precondition #65's ``coadd`` will *require* (mixing species or
        units is a meaningless average, PRD #56 user story 16). At **this** stage it is purely
        informational: nothing here enforces it, so a heterogeneous stack stays fully browsable. An
        **empty** stack is reported homogeneous (vacuously — there is nothing inhomogeneous in it).
        Note this does *not* assert spatial/spectral grid alignment — that is the resampling part of
        the #65 gate, decided by the alignment machinery, not a categorical check."""
        return not self.homogeneity_report().conflicts

    def homogeneity_report(self) -> "HomogeneityReport":
        """A structured description of how (in)homogeneous the members are — the coadd-gate input.

        Reports, per categorical axis (species / transition / bunit), the distinct values the
        members carry, and names the axes that have **more than one** value as ``conflicts``. This
        is the read-only seam #65 turns into a hard precondition: stage 1 only *describes*
        homogeneity (so a consumer can see why a stack is not yet co-addable and which
        :meth:`filter` would fix it); #65's ``coadd`` *raises* on any conflict. An empty stack has
        no conflicts (vacuously homogeneous)."""
        species = _distinct_values(m.species for m in self._members)
        transitions = _distinct_values(m.transition for m in self._members)
        bunits = _distinct_values(str(m.bunit) for m in self._members)
        conflicts = tuple(
            axis
            for axis, values in (
                ("species", species),
                ("transition", transitions),
                ("bunit", bunits),
            )
            if len(values) > 1
        )
        return HomogeneityReport(
            species=species,
            transitions=transitions,
            bunits=bunits,
            conflicts=conflicts,
        )

    # -- stage-2 alignment: explicit, auditable resampling (the coadd preconditions) ----
    # Each is a per-member broadcast over map(): the Stack stays thin and the PHYSICS lives in the
    # reused Cube machinery (convolution / regrid / rest-shift), inheriting its guards (notably the
    # no-super-resolution guard). The point of staging them is that every resampling is a VISIBLE
    # call in the analysis, not a silent step buried inside coadd (PRD #56 user story 14).
    def to_common_beam(self, beam=None, *, save_to_tmp_dir: bool = False) -> "Stack":
        """Convolve every member to a common (largest) beam; return a new :class:`Stack`.

        The spatial-alignment step (PRD #56 user story 14): the members are smoothed to one shared
        resolution so they are spatially comparable before co-addition. With *beam* ``None`` the
        target is the **common beam** of the members — the smallest beam *all* of them can be
        *smoothed* to (``radio_beam``'s common beam), never finer than any input, so no member is
        super-resolved. Pass an explicit *beam* to force a specific (larger) target instead.

        This reuses the existing convolution machinery (:meth:`Cube.convolve_to_beam`) and so
        **inherits its no-super-resolution guard**: a member already coarser than the requested
        *beam* would need deconvolution and raises
        :class:`~astrolyze.core.cube.LossyDirectionError` (astrolyze never invents structure,
        ADR-0003). A member already *at* the common beam is left untouched (no needless convolution).

        ``save_to_tmp_dir`` is forwarded to each convolution (the dask-backed eager-materialisation
        control; see :meth:`Cube.convolve_to_beam`). Raises :class:`ValueError` on an empty stack
        (nothing to bring to a common beam) and if any member lacks a beam (a beam cannot be
        guessed). Records the operation in the stack's alignment provenance."""
        if not self._members:
            raise ValueError(
                "to_common_beam() on an empty stack: there are no members to convolve to a "
                "common beam (filter()/gather first)"
            )
        if beam is not None:
            # An EXPLICIT target: convolve every member to it via convolve_to_beam, which raises
            # LossyDirectionError on a member already coarser than the target (the no-super-
            # resolution guard must fire on an explicit finer-beam request, ADR-0003) — never the
            # silent no-op _convolve_if_needed uses for the coarsest member of an auto common beam.
            target = beam
            aligned = self.map(
                lambda cube: cube.convolve_to_beam(
                    target, save_to_tmp_dir=save_to_tmp_dir
                )
            )
        else:
            # The AUTO common beam is >= every member by construction, so the member already at it
            # is a legitimate no-op (not super-resolution); _convolve_if_needed handles that case.
            target = self._common_beam()
            aligned = self.map(
                lambda cube: cube._convolve_if_needed(
                    target, save_to_tmp_dir=save_to_tmp_dir
                )
            )
        return aligned._with_alignment(self._alignment + (f"to_common_beam({target})",))

    def to_velocity_grid(self, grid=None) -> "Stack":
        """Resample every member onto one common velocity *grid*; return a new :class:`Stack`.

        The spectral-alignment step (PRD #56 user story 14): the members are interpolated onto a
        single shared velocity grid so their channels line up voxel-for-voxel before co-addition,
        reusing the existing regrid machinery (:meth:`Cube.to_velocity_grid`). With *grid* ``None``
        the grid is derived from the members — the union velocity range at the **coarsest** member
        channel width (so no member is up-sampled to a finer spectral resolution it never had —
        the spectral analogue of the no-super-resolution guard, ADR-0003). Pass an explicit *grid*
        (a 1-D velocity ``Quantity``) to force one. Target channels outside a member's native
        coverage are filled with ``NaN`` — an honest "not observed", zero-weighted by :meth:`coadd`.

        Raises :class:`ValueError` on an empty stack. Records the operation in the stack's
        alignment provenance."""
        if not self._members:
            raise ValueError(
                "to_velocity_grid() on an empty stack: there are no members to resample"
            )
        target = grid if grid is not None else self._common_velocity_grid()
        aligned = self.map(lambda cube: cube.to_velocity_grid(target))
        return aligned._with_alignment(
            self._alignment + (f"to_velocity_grid(n={len(target)})",)
        )

    def shift_to_rest(self, v_sys=None) -> "Stack":
        """Shift every member to rest velocity (its line to 0 km/s); return a new :class:`Stack`.

        The velocity-alignment step (PRD #56 user stories 14/17): stacking sources at different
        systemic velocities requires bringing each to rest first, so the same channel samples the
        same *rest* velocity across members. Each member's v_sys is **resolved per member**
        (:meth:`StackMember.resolve_v_sys`) by this precedence:

        1. the **catalog's curated v_sys** (``Record.row.v_sys_kms``, where the catalog carries it);
        2. the **cube's own Metadata** ``systemic_velocity`` (the store schema's optional v_sys);
        3. the **user-supplied** *v_sys* passed here (the explicit fallback).

        A systemic velocity is **never guessed**: if a member has none at any tier this **raises**
        :class:`ValueError` naming the member (PRD #56 user story 17 — velocity alignment must never
        rest on a fabricated v_sys, ADR-0003). The shift itself is an exact coordinate relabel (no
        resampling), delegated to :meth:`Cube.shift_to_rest`. Records the operation (and the
        resolved per-member v_sys) in the stack's alignment provenance.

        Note the per-member resolution means a single *v_sys* you pass here only fills members the
        catalog/store did not already curate — it cannot silently override a curated per-source
        value (a scalar cannot be correct for a multi-source sample)."""
        if not self._members:
            raise ValueError(
                "shift_to_rest() on an empty stack: there are no members to shift to rest"
            )
        shifted = []
        resolved = []
        for member in self._members:
            v = member.resolve_v_sys(supplied=v_sys)
            if v is None:
                raise ValueError(
                    f"shift_to_rest(): no systemic velocity for member {member._label()!r} — "
                    "none is curated in the catalog (v_sys_kms), none is on the cube's Metadata "
                    "(systemic_velocity), and none was supplied. astrolyze never guesses a "
                    "systemic velocity (PRD #56 user story 17, ADR-0003); pass v_sys=... or "
                    "curate it in the catalog/store"
                )
            resolved.append((member._label(), v))
            shifted.append(
                StackMember(
                    cube=member.cube.shift_to_rest(v),
                    object=member.object,
                    survey=member.survey,
                    species=member.species,
                    transition=member.transition,
                    record=member.record,
                )
            )
        aligned = Stack(shifted, self._selection)
        return aligned._with_alignment(
            self._alignment + (f"shift_to_rest({len(resolved)} member(s))",)
        )

    # -- stage-2 co-addition: the homogeneity gate + noise-weighted combine -------------
    def coadd(self, weights: str = "noise") -> "Cube":
        """Combine the members into a single :class:`~astrolyze.core.Cube` — the gated co-addition.

        The end of the staged stacking design (PRD #56 user stories 15/16). It **raises** unless
        the members are homogeneous, so a physically meaningless average is impossible to produce
        silently. The gate (:meth:`homogeneity_report` plus a grid check) requires:

        - one **species**, one **transition**, one **bunit** (the categorical conflicts of
          :meth:`homogeneity_report`) — mixing species or units is a meaningless average;
        - a compatible **spatial grid** (same pixel shape + celestial WCS within tolerance) and
          **spectral grid** (same channel count + velocity axis within tolerance) — members not yet
          resampled onto one grid cannot be added voxel-for-voxel.

        A heterogeneous stack raises a descriptive :class:`CoaddError` naming the mismatch and the
        missing alignment step (``filter`` / ``to_common_beam`` / ``to_velocity_grid``). The only
        path to a valid coadd is the explicit one (ADR-0003 — no silent physics).

        Weighting:

        - ``weights="noise"`` (default) — **inverse-variance** weighting from each member's noise
          companion (loaded from its origin store via
          :meth:`~astrolyze.core.NoiseModel.from_zarr_companion`): per voxel ``w_i = 1/σ_i²``, the
          combined value is ``Σ w_i d_i / Σ w_i`` and the combined noise is ``√(1 / Σ w_i)`` (the
          standard inverse-variance result). A voxel that is ``NaN`` in a member (off its coverage)
          or carries no σ contributes zero weight — never a fabricated value. The σ used is the
          member's **as-published** companion (the noise the corpus curated for that store); it is
          the exact per-voxel σ when the members were not spatially smoothed differently on the way
          in, and the right *relative* weighting (the IVW mean depends on the σ *ratios*) otherwise.
          (Known limitation, issue #65: the companion is not yet re-propagated through a preceding
          :meth:`to_common_beam` convolution — so after heterogeneous smoothing the absolute combined
          σ is conservative. The Cube convolution machinery already supports analytic noise
          propagation; threading it through the stack is follow-up work.)
        - ``weights="uniform"`` — equal weights: the combined value is the (NaN-ignoring) mean and,
          when noise companions are available, the combined noise is ``√(Σ σ_i²)/N`` (uncorrelated
          error propagation of an equal-weight mean).

        The result carries combined provenance in its Metadata (which members, the alignment steps,
        the catalog version) via the provenance seam. Raises :class:`CoaddError` on an empty stack
        or an unknown *weights* mode (never a silent default)."""
        if weights not in ("noise", "uniform"):
            raise CoaddError(
                f"unknown coadd weights mode {weights!r}; use 'noise' (inverse-variance from the "
                "noise companions) or 'uniform' (equal). astrolyze never silently substitutes a "
                "weighting (ADR-0003)"
            )
        self._require_coaddable()
        return _coadd_members(self, weights=weights)

    # -- the coadd gate (the no-silent-physics checkpoint) ------------------------------
    def _require_coaddable(self) -> None:
        """Raise :class:`CoaddError` unless the members are homogeneous AND on one grid.

        The hard precondition #65 turns the read-only homogeneity report into, plus the
        spatial+spectral grid-compatibility check the categorical report does not cover. Every
        failure names the offending axis and the alignment step that fixes it (ADR-0003)."""
        if not self._members:
            raise CoaddError(
                "coadd() on an empty stack: there are no members to combine. Gather/filter to a "
                "non-empty homogeneous subset first."
            )
        # 1) Categorical homogeneity (species / transition / bunit) — reuse the read-only report.
        report = self.homogeneity_report()
        if report.conflicts:
            details = []
            if "species" in report.conflicts:
                details.append(f"species {report.species}")
            if "transition" in report.conflicts:
                details.append(f"transition {report.transitions}")
            if "bunit" in report.conflicts:
                details.append(f"bunit {report.bunits}")
            raise CoaddError(
                "coadd() refuses a heterogeneous stack — a physically meaningless average. "
                f"Conflicting axes: {', '.join(details)}. Narrow to one kind of data first, e.g. "
                "stack.filter(species=..., transition=...) (PRD #56 user story 16, ADR-0003)."
            )
        # 2) Spatial + spectral grid compatibility (the resampling half of the gate).
        self._require_aligned_grids()

    def _require_aligned_grids(self) -> None:
        """Raise :class:`CoaddError` unless every member shares one spatial + spectral grid.

        Compares each member against the first: identical pixel shape, a celestial WCS that maps
        the corner pixels to the same sky within tolerance, and a velocity axis equal within
        tolerance. A mismatch is a stack that has not been resampled onto one grid yet, so it names
        the missing alignment step (``to_common_beam`` then ``to_velocity_grid``)."""
        import numpy as np

        reference = self._members[0].cube
        ref_shape = reference.shape
        ref_corners = _celestial_corners(reference)
        ref_velocity = _member_velocity_axis(reference)
        for member in self._members[1:]:
            cube = member.cube
            if cube.shape != ref_shape:
                raise CoaddError(
                    "coadd() requires members on one grid, but their pixel shapes differ "
                    f"({ref_shape} vs {cube.shape} for {member._label()!r}). Align first: "
                    "to_common_beam() then to_velocity_grid() bring the members onto a shared "
                    "spatial+spectral grid (ADR-0003)."
                )
            corners = _celestial_corners(cube)
            if not np.allclose(corners, ref_corners, rtol=_GRID_RTOL, atol=_GRID_RTOL):
                raise CoaddError(
                    f"coadd() requires one spatial grid, but {member._label()!r} maps its pixels "
                    "to a different sky footprint than the first member (the cutouts are not on a "
                    "common spatial grid). astrolyze does not silently reproject during coadd — "
                    "align the spatial grids first (ADR-0003)."
                )
            velocity = _member_velocity_axis(cube)
            if velocity.shape != ref_velocity.shape or not np.allclose(
                velocity, ref_velocity, rtol=_GRID_RTOL, atol=_GRID_RTOL
            ):
                raise CoaddError(
                    f"coadd() requires one spectral grid, but {member._label()!r} is on a "
                    "different velocity axis than the first member. Resample onto one grid first "
                    "with to_velocity_grid() (and shift_to_rest() to stack across systemic "
                    "velocities) — astrolyze never silently regrids inside coadd (ADR-0003)."
                )

    # -- alignment provenance -----------------------------------------------------------
    @property
    def _alignment(self) -> tuple:
        """The ordered alignment steps applied to reach this stack (carried for coadd provenance)."""
        return getattr(self, "_alignment_steps", ())

    def _with_alignment(self, steps: tuple) -> "Stack":
        """A copy of this stack tagging *steps* as its alignment history (provenance only)."""
        clone = Stack(self._members, self._selection)
        clone._alignment_steps = tuple(steps)
        return clone

    # -- common-beam / common-grid derivation (the alignment-target helpers) ------------
    def _common_beam(self):
        """The members' common beam (the smallest beam all of them can be smoothed to)."""
        import radio_beam

        beams = [member.beam for member in self._members]
        if any(b is None for b in beams):
            missing = [m._label() for m in self._members if m.beam is None]
            raise ValueError(
                "to_common_beam() needs every member to carry a beam; these do not: "
                f"{', '.join(missing)} (a beam cannot be guessed, ADR-0003)"
            )
        if len(beams) == 1:
            return beams[0]
        return radio_beam.commonbeam.commonbeam(radio_beam.Beams(beams=beams))

    def _common_velocity_grid(self):
        """A linear velocity grid spanning the members' union range at their coarsest channel width.

        The coarsest channel width avoids up-sampling any member to a finer spectral resolution it
        never had (the spectral no-super-resolution guard, ADR-0003); the union range keeps every
        member's coverage (channels a member does not reach become NaN, zero-weighted in coadd)."""
        import numpy as np
        import astropy.units as u

        axes = [_member_velocity_axis(m.cube) for m in self._members]
        starts = [a.min() for a in axes]
        stops = [a.max() for a in axes]
        widths = [
            abs(a[1] - a[0]) if a.size > 1 else np.nan for a in axes
        ]  # km/s scalars
        step = float(np.nanmax(widths))
        lo, hi = float(np.min(starts)), float(np.max(stops))
        n = int(np.floor((hi - lo) / step)) + 1
        return (lo + step * np.arange(n)) * (u.km / u.s)

    # -- house-style grid plot (user story 12) ------------------------------------------
    def plot_grid(self, **kwargs):
        """Render every member in the house style — one panel per member; return ``(fig, axes)``.

        The browse-everything view (PRD #56 user story 12): a grid with one panel per member, each
        the member cube's canonical 2-D summary (its moment-0 map, exactly as ``cube.plot()`` shows
        a cube), titled with the member's identity (object / survey / species / transition) so a
        researcher visually compares a source across surveys and lines at a glance. It reuses the
        viz engine's map machinery and the shipped house style (DRY — the same look as
        :func:`~astrolyze.viz.plot_channel_maps`).

        Always works on a **heterogeneous** stack — there is no homogeneity precondition on
        browsing (that gate is #65's, on co-addition). An **empty** stack raises a clear error
        rather than producing an empty figure (nothing to plot is a caller mistake, not silent).

        Keyword arguments (``ncols``, ``cmap``, ``add_beam``, ``backend``) flow to the viz engine;
        see :func:`astrolyze.viz.plot_stack_grid`."""
        from astrolyze import viz

        return viz.plot_stack_grid(self, **kwargs)


@dataclass(frozen=True)
class HomogeneityReport:
    """The categorical homogeneity of a :class:`Stack`'s members (the #65 coadd-gate input).

    Lists the distinct values the members carry on each coadd-relevant categorical axis and names
    the axes with more than one value as :attr:`conflicts`. Read-only at stage 1 (#64): it
    *describes* whether a stack is co-addable; #65's ``coadd`` turns a non-empty :attr:`conflicts`
    into a hard error."""

    species: tuple
    transitions: tuple
    bunits: tuple
    #: The axes carrying more than one distinct value — empty iff the stack is homogeneous.
    conflicts: tuple


def _distinct_values(values) -> tuple:
    """Distinct non-null values, sorted — a stable summary of one categorical member axis."""
    return tuple(sorted({v for v in values if v is not None}))


# ======================================================================================
# Stage-2 co-addition machinery (#65): the homogeneity gate + noise-weighted combine
# ======================================================================================
# Kept beside the Stack it reduces. astrolyze stays thin: numpy does the weighted-mean and
# variance-propagation arithmetic; the value-add is the GATE (no meaningless average) and the
# correct inverse-variance physics (weight = 1/sigma^2, combined sigma = sqrt(1/sum w)). The noise
# is loaded from each member's noise companion via the existing io seam (NoiseModel companion).


class CoaddError(ValueError):
    """Raised when :meth:`Stack.coadd` is asked to combine members that are not co-addable.

    The no-silent-physics checkpoint of the stacking design (PRD #56 user story 16): a coadd is
    refused — never silently produced as a meaningless average — when the members are
    heterogeneous (mixed species / transition / bunit) or not yet on one spatial+spectral grid, or
    when the request itself is ill-formed (empty stack, unknown weights). A ``ValueError`` subclass
    so it reads as the programming/usage error it is; the message names the offending axis and the
    alignment step that fixes it."""


def _celestial_corners(cube):
    """The cube's four celestial-WCS corner sky positions as a ``(4, 2)`` array of (lon, lat) deg.

    A compact fingerprint of the spatial grid for the coadd gate: two cubes share a spatial grid
    iff they have the same pixel shape (checked separately) and map their corner pixels to the same
    sky. Read straight off the celestial WCS, so it catches a shifted/rotated/rescaled grid without
    comparing every pixel."""
    import numpy as np

    celestial = cube._sc.wcs.celestial
    ny, nx = cube.shape[1], cube.shape[2]
    pix_x = [0, nx - 1, 0, nx - 1]
    pix_y = [0, 0, ny - 1, ny - 1]
    world = celestial.pixel_to_world_values(pix_x, pix_y)
    return np.column_stack([np.asarray(world[0]), np.asarray(world[1])])


def _member_velocity_axis(cube):
    """The member's spectral axis as bare km/s values (the coadd spectral-grid fingerprint).

    Uses the cube's authoritative velocity axis where it carries the rest frequency + convention
    (the normal stacking case, after :meth:`Stack.shift_to_rest` / :meth:`Stack.to_velocity_grid`);
    falls back to the raw spectral-axis values so a frequency-only cube still compares by its own
    grid rather than raising here (the categorical gate already ran)."""
    import numpy as np
    import astropy.units as u

    try:
        return np.asarray(cube.velocity_axis().to_value(u.km / u.s), dtype="float64")
    except Exception:
        return np.asarray(cube.spectral_axis.value, dtype="float64")


def _coadd_members(stack: "Stack", *, weights: str):
    """Combine a gated (homogeneous, grid-aligned) stack into one :class:`~astrolyze.core.Cube`.

    Called only after :meth:`Stack._require_coaddable` has passed, so the members share one grid and
    one unit. Builds the combined data + propagated noise per the *weights* mode and hands them back
    through the first member's data seam, so the result carries the beam + WCS + context and gains
    combined provenance on its Metadata."""
    import numpy as np

    members = stack.members
    reference = members[0].cube
    unit = reference.unit

    data = [
        np.asarray(m.cube._data_quantity.to_value(unit), dtype="float64")
        for m in members
    ]
    data = np.stack(data, axis=0)  # (member, z, y, x)

    if weights == "noise":
        sigma = _member_sigmas(stack, reference_shape=reference.shape, unit=unit)
        combined, combined_sigma = _inverse_variance_combine(data, sigma)
    else:  # "uniform"
        sigma = _member_sigmas(
            stack, reference_shape=reference.shape, unit=unit, required=False
        )
        combined, combined_sigma = _uniform_combine(data, sigma)

    result = reference._with_data(combined * unit)
    result.metadata = _coadd_provenance(stack, weights=weights)
    # Attach the propagated combined noise as a companion product where one was produced.
    if combined_sigma is not None:
        result._coadd_sigma = combined_sigma * unit
    from astrolyze.core._base import _emit

    _emit(
        "coadd",
        params={"weights": weights, "n_members": len(members)},
    )
    return result


def _member_sigmas(stack, *, reference_shape, unit, required: bool = True):
    """Per-member per-voxel σ arrays (in *unit*), loaded from each member's noise companion.

    The inverse-variance weighting needs each member's σ; it comes from the noise companion stored
    alongside the member's source store (the L1 ``noise/`` subgroup), loaded through the existing
    :meth:`~astrolyze.core.NoiseModel.from_zarr_companion` seam and reconstructed as a full σ-cube.
    With ``required=True`` (the ``weights="noise"`` path) a member without a loadable companion
    raises :class:`CoaddError` — inverse-variance weighting cannot be faked (ADR-0003). With
    ``required=False`` (the ``weights="uniform"`` path) a missing companion yields ``None`` for that
    member and the combine simply omits its variance contribution.

    Each σ-cube is broadcast/validated to *reference_shape*; a σ companion whose grid does not match
    the (already grid-aligned) data is a real inconsistency and raises."""
    import numpy as np

    sigmas = []
    for member in stack.members:
        sigma = _load_member_sigma(member, unit=unit)
        if sigma is None:
            if required:
                raise CoaddError(
                    f"coadd(weights='noise') needs a noise companion for member "
                    f"{member._label()!r}, but none could be loaded from its origin store "
                    f"({member.origin_store_uri!r}). Inverse-variance weighting cannot be "
                    "fabricated (ADR-0003); use weights='uniform' or ensure each store carries "
                    "its noise/ companion."
                )
            sigmas.append(None)
            continue
        sigma = np.asarray(sigma, dtype="float64")
        if sigma.shape != reference_shape:
            raise CoaddError(
                f"coadd(): the noise companion of member {member._label()!r} has shape "
                f"{sigma.shape}, which does not match the data grid {reference_shape}. The σ "
                "field must be on the same grid as the data (align the noise model too)."
            )
        sigmas.append(sigma)
    return sigmas


def _load_member_sigma(member, *, unit):
    """The member's per-voxel σ as a bare array in *unit*, or ``None`` if no companion is loadable.

    Loads the noise companion from the member's origin store URI (the L1 ``noise/`` subgroup) via
    :meth:`~astrolyze.core.NoiseModel.from_zarr_companion`, converts it to the data unit through the
    noise model's unit hub, and reconstructs the full σ-cube. Returns ``None`` (rather than raising)
    when there is no origin store, no companion group, or the companion is UNRELIABLE — the caller
    decides whether a missing companion is fatal (the ``weights="noise"`` path) or merely skipped
    (``weights="uniform"``)."""
    import numpy as np

    uri = member.origin_store_uri
    if uri is None:
        return None
    try:
        from astrolyze.core import NoiseModel

        model = NoiseModel.from_zarr_companion(uri)
    except Exception:
        # No companion group / unreadable store: treat as "no σ here" (the caller gates on it).
        return None
    if model.unit != unit:
        model = model.to(unit)
    field = np.asarray(model.sigma_cube._data_quantity.to_value(unit), dtype="float64")
    if not np.isfinite(field).any():
        return None  # an UNRELIABLE (all-NaN) companion carries no usable σ.
    return field


def _inverse_variance_combine(data, sigma):
    """Inverse-variance weighted combine of stacked *data* with per-member σ; return (value, σ).

    The standard result (PRD #56 user story 15): per voxel ``w_i = 1/σ_i²``, the combined value is
    ``Σ w_i d_i / Σ w_i`` and the combined variance is ``1 / Σ w_i`` (so combined σ = ``√(1/Σ w_i)``).
    A voxel that is NaN in the data (off a member's coverage) or carries a non-finite/zero σ gets
    **zero weight** — it contributes nothing, never a fabricated value (ADR-0003). A voxel no member
    constrains (every weight zero) is left NaN in both the value and the σ (honest "not measured")."""
    import numpy as np

    data = np.asarray(data, dtype="float64")
    sigma_stack = np.stack([np.asarray(s, dtype="float64") for s in sigma], axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        weight = 1.0 / np.square(sigma_stack)
    # Zero-weight any voxel that is unusable: non-finite weight, non-finite data, or non-positive σ.
    usable = np.isfinite(weight) & np.isfinite(data) & (sigma_stack > 0.0)
    weight = np.where(usable, weight, 0.0)
    values = np.where(usable, data, 0.0)

    weight_sum = np.sum(weight, axis=0)
    weighted = np.sum(weight * values, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        combined = weighted / weight_sum
        combined_sigma = np.sqrt(1.0 / weight_sum)
    # A voxel no member constrains (Σw == 0) is genuinely unmeasured -> NaN, not 0/0 garbage.
    empty = weight_sum <= 0.0
    combined[empty] = np.nan
    combined_sigma[empty] = np.nan
    return combined, combined_sigma


def _uniform_combine(data, sigma):
    """Equal-weight combine of stacked *data*; return (value, σ-or-None).

    The combined value is the NaN-ignoring mean across members. When noise companions are available
    for every member, the combined σ of an N-member equal-weight mean is ``√(Σ σ_i²) / N`` (the
    uncorrelated error propagation of a plain average); when any member lacks a companion the σ is
    not propagated and ``None`` is returned (the value still combines)."""
    import numpy as np

    data = np.asarray(data, dtype="float64")
    combined = np.nanmean(data, axis=0)
    if any(s is None for s in sigma):
        return combined, None
    sigma_stack = np.stack([np.asarray(s, dtype="float64") for s in sigma], axis=0)
    n = sigma_stack.shape[0]
    with np.errstate(invalid="ignore"):
        combined_sigma = np.sqrt(np.nansum(np.square(sigma_stack), axis=0)) / n
    return combined, combined_sigma


def _coadd_provenance(stack, *, weights: str):
    """The combined-product :class:`~astrolyze.io.Metadata`: first member's context + lineage.

    The result inherits the (homogeneous) members' physical context from the first member and gains
    coadd lineage (PRD #56 user story 19): the contributing members' origin store URIs, the
    alignment steps applied, the weighting, and the corpus catalog version — attached on the
    reserved ``provenance`` seam so the header contract is unchanged."""
    from dataclasses import replace

    reference = stack.members[0].cube
    origins = [
        m.origin_store_uri for m in stack.members if m.origin_store_uri is not None
    ]
    provenance = {
        "operation": "coadd",
        "weights": weights,
        "n_members": len(stack.members),
        "member_origins": tuple(origins),
        "alignment_steps": tuple(stack._alignment),
        "catalog_version": stack.selection.catalog_version,
    }
    return replace(reference.metadata, provenance=provenance)


# -- the gather: Collection.stack delegates here (kept beside the container it builds) -------
def gather_stack(collection, position_or_sources, *, size, partial: str) -> Stack:
    """Gather cutouts into a :class:`Stack` — the implementation of :meth:`Collection.stack`.

    Two input shapes (PRD #56 user stories 10 & 11), dispatched on type:

    - **A single position** (a :class:`~astropy.coordinates.SkyCoord`): the covering cubes are
      found via :meth:`Collection.covering` and a :meth:`~astrolyze.core.Cube.cutout` of *size* is
      taken from each — a multi-survey, multi-line view of one target.
    - **A list of sources** (catalog object **names** and/or :class:`~astropy.coordinates.SkyCoord`
      s): each entry is resolved to position(s) and stacked, building a multi-target sample. A
      **name** is looked up in *this collection's catalog* (the names are catalog object names —
      there is no external SIMBAD/Sesame resolution, PRD #56 out-of-scope) and its catalog rows'
      footprint centre(s) are used as the stamp positions; an unknown name raises
      :class:`KeyError`. A bare :class:`~astropy.coordinates.SkyCoord` in the list is gathered via
      ``covering`` exactly like the single-position case.

    Each member carries its origin provenance (stamped by :meth:`Record.open`) and its identity
    (from the covering :class:`Record`). The stack records the :class:`Selection` provenance of the
    request (kind / targets / size / catalog version). A target that nothing covers contributes no
    member (an honest empty contribution), so a stack over an uncovered position is simply empty —
    browsing an empty stack is fine; only ``coadd`` (#65) will object to having nothing to add."""
    from astropy.coordinates import SkyCoord

    catalog_version = collection.catalog_version
    if isinstance(position_or_sources, SkyCoord) and position_or_sources.isscalar:
        members = _members_at_position(
            collection, position_or_sources, size=size, partial=partial
        )
        selection = Selection(
            kind="position",
            targets=position_or_sources,
            size=size,
            catalog_version=catalog_version,
        )
        return Stack(members, selection)

    # A sources list: names and/or SkyCoords. A non-scalar SkyCoord (an array of positions) is also
    # treated as a sources list of its individual coordinates.
    sources = _as_source_list(position_or_sources)
    members = []
    for source in sources:
        members.extend(
            _members_for_source(collection, source, size=size, partial=partial)
        )
    selection = Selection(
        kind="sources",
        targets=position_or_sources,
        size=size,
        catalog_version=catalog_version,
    )
    return Stack(members, selection)


def _as_source_list(position_or_sources) -> list:
    """Normalise the sources argument to a flat list of names / scalar SkyCoords.

    Accepts a Python list/tuple of names and/or SkyCoords, or a non-scalar SkyCoord array (each
    of whose positions becomes one entry). A bare string is treated as a single name (rather than
    iterated character-by-character — a classic foot-gun)."""
    from astropy.coordinates import SkyCoord

    if isinstance(position_or_sources, SkyCoord):
        # A non-scalar SkyCoord -> one scalar SkyCoord per position.
        return list(position_or_sources)
    if isinstance(position_or_sources, str):
        return [position_or_sources]
    return list(position_or_sources)


def _members_for_source(collection, source, *, size, partial) -> list:
    """The stack members for one source entry — a SkyCoord position or a catalog object name."""
    from astropy.coordinates import SkyCoord

    if isinstance(source, SkyCoord):
        return _members_at_position(collection, source, size=size, partial=partial)
    # A catalog object name: resolve it to its catalog rows' footprint centre(s) within THIS
    # collection (no external name resolution — PRD #56 out-of-scope). Each matching row's centre
    # is a stamp position; covering() then re-finds every cube (across surveys) at that centre.
    positions = _name_to_positions(collection, source)
    members = []
    for position in positions:
        members.extend(
            _members_at_position(collection, position, size=size, partial=partial)
        )
    return members


def _name_to_positions(collection, name: str) -> list["SkyCoord"]:
    """The footprint-centre position(s) of catalog object *name* within *collection*.

    The names ``stack()`` accepts are **catalog object names**, resolved against this collection's
    own catalog (PRD #56: no SIMBAD/Sesame in the library). Every catalog row whose ``object`` is
    *name* contributes its footprint centre (``ra_deg`` / ``dec_deg``) as a stamp position; the
    distinct centres are returned so a multi-survey source at one centre is not stamped twice. An
    **unknown name** raises a descriptive :class:`KeyError` (a typo must not silently yield an empty
    sample — ADR-0003), and a known name whose rows carry no centre raises too (it cannot be placed
    on the sky)."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    rows = [record.row for record in collection.records if record.object == name]
    if not rows:
        known = sorted({r.object for r in collection.records if r.object is not None})
        raise KeyError(
            f"no source {name!r} in this collection's catalog; known objects: "
            f"{', '.join(known) if known else '(none)'}. stack() resolves names against the "
            "catalog, not SIMBAD — pass a SkyCoord for an off-catalog position"
        )
    seen: set[tuple[float, float]] = set()
    positions = []
    for row in rows:
        if row.ra_deg is None or row.dec_deg is None:
            continue
        key = (round(row.ra_deg, 9), round(row.dec_deg, 9))
        if key in seen:
            continue
        seen.add(key)
        positions.append(SkyCoord(row.ra_deg * u.deg, row.dec_deg * u.deg))
    if not positions:
        raise KeyError(
            f"source {name!r} has no footprint centre in the catalog (ra_deg/dec_deg are null), "
            "so it cannot be placed on the sky to stack; pass a SkyCoord position instead"
        )
    return positions


def _members_at_position(collection, position, *, size, partial) -> list[StackMember]:
    """A stack member for every cube covering *position* — a cutout of *size* from each.

    The core gather: :meth:`Collection.covering` finds every cube whose WCS contains *position*
    (across surveys), and each is opened and :meth:`~astrolyze.core.Cube.cutout` of the requested
    angular *size* taken around *position*. The cutout carries its origin provenance (from
    :meth:`Record.open`) and the member carries its identity (from the covering :class:`Record`).

    A covering cube whose footprint cannot fit the full stamp around *position* is handled by the
    cutout's own ``partial`` rule (default ``"raise"`` — no silent clipping); pass
    ``partial="trim"`` through :meth:`Collection.stack` to opt into the trimmed stamp instead."""
    members = []
    for record in collection.covering(position).records:
        cube = record.open()
        stamp = cube.cutout(position, size, partial=partial)
        members.append(StackMember.from_record(record, stamp))
    return members


__all__ = [
    "Stack",
    "StackMember",
    "Selection",
    "HomogeneityReport",
    "CoaddError",
    "gather_stack",
]
