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

The Stack carries **selection provenance** (:class:`Selection`): what was asked of
:meth:`~astrolyze.collection.Collection.stack` (the position or the source list, the stamp size,
the catalog version), accessible as :attr:`Stack.selection`. :meth:`~Stack.filter` narrows a stack
to a homogeneous subset (``species=`` / ``survey=`` / ``transition=`` / ``object=``) — the path
from a browse-everything stack to a co-addable one (user story 13) — and preserves that provenance.

Seams left for the alignment + co-addition stage (#65):

- **Homogeneity is read-only here.** :attr:`Stack.is_homogeneous` /
  :meth:`Stack.homogeneity_report` describe whether the members share a species / transition /
  bunit (the notion #65's ``coadd`` will *gate* on), but nothing in this stage *requires* it — a
  heterogeneous stack is fully usable. #65 turns the same report into a hard precondition.
- **Per-member broadcast.** :meth:`Stack.map` applies a per-member ``Cube -> Cube`` operation
  across the stack and returns a new Stack; #65's ``to_common_beam`` / ``to_velocity_grid`` /
  ``shift_to_rest`` are exactly such per-member operations, so they slot onto this seam.
- **Co-addition hooks onto the filtered, homogeneous subset.** The intended path is
  ``stack.filter(species=...).coadd(...)``; #65 adds ``coadd`` (and the alignment methods) onto
  this class, reusing :attr:`members` and the homogeneity report. No ``coadd`` /
  ``to_common_beam`` / ``to_velocity_grid`` / ``shift_to_rest`` exists yet — they are #65.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Callable

if TYPE_CHECKING:
    from astropy.coordinates import SkyCoord

    from astrolyze.core import Cube

    from ._facade import Record


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


__all__ = ["Stack", "StackMember", "Selection", "HomogeneityReport", "gather_stack"]
