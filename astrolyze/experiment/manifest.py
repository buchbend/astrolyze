"""The DB-backed dataset manifest (ADR-0009).

The manifest is the generated registry of what lives in an experiment: one row per dataset,
carrying its identity (where it came from, a citable DOI) and the *full physical provenance*
from the :class:`~astrolyze.io.schema.Metadata` schema. Ingest (#12) generates and keeps it in
sync — it is never hand-edited — and it is shaped to be UI-readable so a future frontend can
render it without rework.

Two design points worth stating:

- **Provenance is the io schema, not a parallel field set.** Columns map 1:1 to ``Metadata``;
  we persist *primitives* (floats + unit strings) and reconstruct the astropy/radio_beam
  objects on read, so the registry round-trips the *value* (a beam, a rest frequency, a
  velocity convention) without coupling the table to those library types. The reconstructed
  :class:`~astrolyze.io.schema.Metadata` rides back on each :class:`DatasetRecord`.
- **Idempotent on source path; backend swappable.** ``register`` keys on the dataset's source
  path (relative to the experiment — the caller resolves that): re-registering updates the
  existing row in place, never duplicating, so re-ingest is a normal step and the row id stays
  stable. The interface is plain SQLAlchemy over a ``db_url`` with no sqlite-only SQL leaking
  out; sqlite is merely the default backend, swappable as scale grows (ADR-0009).
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING

import astropy.units as u
import radio_beam
from sqlalchemy import Float, Integer, String, create_engine, select
from sqlalchemy.engine import make_url
from sqlalchemy.orm import DeclarativeBase, Mapped, Session, mapped_column
from sqlalchemy.pool import StaticPool

from astrolyze.io.schema import SCHEMA_VERSION, Metadata
from astrolyze.units import VelocityConvention

if TYPE_CHECKING:
    from .layout import Experiment

# Fallback when the experiment config declares no manifest db_url. Mirrors the [manifest]
# db_url written by layout.DEFAULT_CONFIG; kept relative so it resolves against the experiment
# root (see for_experiment / _resolve_db_url), never the current working directory.
DEFAULT_DB_URL = "sqlite:///manifest.db"


# -- ORM model -------------------------------------------------------------------------
class _Base(DeclarativeBase):
    pass


class _DatasetRow(_Base):
    """One registered dataset. Provenance columns mirror :class:`Metadata` (primitives +
    unit strings); the astropy/radio_beam objects are reconstructed on read."""

    __tablename__ = "datasets"

    # Identity. ``source_path`` is the idempotency key (unique); the integer primary key is
    # preserved across re-registration (update-in-place), so a dataset's id is stable.
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    source_path: Mapped[str] = mapped_column(
        String, nullable=False, unique=True, index=True
    )
    doi: Mapped[str | None] = mapped_column(String, nullable=True)
    ingested_at: Mapped[datetime] = mapped_column(nullable=False)
    schema_version: Mapped[int] = mapped_column(Integer, nullable=False)

    # Provenance (io Metadata, flattened to primitives).
    object: Mapped[str | None] = mapped_column(String, nullable=True)
    telescope: Mapped[str | None] = mapped_column(String, nullable=True)
    species: Mapped[str | None] = mapped_column(String, nullable=True)
    rest_frequency_hz: Mapped[float | None] = mapped_column(Float, nullable=True)
    velocity_convention: Mapped[str | None] = mapped_column(String, nullable=True)
    bmaj_deg: Mapped[float | None] = mapped_column(Float, nullable=True)
    bmin_deg: Mapped[float | None] = mapped_column(Float, nullable=True)
    bpa_deg: Mapped[float | None] = mapped_column(Float, nullable=True)
    bunit: Mapped[str | None] = mapped_column(String, nullable=True)
    distance_value: Mapped[float | None] = mapped_column(Float, nullable=True)
    distance_unit: Mapped[str | None] = mapped_column(String, nullable=True)
    calibration_error: Mapped[float | None] = mapped_column(Float, nullable=True)
    name_tag: Mapped[str | None] = mapped_column(String, nullable=True)


# Public filter name -> ORM column name. Querying is restricted to these scalar columns so an
# unknown filter is a clear error rather than a silent empty result.
_FILTERABLE = {
    "source_path": "source_path",
    "doi": "doi",
    "object": "object",
    "telescope": "telescope",
    "species": "species",
    "velocity_convention": "velocity_convention",
    "bunit": "bunit",
    "name_tag": "name_tag",
    "schema_version": "schema_version",
}


# -- public read value object ----------------------------------------------------------
@dataclass(frozen=True)
class DatasetRecord:
    """A registered dataset as a caller sees it: identity plus the reconstructed
    :class:`~astrolyze.io.schema.Metadata` provenance. Decoupled from the ORM row so the
    storage backend stays an implementation detail."""

    id: int
    source_path: str
    doi: str | None
    ingested_at: datetime
    schema_version: int
    metadata: Metadata


class Manifest:
    """The dataset registry for one experiment, backed by ``db_url``.

    Constructing it opens (creating if needed) the registry. The interface is backend-generic;
    sqlite is the default (``sqlite:///manifest.db`` in the experiment config), swappable.
    """

    def __init__(self, db_url: str) -> None:
        self._engine = create_engine(db_url, **_engine_kwargs(db_url))
        _Base.metadata.create_all(self._engine)

    @classmethod
    def for_experiment(
        cls, experiment: Experiment, db_url: str | None = None
    ) -> "Manifest":
        """Open the manifest for *experiment*.

        Reads ``manifest.db_url`` from the experiment config (or uses an explicit *db_url*),
        then resolves a **relative sqlite path against the experiment root** — so the registry
        always lives inside the experiment, never wherever the process happens to be running
        (the config default is relative on purpose). Absolute paths, ``:memory:``, and
        non-sqlite backends pass through untouched."""
        if db_url is None:
            db_url = experiment.settings.get("manifest.db_url") or DEFAULT_DB_URL
        return cls(_resolve_db_url(str(db_url), experiment.root))

    def register(
        self, metadata: Metadata, source_path: Path | str, doi: str | None = None
    ) -> DatasetRecord:
        """Register (or re-register) the dataset at *source_path* with its *metadata* and an
        optional *doi*. Idempotent on *source_path*: an existing row is updated in place
        (stable id), a new one inserted otherwise. The latest call's values win."""
        key = _normalise_path(source_path)
        with Session(self._engine) as session:
            row = session.scalar(
                select(_DatasetRow).where(_DatasetRow.source_path == key)
            )
            if row is None:
                row = _DatasetRow(source_path=key)
                session.add(row)
            for column, value in _provenance_columns(metadata).items():
                setattr(row, column, value)
            row.doi = doi
            row.schema_version = SCHEMA_VERSION
            row.ingested_at = _now()
            session.commit()
            session.refresh(row)
            return _row_to_record(row)

    def get(self, id: int) -> DatasetRecord | None:
        """Return the dataset with this id, or ``None`` if it is not registered."""
        with Session(self._engine) as session:
            row = session.get(_DatasetRow, id)
            return _row_to_record(row) if row is not None else None

    def query(self, **filters) -> list[DatasetRecord]:
        """Return datasets matching every supplied equality *filter* (e.g.
        ``query(object="NGC0628")`` or ``query(species="CO21")``). Filtering on an unknown
        column raises :class:`ValueError`."""
        statement = select(_DatasetRow)
        for name, value in filters.items():
            column = _FILTERABLE.get(name)
            if column is None:
                raise ValueError(
                    f"unknown manifest filter: {name!r} "
                    f"(filter on one of {sorted(_FILTERABLE)})"
                )
            statement = statement.where(
                getattr(_DatasetRow, column) == _coerce_filter(name, value)
            )
        with Session(self._engine) as session:
            return [_row_to_record(r) for r in session.scalars(statement).all()]

    def all(self) -> list[DatasetRecord]:
        """Return every registered dataset."""
        with Session(self._engine) as session:
            return [
                _row_to_record(r) for r in session.scalars(select(_DatasetRow)).all()
            ]


# -- engine configuration --------------------------------------------------------------
def _engine_kwargs(db_url: str) -> dict:
    """Engine tuning for the sqlite in-memory case only.

    An ``sqlite://`` (``:memory:``) database lives inside a single DBAPI connection — with the
    default pool each new connection would see an empty schema. ``StaticPool`` keeps one shared
    connection so the tables and rows persist across ``register``/``get``/``query`` on a
    :class:`Manifest`. This is an internal adaptation to the declared backend; the public
    interface stays backend-generic."""
    url = make_url(db_url)
    if url.get_backend_name() == "sqlite" and url.database in (None, "", ":memory:"):
        return {"poolclass": StaticPool, "connect_args": {"check_same_thread": False}}
    return {}


def _resolve_db_url(db_url: str, root: Path) -> str:
    """Anchor a relative sqlite file path in *db_url* to the experiment *root*.

    The config default (``sqlite:///manifest.db``) is relative; left as-is it would resolve
    against the current working directory. We re-anchor it to *root* so the manifest lands in
    the experiment. In-memory (``:memory:``), already-absolute paths, and non-sqlite backends
    are returned unchanged."""
    url = make_url(db_url)
    database = url.database
    if url.get_backend_name() != "sqlite" or not database or database == ":memory:":
        return db_url
    if Path(database).is_absolute():
        return db_url
    absolute = (Path(root) / database).resolve()
    return url.set(database=str(absolute)).render_as_string(hide_password=False)


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _normalise_path(source_path: Path | str) -> str:
    """The idempotency key: a forward-slash string so separators don't fork duplicate rows."""
    return Path(source_path).as_posix()


# -- Metadata <-> primitive columns ----------------------------------------------------
def _provenance_columns(metadata: Metadata) -> dict:
    """Flatten the io :class:`Metadata` provenance to storable primitives (None stays None)."""
    beam = metadata.beam
    distance = metadata.distance
    return {
        "object": metadata.object,
        "telescope": metadata.telescope,
        "species": metadata.species,
        "rest_frequency_hz": (
            metadata.rest_frequency.to_value(u.Hz)
            if metadata.rest_frequency is not None
            else None
        ),
        "velocity_convention": (
            VelocityConvention(metadata.velocity_convention).value
            if metadata.velocity_convention is not None
            else None
        ),
        "bmaj_deg": beam.major.to_value(u.deg) if beam is not None else None,
        "bmin_deg": beam.minor.to_value(u.deg) if beam is not None else None,
        "bpa_deg": beam.pa.to_value(u.deg) if beam is not None else None,
        "bunit": str(metadata.bunit) if metadata.bunit is not None else None,
        "distance_value": distance.value if distance is not None else None,
        "distance_unit": str(distance.unit) if distance is not None else None,
        "calibration_error": (
            float(metadata.calibration_error)
            if metadata.calibration_error is not None
            else None
        ),
        "name_tag": metadata.name_tag,
    }


def _row_to_metadata(row: _DatasetRow) -> Metadata:
    """Reconstruct the io :class:`Metadata` (astropy/radio_beam types) from stored primitives."""
    beam = None
    if row.bmaj_deg is not None:
        beam = radio_beam.Beam(
            major=row.bmaj_deg * u.deg,
            minor=row.bmin_deg * u.deg,
            pa=row.bpa_deg * u.deg,
        )
    distance = None
    if row.distance_value is not None:
        distance = row.distance_value * u.Unit(row.distance_unit)
    return Metadata(
        object=row.object,
        telescope=row.telescope,
        species=row.species,
        rest_frequency=(
            row.rest_frequency_hz * u.Hz if row.rest_frequency_hz is not None else None
        ),
        velocity_convention=(
            VelocityConvention(row.velocity_convention)
            if row.velocity_convention is not None
            else None
        ),
        beam=beam,
        bunit=u.Unit(row.bunit) if row.bunit is not None else None,
        distance=distance,
        calibration_error=row.calibration_error,
        name_tag=row.name_tag,
    )


def _row_to_record(row: _DatasetRow) -> DatasetRecord:
    return DatasetRecord(
        id=row.id,
        source_path=row.source_path,
        doi=row.doi,
        ingested_at=row.ingested_at,
        schema_version=row.schema_version,
        metadata=_row_to_metadata(row),
    )


def _coerce_filter(name: str, value):
    """Normalise a filter value to its stored form (enum -> ``.value``, path -> posix str)."""
    if name == "velocity_convention" and value is not None:
        return VelocityConvention(value).value
    if name == "source_path" and value is not None:
        return _normalise_path(value)
    return value


__all__ = ["Manifest", "DatasetRecord"]
