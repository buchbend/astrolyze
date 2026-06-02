# I/O & provenance: FITS header is authoritative, filename is a synced browsable projection

astrolyze has carried metadata two ways: v1/astrolyze3 encoded it in the **filename**
(`source_telescope_species_fluxunit_resolution_...`, parsed on load, sqlite-backed);
astrolyze2 pivoted to **FITS header keywords**. Each has a virtue — headers are
self-describing and robust; filenames are great for quick visual browsing of a data
directory. We keep both, with a clear authority order.

**Decision — (a) header-keyword convention, with the filename kept in sync:**

- **The FITS header is the single source of truth.** astrolyze defines a metadata schema in
  the header: rest frequency + **velocity convention** (mandatory per ADR-0003), beam
  (BMAJ/BMIN), distance, species, calibration error, object, name-tag, and provenance.
  **The header is what we read** on load. Pure FITS, no required sidecar.
- **The filename is a derived, human-browsable projection of the header** (astrolyze-1
  style, for quick directory browsing). It is generated/updated **from the header on
  ingest**, and on every **header write we check filename↔header consistency** and update
  the filename to match. The filename never overrides the header; it mirrors it.

**Loader strictness — (ii) lazy enforcement.** When a loaded FITS lacks mandatory context
(rest frequency / velocity convention — e.g. real archival PHANGS/SEDIGISM headers), the
**load succeeds** but the object is "incomplete": any operation that needs the missing
context (unit conversion, velocity work) **raises a clear error** telling the user to set
it. Never silently guess (ADR-0003), never refuse to open a real-world file.

Rejected for v1: (b) header + sidecar metadata file (two files to sync), and (c) a full
in-memory `Provenance` object serialized to FITS HISTORY — valuable for downstream
heterogeneity tracking, but YAGNI now. Leave a `Provenance` seam so (c)
can layer on without breaking the header contract.

## Consequences

- One self-describing file; no header/sidecar sync problem. Filename stays informative and
  consistent because it is regenerated, not hand-maintained.
- **Ingest renames files** (surprising if unexpected): bringing an external cube into
  astrolyze derives the canonical name from its header. Document this clearly.
- The exact header keyword schema and the filename field grammar are an **open detail** to
  pin during build (harvest v1's field set; reconcile with astrolyze2's keyword names).
- Round-trip (load → save → load) must preserve all schema fields — a test obligation.
