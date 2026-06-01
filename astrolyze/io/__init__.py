"""I/O and the metadata schema.

The FITS header is authoritative; the filename is a derived, browsable projection of it.
Loading is lazy: an incomplete header still opens, but operations needing the missing
context raise. Implemented in issue #0003.
"""
