r"""
The ``maps`` module aims to make handling of astronomical maps (mainly in Radio
and Infrared astronomy) easier. So far it supports maps in the FITS, GILDAS and
MIRIAD data formats. It wraps functions of GILDAS and MIRIAD and provides own
new functions written with python on the basis of pyfits. 

It is using a 'Name Convention' for ease of use. Meaning that the file name
already includes basic information about the map it contains. A Name that
follows this 'Convention' is eg: M33_30m-HERA_CO21_Ta*_12_cube.fits

All items **MUST** be seperated by an undescore (_) and **HAVE** to include at
minimum the following properties:: 

1. source
2. telescope
3. wavelength OR frequency OR lineName
4. flux unit
5. resolution

These items are transferred to variables of the Map class that is used to
handle the maps. See below or in the tutorial for more explanation. 
Additionaly the Map Class recognizes all following items as comments. In the
name example above cube would be a comment. Comments are not transfered to
individual internal variables of the map objects but are passed on as a list to
the single variable comments.

The last item is followed by the files extension:

* .fits -> FITS 
* .gdf, .mean, .velo, .width, .lmv -> GILDAS
* nothing -> MIRIAD (Miriads file format uses directories to store the data.)

Maps that are not following this name convention are not supported. This is
to assure that all parts of the program work, since they strongly depend on the
parameters passed on by the name, as is explained below or in the tutorial.

Altough this is somewhat redundant to the header information of the files, it 
has been decided to go that way since unfortunately not all headers are kept up
to date and manipulating the file name is easier to do. 

The maps module tries to keep track if a variable that should go into the
header of a fits file is changed and up-dates the header subsequently (Maybe
not true in all cases, though.).

Using this name convention has another benefit since it makes the life of your
fellow astronomers easier when they have to work with your data since they
readily know their most important basic properties.
"""
#__all__ = ["mapClassMain", "mapClassFits", "mapClassGildas", "mapClassMiriad"]
# TODO: Decide on how astrolyze is initated best.
from main import *
from fits import *
from gildas import *
from miriad import *
from stack import *
import tools as mtools
