'''
The Map class is the parent class for the handling of astronomical map files 
in the FITS, GILDAS and MIRIAD data formats.
It is meant to make handling maps easier by wrapping functions of GILDAS and 
MIRIAD and providing new functions in python on the basis of pyfits.
It is using a 'Name Convention' for ease of use. Meaning that the file name 
already includes basic information about the map it contains. A Name that 
follows this 'Convention' is eg:

M33_30m-HERA_CO21_Ta*_12_cube.fits

All items **MUST** be seperated by an undescore (_) and have to include at 
minimum the following properties::

1. source
2. telescope
3. wavelength OR frequency OR lineName
4. flux unit
5. resolution

Additionaly the map Class recognizes all following items as comments. In the 
name example above the comment be "cube". Comments are not transfered to 
internal variables of the map objects.

The last item is followed by the files extension:

* .fits -> FITS 
* .gdf -> GILDAS
* nothing -> MIRIAD (Miriads file format uses directories to store the data)

Maps that are not following this name convention are not supported 
to assure that all parts of the program work, since they mostly 
depend on the items set as will be explained below.

Also it makes the life of your fellow astronomers easier when they have to work with your data since they directly know their basic properties.
'''
#__all__ = ["mapClassMain", "mapClassFits", "mapClassGildas", "mapClassMiriad"]

from main import *
from fits import *
from gildas import *
from miriad import *
