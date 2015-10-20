===========================
 Installation of astrolyze
===========================


astrolyze is only tested on Linux/Ubuntu so far.

Dependencies
============


Python
------

astrolyze depends on the following python packages

::

  numpy
  pyfits
  matplotlib
  scipy
  pywcs
  pysqlite2
  generaltools
  
In Ubuntu (and thus most probable also Debian) these dependencies can be
installed via::

  sudo apt-get install python-matplotlib python-pywcs python-scipy python-numpy
  python-pysqlite2 python-pyfits

`generaltools` can be installed via pypi::

  sudo pip install generaltools
  
.. warning:: 

   It does not work with python3 yet. Since the Gildas Python extension is not yet available for python 3!!

Gildas
------

To be able to use GILDAS from within astrolyze it is enough to have a working
GILDAS installation compiled with the python support. The actual GILDAS version
and installation instructions can be found here:

http://www.iram.fr/IRAMFR/GILDAS/


Miriad
------

Also Miriad just has to be installed and working. At the moment only the smooth
function of miriad is used in astrolyze but it is worth installing it since it
it a common task. Miriad and instructions for installation can be found here:


ftp://ftp.atnf.csiro.au/pub/software/miriad/INSTALL.html


Installation
============

Astrolyze is available via pypi, use::

  sudo pip install astrolyze --user

to install the package. 

In the current state the setup does not allow a custom installation path. The
files will be installed depending on your system configuration in either of the
following locations::

 /usr/lib/pythonX.Y/site-packages
 /usr/local/lib/pythonX.Y/site-packages

Here X and Y are the major and minor number your python installations.


Configuration of the (optional) Database
========================================

The setup.py script will generate a ``parameter.db`` sqlite database containing
information about the maps and files opened in Astrolyze. From the keyword of
the ''Naming Convention''.

In the current state  astrolyze reads in

 - Additional informations of the source
 - Frequencies and wavelengths
 - Calibration error for specific telescopes

The database has to be populated by the user. This is done via the three text
files in the ``/home/USER/.astrolyze/cfg`` folder that contains:

galaxy_parameter.txt:

.. literalinclude:: ../../cfg/galaxy_parameter.txt

line_parameter.txt:

.. literalinclude:: ../../cfg/line_parameter.txt

calibration_error.txt:

.. literalinclude:: ../../cfg/calibration_error.txt

The Names of the source telescope and lines have to be exactly how they are
used in the map names. However the writing can be an arbitrary mix of upper an
lower case characters. Internally Astrolyze converts them to upper case before
comparing.

To generate the database from these files you have to run the script::

  > ~/.local/bin/astrolyze_opt_db_setup.py

installed in `.local/bin`.

