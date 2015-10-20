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

   Astrolyze does not work in python3 yet. Since the Gildas Python extension is not yet available for python 3!!

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

When installed with the --user flag the package will be installed in::

    /home/USERNAME/.local/lib/python2.7/site-packages/astrolyze

The configuration files for `astrolyze` are installed in::

    /home/USERNAME/.astrolyze/

Further a script to configure the database for additional information is installed in::

    /home/USERNAME/.local/bin/


Configuration of the (optional) Database
========================================

Astrolyze can store META Information for Sources, Lines and Calibration
Uncertainties in a database. Based on the information deduced from the naming
convention the additional information is loaded automatically when a map is
opened.

`astrolyze` can read:

 - Additional informations of the source
 - Frequencies and wavelengths
 - Calibration error for specific telescopes

The database has to be populated after installation This is done via the three text
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

