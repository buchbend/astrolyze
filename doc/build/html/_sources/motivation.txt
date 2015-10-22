.. _label_motivation:

========================================
Motivation - What is astrolyze all about
========================================

First, before delving into the details of the astrolyze package, here a few
examples are given of what astrolyze is all about. The following is a short
teaser of some of the features of astrolyze; a thorough introduction with more
features and possibilities are given in the :ref:`label-naming-convention`.

Inter-operating Fits, Gildas and Miriad
---------------------------------------

A central concept of `astrolyze` is to enable to use different astronomically map formats and
different programs from within python seamlessly. Say, for example, that you
have a Fits image called ``M33_30m_12CO10_Tmb_12.fits`` and that you want to
smooth it in miriad to 40 arcsec resolution, re-project it with Gildas to a new
central coordinate and finally convert it back to fits-format:

This is how you would do it with astrolyze::

  from astrolyze import *
  map_ = FitsMap('M33_30m_12CO10_Tmb_12.fits')
  map_ = map_.toMiriad()
  map_ = map_.smooth(40)
  map_ = map_.toGildas()
  map_ = map_.reproject(coordinate=['01:34:50.890', '+31:08:28.03'])
  map_ = map_.toFits()

Please note the special format of the map name. The format follows a certain
``naming convention`` that **has** to be used with astrolyze. The reason for the
naming conventions and it's internal logic is explained below.

.. note::
   
   As a side note I use an underscore for the ``map_`` variable, because
   otherwise the python function ``map`` is overwritten which may lead to
   problems.

Changing map units
------------------

With astrolyze the units of a map can be quickly transformed between common
units used in (Radio-) Astronomy (as far as the conversion was implemented
already). Take for example again the ``M33_30m_12CO10_Tmb_12.fits`` map.
Following the naming convention this map is in main beam temperature
(Tmb). Changing its units to ``JyB`` is as easy as::

  from astrolyze import *
  map_ = FitsMap('M33_30m_12CO10_Tmb_12.fits')
  map_ = map_.change_unit('JyB')

Another side note: It's only that easy when astrolyze is set-up with the
database so that astrolyze knows about the frequency/wavelenght of the line in
the map.

.. warning::
   Please double check the results of the unit conversion.  

Working with stacks of images
-----------------------------

The previous examples demonstrated how single maps can be used in astrolyze.
It is further possible to work on a stack of images and perform tasks on all of
them.

To create a stack, all files that are going to be in the stack have to be
located in one folder (with possible sub-folders NOTE: The functionality with
sub-folders is not thoroughly tested, yet.)
A stack is initialized as follows::

  from astrolyze import *
  example_stack = Stack('path_to_folder')

The maps can be a mix of GILDAS, Fits and MIRIAD maps. The Instance of the
Stack object (here: ``example_stack``) contains a variable called stack which
is a list with Instances of the corresponding maps Objects (GildasMap, FitsMap
and MiriadMap).

The stack module provides several tools to ``unify`` the stack for further
analysis. The maps can be all re-gridded and re-projected to the same central
coordinates, pixel-sizes and dimensions as a given template image via::

  example_stack.unify_dimensions(template='path_to_template_file',
                                 folder='path_to_output_folder') .

Also the maps can be all smoothed to the same resolution, by default this is
the largest resolution found in the stack but can also be given manually::

  example_stack.unify_resolutions(folder='path_to_output_folder') .

Astrolyze includes unit conversions that can be used to change all maps to the
same resolutions as long as the input and output units are programmed. For
the stack::

  example_stack.unify_units(folder='path_to_output_folder') ,

can be used.

Producing SEDs
--------------

Build on-top of the stack module astrolyze, contains the ``sed`` module
which allows to analyze and plot **dust-seds**. SEDs can be read out for an
arbitrary number of positions or for the entire maps. In the latter case
temperature, mass, beta and $chi^2$ maps will be created.

When all maps have the same resolution and dimension (i.e. pixel size and
number) producing temperature maps can be done as follows::

  from astrolyze import *
  sed = SedStack(folder='path_to_input_folder', full_map=True,
         output_folder='path_to_output_folder') .

To generate SEDs at given coordinates it is easiest to provide a separate file
(e.g. ``coordinates.txt``) with the names and coordinates of the positions to be
extracted as follows::

  source_1    1:34:7.00     +30:47:52.00
  source_2    1:33:55.80    +30:43:2.00
  source_3    1:33:52.40    +30:39:18.00
    .              .             .
    .              .             .
    .              .             .

Then a stack of seds can be created::

  from astrolyze import *
  seds = SedStack(folder='path_to_input_folder', flux_acquisition='pixel')

By default the SED is also directly fitted. One can produce a quick preview
plot of the SEDs via::

  for i in seds.sed_stack:
      i.create_figure()


Not only images ...
-------------------

Last but not least astrolyze is also able to work with 30m class spectra from
within python based on the same principles used to work with images/maps. The
implementation makes extensive use of pyGildas.  For example if you have a
file with the spectra of a cube, e.g. ``M33_30m_12CO10_Tmb_21_cube.30m`` you
can extract the spectra that corresponds closest to a given coordinate as
follows::

  from astrolyze import *
  spectra = ClassSpectra('M33_30m_12CO10_Tmb_21_cube.30m')
  coordinate = ['1:34:7.00', '+30:47:52.00']
  spectrum = spectra.get_spectra_from_cube(coordinate)
  # To show the spectrum in the Class window
  spectrum.quick_view()


