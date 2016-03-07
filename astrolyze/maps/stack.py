# Copyright (C) 2012, Christof Buchbender
# BSD License (License.txt)

# standart library
import math
import os
import string
import sys
import matplotlib.pyplot as plt
import numpy as np

# from generaltools
from generaltools import log_tools

# from astrolyze
import main
import fits
import gildas
import miriad

import astrolyze.functions.constants as const

USER = os.getenv("USER")

class Stack(object):
    r""" Class to work with several maps at once.

    Allows to treat all images inside a folder simultaneously. This allows to
    to perform the same transformations on all maps.

    The images in the input folder can have arbitrary formats, units,
    resolutions and other properties. The Stack class can help to unify the
    different parameters of the maps to help comparing them.


    Parameters
    ----------

    folder : str
        The folder containing the images that are to be treated as a stack
    data_format : str
        If not None filter images by format.

    Examples
    --------

    The stack class can be used like builtin objects that are iterable e.g.:

        >>> st = Stack("folder_name")
        >>> list_map_names = [i.map_name for i in st]


    .. note::

      `Stack` is also the basis for the :py:mod:`astrolyze.sed.Sed` package.
       """
    def __init__(self, folder, data_format=None):
        r"""
        Initialization of the stack. Reading in the maps and creating a list
        of ``Map`` class objects, depending on the input format.
        """
        self.log = log_tools.init_logger(
            directory="/home/{}/.astrolyze/".format(USER),
            name="astrolyze"
        )
        self.gildas_formats = ['gdf', 'mean', 'velo', 'width', 'lmv',
                               'lmv-clean']
        self.fits_formats = ['fits']
        self.miriad_formats = ['']
        self.folder = folder
        self.get_list(data_format=data_format)
        print self.list
        self.stack = []
        for i in self.list:
            self.stack += [self.get_map_format(i)]
        self.resolutions = []
        self.units = []
        self.dataFormats = []
        for i in self.stack:
            print i
            self.resolutions += [i.resolution]
            self.units += [i.fluxUnit]
            self.dataFormats = [i.dataFormat]
        self.log = log_tools.init_logger(
            directory="/home/{}/.astrolyze/".format(USER),
            name="astrolyze"
        )

    def __repr__(self):
        string = "Stack(\"{}\")\nContains:\n".format(self.folder)
        for map_ in self.stack:
            string += "{}\n".format(map_)
        return string

    def __len__(self):
        return len(self.stack)

    def __getitem__(self, position):
        return self.stack[position]

    def __iter__(self):
        return (i for i in self.stack)

    def get_map_format(self, map_name):
        r""" This function creates returns the correct ``GildasMap``,
        ``FitsMap`` or ``MiriadMap`` object without that the data format of the
        map has to be given.

        Parameters
        ----------

        map_name : string
           Path and name of the input map.

        Returns
        -------

        Either a ``GildasMap``, ``FitsMap`` or ``MiriadMap` Either a
        ``Map`` class object depending on the data format of the input map.
        """
        map_name_split = map_name.strip().split('.')
        dataFormat = map_name_split.pop()
        if dataFormat in self.gildas_formats:
            return gildas.GildasMap(map_name)
        if dataFormat in self.fits_formats:
            return fits.FitsMap(map_name)
           # check for miriad maybe not perfect.
        if len(map_name.split('_')) >= 5 and os.path.isdir(map_name):
            return miriad.MiriadMap(map_name)

    def get_list(self, data_format=None, depth=False):
        r"""
        Loading a list of files in all sub-folders.

        Parameters
        ----------

        folder : string
            The path to the folder that has to be parsed.
        data_format : string
            Search for specific files containing the string, e.g.
            '.fits'
        depth : integer
            The steps into the sub-folder system. Defaults to maximum depth.

        Returns
        -------

        final_list : array
            Array with the string to the files in the folder and sub folders.
        folder_list :
            Array with the strings to the folders. Only if depth is set.
        """

        def get_sub_folder(folder):
            r"""
            Returns a list with the contents of a folder.
            """
            miriad_subfolder = ['header', 'mask', 'image', 'history']
            os.system('ls ' + str(folder) + ' > temp')
            filein = open('temp').readlines()
            list = []
            for i in filein:
                # We have to exclude sub_folders that are miriad file
                # sub-folders.
                if (i.strip() in miriad_subfolder
                    and len(folder.split('_')) >= 5):
                    pass
                else:
                    list += [str(folder) + '/' + str(i.strip())]
            os.system('rm -r temp')
            return list
        self.folder_list = [self.folder]
        self.list = []
        # Variable to steer the depth of the search.
        x = 1
        execute = True
        while execute:
            # Variables to check if the algorithm found folder
            # and or files in the last iterations. If only files are
            # found the loop is stopped.
            folder_check = False
            file_check = False
            new_folder_list = []
            for folder in self.folder_list:
                list = get_sub_folder(folder)
                for sub_folder in list:
                    if os.path.isdir(sub_folder):
                        new_folder_list += [sub_folder]
                        folder_check = True
                        if len(sub_folder.split('_')) >= 5:
                            self.list += [sub_folder]
                    if os.path.isfile(sub_folder):
                        if (data_format is not None and data_format
                            in sub_folder):
                            if sub_folder not in self.list:
                                self.list += [sub_folder]
                        if (data_format is None):
                                self.list += [sub_folder]
                        file_check = True
            if not folder_check and file_check:
                execute = False
            x = x + 1
            if x > depth and depth:
                execute = False
            if x > 1000:
                print 'Not all folders contain files.'
                execute = False
            # Replace the list of folders from last iteration with the
            # list of folders of the next.
            self.folder_list = new_folder_list

    def copy_structure(self, old_prefix, new_prefix):
        r"""
        Copies a folder structure from old_prefix to new_prefix. To assure all
        folders exists before working with or copying data.

        Parameters
        ----------

        list : list
            A list containing the relative or absolute paths to files.
        old_prefix : string
            The old path to the folder structure that has to be copied. Has to
            actually appear in all the strings in list.
        new_prefix : string
            The path to where the folder structure is to be copied.

        Notes
        -----

        This is useful if one is working on many files stored in several
        sub-folders retrieved using :py:func:`get_list`.

        Examples
        --------

        Say the folder structure is like this

        >>> ls ../modified/
        co10/
        co21/
        >>> ls ../modfied/co10/
        map1/
        map2/
        >>> ls ../modified/co21/
        map1/
        map2/

        This can be copied to say ../even_more_modified by doing as follows:

        >>> from astrolyze.maptools import *
        >>> list = maptools.get_list(../modified)
        >>> maptools.copy_structure(list, old_prefix='../modified',
        >>>                         new_prefix='../even_more_modified')
        """
        for item in list:
            item = item.replace(old_prefix, new_prefix)
            folder_parts = item.split('/')[:-1]
            string = folder_parts[0]
            folders = [string]
            for i in folder_parts[1:]:
                print i
                string = string + '/' + i
                folders += [string]
            for i in folders:
                print i
                if not os.path.isdir(i):
                    print 'mkdir ' + i
                    os.system('mkdir ' + i)
                else:
                    pass

    def unify_resolutions(self, folder=None, resolution=False):
        '''
        Smoothing all maps to the same resolution; either the maximum
        resolution found in the stack or a given resolution.

        Parameters
        ----------

        resolution : float or list
            This may be either:

             1. A list with three entries, i.e.
                [[minor_fwhm], [major_fwhm], [position_angle]]
             2. A list with two entries, position_angle defaults to 0, i.e.
                [[minor_fwhm], [major_fwhm]]
             3. A float. Same minor, major fwhm pa=0

        folder : string
            The path tot the folder in which the files are to be stored.

        Notes
        -----

        The position angle of the output images is fixed to ``Zero`` and can
        currently not be modified.

        Depending on the map units different scaling normalizations have to be
        used after smoothing such that the output units are correct. This
        function tries to deduce the scaling by itself based on the unit that
        is given in the map name. Be sure that this is correct otherwise the
        flux values may be wrong in the output images.
        '''
        folder = check_folder_state(folder)
        # Deduce the maximum resolution of the maps in the Stack.
        max_major_fwhm = max(self.resolutions, key=lambda x: x[0])[0]
        max_minor_fwhm = max(self.resolutions, key=lambda x: x[1])[1]
        # Not well implemented the position angle is fixed.!!
        pa = 0
        new_stack = []
        # This is executed if no specific output resolution is given.
        if not resolution:
            for map_ in self.stack:
                # Check the necessary scaling.
                if map_.fluxUnit.upper() in ['JYB',]:
                    scaling = ''
                if map_.fluxUnit.upper() in ['TMB', 'T', 'KKMS']:
                    scaling = '0.0'
                old_resolution = map_.resolution
                if map_.resolution != 'uk':
                    if folder is not None:
                        map_.prefix = folder
                    if map_.dataFormat not in self.miriad_formats:
                        if folder is not None:
                            map_.prefix = folder
                        map_ = map_.toMiriad()
                    elif map_.dataFormat in self.miriad_formats:
                        if folder is not None:
                            os.system('cp {0} {1}'.format(map_.map_name,
                                                          folder))
                            map_.prefix = folder
                        map_.map_name = map_.returnName()
                    old_map = map_.map_name
                    if (float(max_major_fwhm) != float(map_.resolution[0]) and
                        float(max_minor_fwhm) != float(map_.resolution[1])):
                        map_ = map_.smooth([max_major_fwhm, max_minor_fwhm,
                                            pa], scale=scaling)
                        os.system('rm -rf ' +
                                  str(map_.returnName(dataFormat='fits')))
                        map_ = map_.toFits()
                        new_stack += [map_]
                        os.system('rm -rf ' + map_.returnName(dataFormat=None))
                        os.system('rm -rf ' +
                                  map_.returnName(dataFormat=None,
                                                  resolution=old_resolution))
                    elif (float(max_major_fwhm) == float(map_.resolution[0])
                          and
                          float(max_minor_fwhm) == float(map_.resolution[1])):
                        os.system('rm -rf ' +
                                  str(map_.returnName(dataFormat='fits')))
                        map_ = map_.toFits()
                        new_stack += [map_]
                        os.system('rm -rf {}'.format(
                            map_.returnName(dataFormat=None))
                        )
        if resolution:
            for map_ in self.stack:
                try:
                    if (resolution[0] > map_.resolution[0] or resolution[1] >
                        map_.resolution[1]):
                        print (map_.mapName + ' : resolution higher than new '
                               'one. --> Skip! ')
                        continue
                except:
                    try:
                        if (resolution > map_.resolution[0] or resolution >
                            map_.resolution[1]):
                            print (map_.mapName + ' : resolution higher than '
                                   'new one. --> Skip!')
                            continue
                    except:
                        pass
                if folder is not None:
                    map_.prefix = folder
                # change to miriad and save in new folder
                map_ = map_.toMiriad()
                if map_.fluxUnit.upper() in ['JYB',]:
                    scaling = ''
                if map_.fluxUnit.upper() in ['TMB', 'T', 'KKMS']:
                    scaling = '0.0'
                map_ = map_.smooth(resolution, scale=scaling)
                os.system('rm -rf ' + map_.returnName(dataFormat='fits'))
                map_ = map_.toFits()
                new_stack += [map_]
                os.system('rm -rf ' +  map_.returnName(dataFormat=None))
        return new_stack

    def unify_units(self, unit='JyB', folder=None, debug=True):
        r"""
        Change all maps in the stack to the same unit.

        Parameters
        ----------

        unit : string
            See :py:func:`astrolyze.maps.fits.FitsMap.change_units`
            for details.
        folder : string
            The target folder. By default the maps
            are put into their current folder.
        """
        # Assure that the folder-string has a '\' at the end.
        folder = check_folder_state(folder)
        # The stack has to be updated thus create a new container.
        new_stack = []
        for map_ in self.stack:
            # The map format has to be fits for this to work.
            if map_.dataFormat not in self.fits_formats:
                if folder is not None:
                    map_.prefix = folder
                map_ = map_.toFits()
            elif map_.dataFormat in self.fits_formats:
                map_.prefix = folder
                map_ = map_.update_file()
            old_map = map_.map_name
            print 'old', old_map
            map_ = map_.change_unit(unit, debug=debug)
            print map_.map_name, old_map
            if old_map.strip() != map_.map_name.strip():
                os.system('rm -rf {}'.format(old_map))
            new_stack += [map_]
        return new_stack

    def unify_dimensions(self, template=None, folder=None):
        r""" Reproject all maps to the same central coordinates and map
        dimensions.  All properties of one template map are copied to all the
        other maps using the ``"reproject"`` task of GILDAS.

        Parameters
        ----------

        template : string
            Path to the map that will serve as a template, it may be one of
            the input maps of the stack.
        folder : string
            Path to the folder where the output maps are stored.
        """
        # Assure that the folder-string has a '\' at the end.
        folder = check_folder_state(folder)
        # Load the template map and convert it to Gildas if necessary.
        template = self.get_map_format(template)
        if template.dataFormat not in self.gildas_formats:
            template = template.toGildas()
        # The stack has to be updated thus create a new container.
        new_stack = []
        # Iteratively convert all maps.
        for map_ in self.stack:
            if map_.dataFormat not in self.gildas_formats:
                if folder is not None:
                    map_.prefix = folder
                map_ = map_.toGildas()
            elif map_.dataFormat in self.gildas_formats:
                os.system('cp {0} {1}'.format(map_.map_name, folder))
                if folder is not None:
                    map_.prefix = folder
                map_.map_name = map_.returnName()
            old_map = map_.map_name
            map_ = map_.reproject(template=template.map_name)
            map_ = map_.toFits()
            while len(map_.data) == 1:
                map_.data = map_.data[0]
            max_value = max(
                map_.data[np.where(np.invert(np.isnan(map_.data)))]
            )
            min_value = min(
                map_.data[np.where(np.invert(np.isnan(map_.data)))]
            )
            map_.header["NAXIS"] = 2
            for i in [3, 4]:
                del map_.header['NAXIS{}'.format(i)]
                del map_.header['CDELT{}'.format(i)]
                del map_.header['CRPIX{}'.format(i)]
                del map_.header['CRVAL{}'.format(i)]
                del map_.header['CTYPE{}'.format(i)]
                del map_.header['CROTA{}'.format(i)]
            map_.header['DATAMAX'] = max_value
            map_.header['DATAMIN'] = min_value
            map_ = map_.update_file()
            new_stack += [map_]
            os.system('rm -rf {}'.format(old_map))
            os.system('rm -rf {}'.format(map_.returnName(dataFormat='gdf')))
        return new_stack

    def unify_projections(self, coordinate=None, angle=None, folder=None):
        r"""
        Changing the central coordinate and the rotation angle.
        """
        # Assure that the folder-string has a '\' at the end.
        folder = check_folder_state(folder)
        # The stack has to be updated thus create a new container.
        new_stack = []
        # Iteratively convert all maps.
        for map_ in self.stack:
            if map_.dataFormat not in self.gildas_formats:
                if folder is not None:
                    map_.prefix = folder
                map_ = map_.toGildas()
            elif map_.dataFormat in self.gildas_formats:
                os.system('cp {0} {1}'.format(map_.map_name, folder))
                if folder is not None:
                    map_.prefix = folder
                map_.map_name = map_.returnName()
            old_map = map_.map_name
            if coordinate:
                map_ = map_.reproject(coord=coordinate)
            if angle:
                map_ = map_.goRot(angle)
            map_ = map_.toFits()
            new_stack += [map_]
            os.system('rm -rf {}'.format(old_map))
            os.system('rm -rf {}'.format(map_.returnName(dataFormat='gdf')))
        return new_stack

    def unify_formats(self, target_format='fits', folder=None):
        r""" Converts all maps to the ``target_format`` which has to be one
        that is known to astrolyze. Please see TODO for a list of kwon formats.
        Default is fits.

        Parameters
        ----------

        target_format : string
            The format to which all maps are converted.
        folder : string
            Path to the folder where the new maps are stored.
            If ``None`` the current folder will be used.
        """
        folder = check_folder_state(folder)
        # Check the target_format and allow alternative format strings.
        if target_format in self.gildas_formats:
            target_format = self.gildas_formats
        if target_format in self.fits_formats:
            target_format = self.fits_formats
        if target_format in self.miriad_formats:
            target_format = self.fits_formats
        # Iterate over the Stack and update it.
        new_stack = []
        for i in self.stack:
            if i.dataFormat not in target_format:
                if folder is not None:
                    i.prefix = folder
                if target_format == self.fits_formats:
                    map_ = i.toFits()
                if target_format == self.gildas_formats:
                    map_ = i.toGildas()
                if target_format == self.miriad_formats:
                    map_ = i.toMiriad()
            elif i.dataFormat in target_format:
                if folder is not None:
                    os.system('cp {} {}'.format(i.map_name, folder))
                    map_ = self.get_map_format(i.returnName(prefix=folder))
                if folder is None:
                    map_ = self.get_map_format(i.map_name)
            new_stack += [map_]
        return new_stack

    def pixel_pixel_compare(self, folder=None, plot=False, tol=1e6):
        r""" Producing a pixel-to-pixel comparison for all combinations
        possible between the maps of the stack.

        Parameters
        ----------

        folder : string
           Path to the folder where the text files with the pixel
           to pixel comparisons are stored

        plot : [True | False]
           Decides whether to produce pixel-to-pixel plots directly.

        tol : float
           The tolerance for the maximum difference between the
           values of the pixel of two compared maps. Default to 1e6.

           Returns
           -------

           Created text files in the specified folders that contain two columns
           with the pixel-to-pixel comparisons. These can be used
        """

        def _scatter_plot(folder):
            r""" This plots the pixel_to_pixel comparisons.

            Parameters
            ----------

            folder : string
                Path to the folder where the images should be stored.
            """
            folder = check_folder_state(folder)
            colors = ['green', 'red', 'black', 'yellow', 'blue', 'navy']
            os.system('ls > temp')
            filein = open('temp').readlines()
            if 'scatterPlots\n' in filein:
                os.system('rm -rf scatterPlots/*.*')
                pass
            else:
                os.system('mkdir scatterPlots')
            os.system('rm temp')
            for j in list:
                plt.clf()
                os.system('ls ' + str(j) + ' > temp')
                filein = open('temp').readlines()
                a = 0
                for i in filein:
                    name = (i.split('.')[0].split('_')[0] + '_' +
                            i.split('.')[0].split('_')[1])
                    region = i.split('.')[0].split('_')[2]
                    filein = open(j + '/' + i.strip()).readlines()
                    x = []
                    y = []
                    for k in filein:
                        x += [k.split()[0]]
                        y += [k.split()[1]]
                    plt.plot(x, y, 'x', color=colors[a], label=region)
                    a = a + 1
                plt.legend(numpoints=1, loc='lower right')
                plt.xlabel(str(name.split('_')[0]))
                plt.ylabel(str(name.split('_')[1]))
                plt.savefig('scatterPlots/' + name + '.eps', dpi=None,
                            facecolor='w', edgecolor='w',
                            orientation='portrait', papertype='a4',
                            format='eps', transparent=False,
                            bbox_inches='tight')

        folder = check_folder_state(folder)
        self.unify_formats()
        for i in range(len(list)-1):
            item = list.pop()
            i1 = item.strip().split('.')
            dataFormat = i1.pop()
            print dataFormat
            for j in list:
                i1 = j.strip().split('.')
                print j
                dataFormat = i1.pop()
                print dataFormat
                if dataFormat == 'fits':
                    map = mapClassFits.fitsMap(j)
                    map = map.toMiriad()
                    j = map.mapName
                if dataFormat == 'gdf':
                    map = mapClassGildas.gildasMap(j)
                    map = map.toMiriad()
                    j = map.mapName
                first = mapClassMiriad.miriadMap(str(item))
                second = mapClassMiriad.miriadMap(str(j))
                os.system('mkdir ' + str(folder) + '/' + str(one.species) +
                          '_' + str(two.species))
                os.system('imcmp in1=' + str(item) + ' in2=' + str(j) +
                          'log='+str(folder) + '/' + str(one.species) +
                          '_'+str(two.species) + '/' + str(one.species) + '_' +
                          str(two.species) + '_' + str(one.comments[0]) +
                          '.log tol='+str(tol))


# Some help functions


def check_folder_state(folder):
    if folder is not None and '/' not in folder[-1]:
        folder = folder + '/'
    # If it does not exist we create it.
    if folder is not None and not os.path.isdir(folder):
        os.system('mkdir {}'.format(folder))
    return folder
