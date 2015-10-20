# Copyright (C) 2012, Christof Buchbender
# BSD License (License.txt)
import math
import os
import string
import sys
import matplotlib.pyplot as plt

import astrolyze.maps.main
import astrolyze.maps.fits
import astrolyze.maps.gildas
import astrolyze.maps.miriad

import astrolyze.functions.constants as const


def get_list(folder, data_format=None, depth=False):
    r"""
    Loading a list of files in all subfolders.

    Parameters
    ----------

    folder : string
        The path to the folder that has to be parsed.
    data_format : string
        Search for specific files containing the string, e.g.
        '.fits'
    depth : integer
        The steps into the subfolder system. Defaults to maximum depth.

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
        os.system('ls ' + str(folder) + ' > temp')
        filein = open('temp').readlines()
        list = []
        for i in filein:
            list += [str(folder) + '/' + str(i.strip())]
        os.system('rm -r temp')
        return list

    folder_list = [folder]
    final_list = []
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
        for folder in folder_list:
            list = get_sub_folder(folder)
            for sub_folder in list:
                if os.path.isdir(sub_folder):
                    new_folder_list += [sub_folder]
                    folder_check = True
                if os.path.isfile(sub_folder):
                    if ((data_format is not None and data_format
                         in sub_folder)):
                        if sub_folder not in final_list:
                            final_list += [sub_folder]
                    if (data_format is None):
                            final_list += [sub_folder]
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

        folder_list = new_folder_list
    if not depth:
        return final_list
    if depth:
        return final_list, folder_list


def copy_structure(list, old_prefix, new_prefix):
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

    This is usefull if one is working on many files stored in several
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


def unifyProjection(list, coordinate, angle=None, prefix='repro'):
    os.system('ls > temp')
    filein = open('temp').readlines()
    if folder + '\n' in filein:
        os.system('rm -rf ' + folder + '/*.*')
        pass
    else:
        os.system('mkdir ' + folder)
    os.system('rm temp')
    for i in list:
        i1 = i.strip().split('.')
        dataFormat = i1.pop()
        if dataFormat == 'fits':
            map = FitsMap(i)
            map = map.toGildas(prefix=folder)
        if dataFormat == 'gdf':
            map = GildasMap(i)
            map = map.toGildas(prefix=folder)
        if dataFormat != 'gdf':
            if dataFormat != 'fits':
                map = MiriadMap(str(i.strip()))
                map = map.toGildas(prefix=folder)
        map = map.reproject(coord=coordinate)
        if angle:
            map.goRot(angle)
        map.toFits()
    os.system('rm -rf ' + folder + '/*.gdf')


def unifyUnits(list, folder='units'):
    '''
    NOT READY YET!
    '''
    os.system('ls > temp')
    filein = open('temp').readlines()
    if folder + '\n' in filein:
        os.system('rm -rf ' + folder + '/*.*')
        pass
    else:
        os.system('mkdir ' + folder)
    os.system('rm temp')
    for i in list:
        i1 = i.strip().split('.')
        dataFormat = i1.pop()
        if dataFormat == 'fits':
            maps += [fits.FitsMap(str(i.strip()))]
        if dataFormat == 'gdf':
            maps += [gildas.GildasMap(str(i.strip()))]
        if dataFormat != 'gdf':
            if dataFormat != 'fits':
                maps += [miriad.MiriadMap(str(i.strip()))]
    for i in maps:
        if i.fluxUnit == 'JyB':
            if finalUnit == 'JyB':
                os.system('cp -r ' + str(i.map_name) + ' ' + folder + '/')
                pass
            if finalUnit == 'Tmb':
                i = i.toFits()
                i.data = i.data
    pass


def unifyMaps(list, tinMap, folder='reg'):
    '''
    changes the dimensions and pixel sizes off all maps to that of a template
    map.
    '''
    os.system('ls > temp')
    filein = open('temp').readlines()
    if folder + '\n' in filein:
        os.system('rm -rf ' + folder + '/*.*')
        pass
    else:
        os.system('mkdir ' + folder)
    os.system('rm temp')
    maps = []
    for i in list:
        print i
        i1 = i.strip().split('.')
        dataFormat = i1.pop()
        if dataFormat == 'fits':
            maps += [fits.FitsMap(str(i.strip()))]
        if dataFormat == 'gdf':
            maps += [gildas.GildasMap(str(i.strip()))]
        if dataFormat != 'gdf':
            if dataFormat != 'fits':
                maps += [miriad.MiriadMap(str(i.strip()))]
    i1 = tinMap.strip().split('.')
    dataFormat = i1.pop()
    if dataFormat == 'fits':
        tinMap1 = fits.FitsMap(tinMap)
        tinMap1 = tinMap1.toMiriad()
    if dataFormat == 'gdf':
        tinMap1 = gildas.GildasMap(tinMap)
        tinMap1 = tinMap1.toMiriad()
    if dataFormat != 'gdf':
        if dataFormat != 'fits':
            tinMap1 = miriad.MiriadMap(tinMap)
    for i in maps:
        i.prefix = str(folder) + '/'
        mirMap = i.toMiriad()
        old = mirMap.map_name
        mirMap.regrid(tin=tinMap1.map_name)
        os.system('rm -rf ' + str(old))
        old = mirMap.map_name
        mirMap = mirMap.toFits()
        os.system('rm -rf ' + str(addOn.myStrip(mirMap.map_name, 5)))
    os.system('rm -rf ' + tinMap1.map_name)


def scatterPlot(list):
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
        plt.savefig('scatterPlots/' + name + '.eps', dpi=None, facecolor='w',
                    edgecolor='w', orientation='portrait', papertype='a4',
                    format='eps', transparent=False, bbox_inches='tight')


def unifyResolution(liste, resolution=False, folder='smooth', scaling=''):
    '''
    Approved.
    '''
    os.system('ls > temp')
    filein = open('temp').readlines()
    if folder + '\n' in filein:
        os.system('rm -rf ' + folder + '/*')
        os.system('rm -rf ' + folder + '/*.*')
        pass
    else:
        os.system('mkdir ' + folder)
    os.system('rm temp')
    maps = []
    filein = liste
    for i in filein:
        i1 = i.strip().split('.')
        dataFormat = i1.pop()
        if dataFormat == 'fits':
            maps += [fits.FitsMap(str(i.strip()))]
        if dataFormat == 'gdf':
            maps += [gildas.GildasMap(str(i.strip()))]
        if dataFormat != 'gdf':
            if dataFormat != 'fits':
                maps += [miriad.MiriadMap(str(i.strip()))]
    resolutions = []
    if not resolution:
        for i in maps:
            maxRes = max(i.resolution[0:2])
            resolutions += [float(maxRes)]
        for i in maps:
            i.prefix = folder + '/'
            mirMap = i.toMiriad()
            if float(max(resolutions)) != float(mirMap.resolution):
                sm = mirMap.resolution[0]
                mirMap.smooth(max(resolutions), mirMap.resolution, scaling)
                os.system('rm -rf ' + str(mirMap.returnName(resolution=sm)))
            mirMap.toFits()
    else:
        for i in maps:
            if i.resolution is list:
                maxRes = max(i.resolution[0:2])
            else:
                maxRes = i.resolution
            if float(resolution) > float(max(maxRes)):
                # change the intrinsic folder
                i.prefix = folder + '/'
                # change to miriad and save in new folder
                i = i.toMiriad()
                # store old resolution
                sm = i.resolution
                # smooth to new resolution
                i.smooth(resolution, i.resolution[0], scale=scaling)
                os.system('rm -rf ' + str(i.returnName(resolution=sm)))
                i = i.toFits()
                os.system('rm -rf ' + str(addOn.myStrip(i.map_name, 5)))


def unifyPixelSize(pixSize, liste, folder):
    maps = []
    filein = open(liste).readlines()
    for i in filein:
        i1 = i.strip().split('.')
        if i1[1] == 'fits':
            maps += [fits.FitsMap(str(i.strip()))]
        if i1[1] == 'gdf':
            maps += [gildas.GildasMap(str(i.strip()))]
        if len(i1) == 1:
            maps += [miriad.MiriadMap(str(i.strip()))]
    pixSizes = []
    for i in maps:
        mirMap = i.toMiriad()
        mirMap.regridMiriadToArcsec(pixSize, scale='')
