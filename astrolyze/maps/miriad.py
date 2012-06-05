# Copyright (C) 2012, Christof Buchbender
# BSD Licencse
import math
import os
import string
import sys
import pyfits
import pywcs

from pysqlite2 import dbapi2 as sqlite
from scipy.ndimage import gaussian_filter

# My own modules

import main
import fits
import gildas

import astrolyze.functions.constants as const

class MiriadMap(main.Map):

    def __init__(self, mapName, nameConvention=True):
        main.Map.__init__(self, mapName, nameConvention)
        if self.dataFormat is not '':
            print 'Exiting: not the right format'
            print map.resolution
            sys.exit()

    def toFits(self):
        r"""
        Converts the actual map to a Fits map.

        Returns
        -------

        FitsMap Object.

        Examples
        --------

        With:

        >>> map = miriadMap('M33_MIPS_24mum_JyB_5')
        >>> map = map.toFits()

        it is possible to continue working with the Fits map, using
        :class:`maps.fits.FitsMap` class.
        """
        os.system('rm '+str(self.returnName(dataFormat='fits')))
        print 'rm '+str(self.returnName(dataFormat='fits'))
        string = ('fits in=' + str(self.mapName) + ' out=' +
                  str(self.returnName(dataFormat='fits')) + ' op=xyout')
        print string
        os.system(string)
        self.fitsName = str(self.returnName())+'.fits'
        return fits.FitsMap(str(self.returnName())+'.fits')

    def toGildas(self):
        r"""
        Converts the actual map to a Gildas map.

        Returns
        -------

        GildasMap Object.

        Examples
        --------

        With:

        >>> map = miriadMap('M33_MIPS_24mum_JyB_5')
        >>> map = map.toGildas()

        it is possible to continue working with the Fits map, using
        :class:`maps.gildas.GildasMap` class.
        """
        self.toFits()
        os.system('rm '+str(self.returnName())+'.gdf')
        self.convFile = open('temp.greg','w')
        self.convFile.write('fits '+str(self.returnName(dataFormat = 'fits')) +
                            ' to ' + str(self.returnName(dataFormat = 'gdf')) +
                            '\nexit\n')
        self.convFile.close()
        os.system('greg -nw @temp.greg')
        self.gildasName = str(self.mapName)+'.gdf'
        os.system('rm temp.greg')
        return gildas.GildasMap(self.gildasName)

    def toMiriad(self):
        r"""
        Copies the actual map changing the name such that it takes changes in
        keywords into account.

        Returns
        -------

        MiriadMap Object.

        Examples
        --------

        With:

        >>> map = miriadMap('M33_MIPS_24mum_JyB_5')
        >>> map = map.toMiriad()

        it is possible to continue working with the Miriad map, using
        :class:`maps.gildas.MiriadMap` class.
        """
        os.system('cp -rf '+str(self.mapName)+' '+str(self.returnName()))
        print 'cp -rf '+str(self.mapName)+' '+str(self.returnName())
        self.mapName=self.returnName()
        return self

    def smooth(self, new_resolution, old_resolution=None, scale='0.0'):
        r"""
        Smooths a miriad map to the new resolution. TODO Update to
        new reoslution scheme!!!!!

        Parameters
        ----------

        new_resolution: float or list
            The resolution in of the smoothed image. 
            Can be a:
                * float: Output beam has same major and minor axis [arcsec] and
                         the position angle (PA) [degrees] is 0.
                * A list with two entries:
                        The major and minor axis. PA is 0. 
                        E.g. [major_axis, minor_axis ]
                * A list with three entries:
                        [major_axis, minor_axis, PA] 
        old_resolutio: float
            If None the self.resolution information is taken into account. 
        scale:string

        Returns
        -------

        MiriadMap Object:
            The smoothed image.
        """
        if old_resolution == None:
            oldMajor = self.resolution[0]
            oldMinor = self.resolution[1]
            pa = self.resolution[2]
        if old_resolution != None:
            if old_resolution is list:
                if len(old_resolution) == 2:
                    oldMajor = old_resolution[0]
                    oldMinor = old_resolution[1]
                    pa = 0
                if len(old_resolution) == 3:
                    oldMajor = old_resolution[0]
                    oldMinor = old_resolution[1]
                    pa = old_resolution[2]
            if old_resolution is not list:
                oldMajor = old_resolution
                oldMinor = old_resolution
                pa = 0
        if (float(oldMajor) > float(new_resolution) or float(oldMinor) >
           float(new_resolution)):
            print 'Error: Old Resolution bigger than new one!'
        # calculate the fwhm for the convolving gaussian    
        fwhmMajor = math.sqrt(float(new_resolution)**2-float(oldMajor)**2)
        fwhmMinor = math.sqrt(float(new_resolution)**2-float(oldMinor)**2)
        print fwhmMajor, fwhmMinor
        os.system('rm -rf '+str(self.returnName(resolution=[
                  float(new_resolution), float(new_resolution), 0.0])))
        if scale != '':
            executeString = ('smooth in='+str(self.mapName)+' '
        'out='+str(self.returnName(resolution=[float(new_resolution),
                   float(new_resolution), 0.0]))+' '
        'fwhm='+str('%.2f'%(fwhmMajor))+', '+str('%.2f'%(fwhmMinor))+' '
        'pa='+str(pa)+' scale='+str(scale))
            print executeString
            os.system(executeString)
        else:
            executeString =  ('smooth in='+str(self.mapName)+' '
        'out='+str(self.returnName(resolution=[float(new_resolution),
                   float(new_resolution), 0.0]))+' '
        'fwhm='+str('%.2f'%(fwhmMajor))+', '+str('%.2f'%(fwhmMinor))+' '
        'pa='+str(pa)+' scale='+str(scale))
            print executeString
            os.system(executeString)
        # set the mapName and the resolution to the new values<
        self.mapName = self.returnName(resolution=[float(new_resolution), 
                                                   float(new_resolution), 0.0])
        self.resolution = [float(new_resolution), float(new_resolution), 0.0]
        return MiriadMap(self.returnName())

    def _moment(self, iN='', region='', out='', mom='0', axis='', 
               clip='', rngmsk='', raw=''):
        '''
        Wrap around MIRIADs moment task. 
        keywords are as in miriad. 
        By default (-> if you give no arguments to the function)
        it creates the zeroth moment of the map
        '''
        fileout=open('miriad.out','a')
        string = 'moment '
        if iN == '':
            iN= self.mapName
        string += 'in='+str(iN)+' '
        if region != '':
            string += 'region='+str(region)+' '
        if out=='':
            self.newComments = []
            if 'cube' in self.comments:
                print 'yes1'
                for i in self.comments:
            
                    if str(i) == 'cube':
                        print 'yes1'
                        self.newComments += ['mom'+str(mom)]
                    else:
                        self.newComments += [str(i)]
                self.comments = self.newComments
            else:
                self.comments += ['mom'+str(mom)]
            print self.comments
            out = self.returnName()
        string += 'out='+str(out)+' '
        string += 'mom='+str(mom)+' '
        if axis != '':   
            string += 'axis='+str(axis)+' '
        if clip != '':   
            string += 'clip='+str(clip)+' '
        if rngmsk != '':   
            string += 'rngmsk='+str(rngmsk)+' '
        if raw != '':   
            string += 'raw='+str(raw)+' '

        os.system('rm -rf '+str(self.returnName()))
        print string
        os.system(string)
        fileout.write(string+'\n')
        self.mapName = out
        fileout.close()  


    def regrid(self, iN='', out='', axes='1,2', tin='', desc='', 
               options='', project='', rotate='', tol=''):
        fileout=open('miriad.out','a')
        string = 'regrid '
        if iN == '':
            iN= self.mapName
        string += 'in='+str(iN)+' '
        if out=='':
            self.comments += ['regrid']
            out= self.returnName()
        string += 'out='+str(out)+' '
        if axes != '':
            string += 'axes='+str(axes)+' '
        if tin != '':
            string += 'tin='+str(tin)+' '
        if desc != '':   
            string += 'desc='+str(desc)+' '
        if options != '':   
            string += 'options='+str(options)+' '
        if project != '':   
            string += 'project='+str(project)+' '
        if  rotate != '':   
            string += 'rotate='+str(rotate)+' '
        if tol != '':   
            string += 'tol='+str(tol)+' '

        os.system('rm -rf '+str(self.returnName()))
        print string
        os.system(string)
        fileout.write(string+'\n')
        self.mapName = out
        fileout.close()



    def ellipseMask(self, pa, incline, radius, coord, out, 
                    pix_or_coord='coord', logic='lt'):
        tempFitsMap = self.toFits()
        xy = tempFitsMap.sky2xy(coord)
        x0 = str(int(floor(float(xy[0]))))
        y0 = str(int(floor(float(xy[1]))))
        os.system('rm -rf '+out) 
        print '#################'
        print self.inclination
        print self.pa
        print '#################'
        sineCosArg = str(float(2*math.pi/360*(-90+float(self.pa)))) 
        inclinationRadian = str(2*math.pi/360*self.inclination)
        os.system('maths \'exp=<'+self.mapName+'>\'  \'mask=sqrt((((x-'+x0+')*'
                  'cos('+sineCosArg+'))-((y-'+y0+')*sin('+sineCosArg+')))**2+'
                  '(((x-'+x0+')*sin('+sineCosArg+'))+((y-'+y0+')*cos'
                  '('+sineCosArg+')))**2/(cos('+inclinationRadian+')**2)).'+
                  logic+'.'+str(radius*60/math.sqrt((
                  float(tempFitsMap.header['CDELT1'])/(1./60/60))**2))+'\' '
                  'out='+out+' xrange=0,'+str(tempFitsMap.header['NAXIS1'])+' '
                  'yrange=0,'+str(tempFitsMap.header['NAXIS2']))
        self.mapName = out


    def _regridMiriadToArcsec(self, value, JyB_KkmS='KkmS'):
        fitsFile = self.toFits()
        self.naxis1 = float(fitsFile.header['naxis1'])
        self.naxis2 = float(fitsFile.header['naxis2'])
        self.cdelt1 = float(fitsFile.header['CDELT1'])
        self.cdelt2 = float(fitsFile.header['CDELT2'])
        self.crval1 = float(fitsFile.header['CRVAL1'])
        self.crval2 = float(fitsFile.header['CRVAL2'])
        self.crpix1 = float(fitsFile.header['CRPIX1'])
        self.crpix2 = float(fitsFile.header['CRPIX2'])
        
        print (self.naxis1, self.naxis2, self.cdelt1, self.cdelt2, self.crval1, 
               self.crval2, self.crpix1, self.crpix2)

        self.cdelt1arcs = float(self.cdelt1)/(value*(1./60/60))
        self.cdelt2arcs = float(self.cdelt1)/(value*(1./60/60))
        self.naxis1New = self.naxis1*math.sqrt(self.cdelt1arcs**2)
        self.crpix1New = self.naxis1New/(self.naxis1/self.crpix1)
        self.naxis2New = self.naxis2*math.sqrt(self.cdelt2arcs**2)
        self.crpix2New = self.naxis2New/(self.naxis2/self.crpix2)
    
        self.newPix = float(value)/60/60
        print 'oldPixel',self.cdelt1*60*60,self.cdelt2*60*60
        print 'axis1',self.naxis1,'->',self.naxis1New
        print 'axis2',self.naxis2,'->',self.naxis2New
        
        self.regrid(desc=str(self.crval1*2*math.pi/360)+','+
                    str(self.crpix1New)+','+
                    str(float(value)*-2*math.pi/360/60/60)+','+
                    str(self.naxis1New)+','+str(self.crval2*2*math.pi/360)+','+
                    str(self.crpix2New)+','+
                    str(float(value)*2*math.pi/360/60/60)+','+
                    str(self.naxis2New))

        if JyB_KkmS == 'KkmS':
            os.system('rm -rf '+str(self.returnName())+'_norm')
            os.system('maths \'exp=<'+str(self.returnName())+'>/'+str(self.cdelt1arcs*self.cdelt2arcs)+'\' out='+str(self.returnName())+'_norm')
            self.comments += ['norm']
            self.mapName = self.returnName()

