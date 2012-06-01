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

from astrolyze.setup.paths import prefix
import astrolyze.functions.constants as const


class MiriadMap(main.Map):

    def __init__(self, mapName, nameConvention=True):
        main.Map.__init__(self, mapName, nameConvention)
        if self.dataFormat is not '':
            print 'Exiting: not the right format'
            print map.resolution
            sys.exit()

    def toFits(self):
        os.system('rm '+str(self.returnName())+'.fits')
        print 'rm '+str(self.mapName)+'.fits'
        string = ('fits in='+str(self.mapName)+' out='+str(self.returnName())+'.'
                 'fits op=xyout')
        print string
        os.system(string)
        self.fitsName = str(self.returnName())+'.fits'
        return fits.FitsMap(str(self.returnName())+'.fits')

    def toGildas(self):
        self.toFits()
        os.system('rm '+str(self.returnName())+'.gdf')
        self.convFile = open('temp.greg','w')
        self.convFile.write('fits '+str(self.mapName)+'.fits to '+
                            str(self.returnName())+'.gdf\n'
                            'exit\n')
        self.convFile.close()
        os.system('greg -nw @temp.greg')
        self.gildasName = str(self.mapName)+'.gdf'
        os.system('rm temp.greg')
        return gildas.GildasMap(self.gildasName)

    def toMiriad(self):
        os.system('cp -rf '+str(self.mapName)+' '+str(self.returnName()))
        print 'cp -rf '+str(self.mapName)+' '+str(self.returnName())
        self.mapName=self.returnName()
        return self

    def smooth(self, newRes, oldRes=None, scale='0.0'):
        if oldRes == None:
            oldMajor = self.resolution[0]
            oldMinor = self.resolution[1]
            pa = self.resolution[2]
        if oldRes != None:
            if oldRes is list:
                if len(oldRes) == 2:
                    oldMajor = oldRes[0]
                    oldMinor = oldRes[1]
                    pa = 0
                if len(oldRes) == 3:
                    oldMajor = oldRes[0]
                    oldMinor = oldRes[1]
                    pa = oldRes[2]
            if oldRes is not list:
                oldMajor = oldRes
                oldMinor = oldRes
                pa = 0
        if float(oldMajor) > float(newRes) or float(oldMinor) > float(newRes):
            print 'Error: Old Resolution bigger than new one!'
        # calculate the fwhm for the convolving gaussian    
        fwhmMajor = math.sqrt(float(newRes)**2-float(oldMajor)**2)
        fwhmMinor = math.sqrt(float(newRes)**2-float(oldMinor)**2)
        print fwhmMajor, fwhmMinor
        os.system('rm -rf '+str(self.returnName(resolution=[
                                                float(newRes), float(newRes), 
                                                0.0])))
        if scale != '':
            executeString = ('smooth in='+str(self.mapName)+' '
        'out='+str(self.returnName(resolution=[float(newRes), 
                                               float(newRes), 0.0]))+' '
        'fwhm='+str('%.2f'%(fwhmMajor))+', '+str('%.2f'%(fwhmMinor))+' '
        'pa='+str(pa)+' scale='+str(scale))
            print executeString
            os.system(executeString)
        else:
            executeString =  ('smooth in='+str(self.mapName)+' '
        'out='+str(self.returnName(resolution=[float(newRes), 
                                               float(newRes), 0.0]))+' '
        'fwhm='+str('%.2f'%(fwhmMajor))+', '+str('%.2f'%(fwhmMinor))+' '
        'pa='+str(pa)+' scale='+str(scale))
            print executeString
            os.system(executeString)
        # set the mapName and the resolution to the new values<
        self.mapName = self.returnName(resolution=[float(newRes), 
                                                   float(newRes), 0.0])
        self.resolution = [float(newRes), float(newRes), 0.0]
        return mapClassMiriad.MiriadMap(self.returnName())

    def moment(self, iN='', region='', out='', mom='0', axis='', 
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


    def regridMiriadToArcsec(self, value, JyB_KkmS='KkmS'):
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

