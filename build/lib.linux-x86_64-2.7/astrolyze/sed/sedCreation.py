import numpy as np
from astrolyze.setup.paths import prefix

def getSED(fileList, coordList, tableFolder=''):
    '''
    This function reads SEDs from the fits maps given in the file fileList
    which has to have the followin format:
    
    frequency1 error1 pathname1
    frequency2 error2 pathname2
    .
    .
    .
    
    The fluxes are read at the coordinates that are specified in the file 
    coordList with the format
    
    nameOfPosition1 RAcoordInDeg1 DECcoordinDeg1    
    nameOfPosition2 RAcoordInDeg2 DECcoordinDeg2
    .
    .
    THIS IS OUTDATED MOVE TO astrolyze.mapclass!!!!!!!
    .
    '''
    if type(fileList) == str:
        files = open(fileList).readlines()
    else:
        files = fileList
    coords = coordList
    if type(coordList) == str:
        coords = open(coordList).readlines()
    
    frequencies = []
    errors = []
    filePaths = []
    sumMode = []
    for i in files:
        if type(fileList) == str:
            i=i.split()
        frequencies += [i[0]]
        errors += [i[1]]
        filePaths += [prefix.sys+i[2]]
        try:
            sumMode += [i[3]]
        except:
            sumMode += ['fluxread'] 
    
    coordsRA = []
    coordsDEC = []
    positionNames = []
    for i in coords:
        if type(coordList) == str:
            i=i.split()
        positionNames += [i[0]]
        coordsRA += [i[1]]
        coordsDEC += [i[2]]
    returnNames = []
    returnFinal = []
    Fluxes=[]
    Errors=[]
    Frequencies=[]

    for i in range(len(coordsRA)):
        fileout=open(str(tableFolder)+'/'+str(positionNames[i])+'.tab','w')
        returnNames += [positionNames[i]]
        x=0
        Fluxes+=[[]]
        Errors+=[[]]
        Frequencies+=[[]]
        
        for j in range(len(filePaths)):
            try:
                if sumMode[x] == 'fluxread':
                    flux = fluxread(filePaths[j],coordsRA[i],coordsDEC[i],'Coord')
                    print flux
                    if invert(isnan(flux)): # checks that the flux is not a 'NAN', means -> if flux not 'nan'
                        fileout.write(str(frequencies[j])+' '+str(float(fluxread(filePaths[j],coordsRA[i],coordsDEC[i],'Coord')))+' '+str(float(fluxread(filePaths[j],coordsRA[i],coordsDEC[i],'Coord'))*float(errors[j]))+'  \n')
                        Frequencies[i]+= [float(frequencies[j])]
                        Fluxes[i] += [float(fluxread(filePaths[j],coordsRA[i],coordsDEC[i],'Coord'))]
                        Errors[i]+= [float(fluxread(filePaths[j],coordsRA[i],coordsDEC[i],'Coord'))*float(errors[j])] 
                if sumMode[x] != 'fluxread':
                    map = fitsMap(filePaths[j])
                    flux = fluxAperture([coordsRA[i],coordsDEC[i]],float(sumMode[x]))
                    print flux
                    if invert(isnan(flux)): # checks that the flux is not a 'NAN', means -> if flux not 'nan'
                        fileout.write(str(frequencies[j])+' '+str(float(fluxAperture([coordsRA[i],coordsDEC[i]],float(sumMode[x]))[0]))+' '+str(float(fluxAperture([coordsRA[i],coordsDEC[i]],float(sumMode[x]))[0]*float(errors[j]))+'  \n'))
                        Frequencies[i]+= [float(frequencies[j])]
                        Fluxes[i] += [float(fluxAperture([coordsRA[i],coordsDEC[i]],float(sumMode[x]))[0])]
                        Errors[i]+= [float(fluxAperture([coordsRA[i],coordsDEC[i]],float(sumMode[x]))[0])*float(errors[j])] 

            except:
                print filePaths[j],'geht nicht'
                pass

            x += 1
        
        fileout.close()
    x=0
    
    for i in returnNames:
        SED=asarray([Frequencies[x],Fluxes[x],Errors[x]])
        returnFinal+=[[returnNames[x],SED]]
        #returnFinal+=[[returnNames[x],[Frequencies[x],Fluxes[x],Errors[x]]]]
        x=x+1   
    
    return returnFinal

def getFluxes(name,tableFolder=''):
    '''
    reads the data written by getSED and provides an array that is undestood 
    by the fitting Algorithms astroFunctions.py
    '''
    filein = open(str(tableFolder)+'/'+str(name)+'.tab').readlines()
    data = [[],[],[]]
    for i in filein:
        
        i = i.split()
        data[0] += [float(i[0])]
        data[1] += [float(i[1])]
        data[2] += [float(i[2])]

    return asarray(data)
                













































































        

        
    
    


