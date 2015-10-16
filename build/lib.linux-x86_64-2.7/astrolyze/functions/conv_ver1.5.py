# Version 1.5
# Christof Buchbender

import math, string, os ,sys,getopt

print 'Version 1.5'

# Prints out the information of the command line suffixes that are possible
def usage():
    print "\nTool to convert between Offsets in real (telescope) angles\nand absolute Coordinates in the equatorial System. An additional rotation of \nthe reference frame is possible. \n\n   -h; --help  :   print this help \n  -r; --relative  :   Default. Changes from offsets to absolsute coordinates\n                default input/output files are:     \n              offsets.txt/calcAbsCoords.txt\n     -a; --absolute  :   Calculates offsets from absolute coordinates\n              default input/output files are:\n               absCoords.txt/calcOffsets.txt\n     -i; --input     :   allows to specify input file name \n    -o; --output    :   allows to specify output file name\n    -c; --cli   :   command line interface for one conversion only.\n               The conversion mode has to be assigned with\n               -r/--relative or -a/--absolute\n\n"

def usageRelativeCommand():
    print '++++Command Line Usage++++\n\nReference Coordinates:\n\n The Reference coordinates in RA DEC has to be given as follows:\n 01:35:18.8 +30:53:23.4\n The Dec coordinate value can hava a positive \"+\" or a negative \"-\" sign.\n If the sign is missing it is interpretet as positive. \n\nAngle:\n The rotation angle ihas to be given in degrees\n\nOffsets: \n Units can be given after the definition of the offsets. \n E.g. for an offset of 3 and 4 arcseconds in RA and Dec\n the offsets can be specified as follows: \n 3 4 arcsec\n Accepted Units are for:  \n\n arcminutes : arcmin \' \n arcseconds : arcsec \'\' \" \n degrees    : deg \n hours      : hours  h \n minutes    : minutes min m \n seconds    : seconds sec s \n radian     : radian rad'

# If the default files are not present and no file is specified with -i/--input, a default file with content that explaines the nessecary structure of the input files are created.

def exampleFiles(which):
    if which == 'inputRel':
        example = open('offsets.txt','w')
        example.write('01:34:11.8, +30:50:23.4      #Central coordinates, has to be first row!\n0               #angle in degrees if reference system is..\n                #..rotated against the equatorial system..\n                #..has to be second row! \n3 1                  #offsets in RA and Dec in arcminutes\n\n#the default units, which are arcminutes, can be changed to other units.\n\n240 23 arcsec           #possible entries are: arcsec,",'' \n0.05 0.03 deg          #degrees, can be: deg\n2 4 hours            #hourangle, can be: hours or h \n2 3 min                #minutes, can be: minutes, min or m\n2 3 sec                #seconds, can be: seconds, sec or s \n1.16e-3 4.43e-3 rad       #radian, can be: rad, radian ')
        print '\nInput file not found. \nThe default file, \'offsets.txt\', has been created,\nwhich describes the nessecary structure of the \ninput file for the conversion between offsets and coordinates.\n'
        sys.exit(2)
    elif which == 'inputAbs':
        example = open('absCoords.txt','w')
        example.write('01:34:11.8, +30:50:23.4      #central Coordinates, has to be first row \n0               #angle if reference system is rotated.. \n              #..against the equatorial system..\n                #..has to be second row\n01:31:11.8, +30:52:23.4        #first coordinate\n01:35:18.8, +30:53:23.4      #second coordinate and so on')
        print '\nInput file not found. \nThe default file \'absCoords.txt\' has been created \nwhich describes the nessecary structure of the input file\nfor the conversion between absolute coordinates and offsets in real angle.\n'

        sys.exit(2)

def main(argv):
    # Default Usage Converting Offsets to absolute coordinates
    inputFile = 'offsets.txt'
    outputFile = 'calcAbsCoords.txt'
    offset = 1          #switches for the different conversion types
    absolute = 0
    which = 'inputRel'
    helpRelative = 0
    distanceCalc = None
    #Reading the parameters from the command line
    try:
        # first list are the possible short options, letters with dopple point afterward need an argument
        opts, args = getopt.getopt(argv, "hci:o:rad:e:D:", ["help","cli" "input=","output=","relative","absolute","deg","eqa","distance"])
        if '-r' in opts and '-a' in opts:
            usage()
            sys.exit(2)
        if '--relative' in opts and '--absolute' in opts:
            usage()
            sys.exit(2)

    except getopt.GetoptError:      # display the help and exit the programm if wrong parameters are given
        usage()
        sys.exit(2)

    remember=0

    for opt, arg in opts:           #define the reaction to the different parameters

        if opt in ("-h","--help"):
            help = 1
        elif opt in ('-i','--input'):
            inputFile = arg
        elif opt in ('-o','--output'):
            outputFile = arg
        elif opt in ('-c','--cli'):
            offset = 0
            absolute = 0
            remember=1
        elif opt in ('-d','--deg'):
            argnew = arg.split(',')
            print argnew
            print equatorial2Deg(argnew)
        elif opt in ('-e','--eqa'):
            argnew = arg.split(',')
            print argnew
            print deg2Equatorial(argnew)
        elif opt in ('-D','--Distance'):
            distanceCalc = 1
            absolute     = 0
            offset       = 0
            Distance = float(arg)
            inputFile = 'absCoords.txt'
            outputFile = 'calcOffsets.txt'

        elif opt in ('-r','--relative'):
            which='inputRel'
            helpRelative =1
            offset = 1
            if remember == 1:
                offset=0
            absolute = 0
        elif opt in ('-a','--absolute'):
            which='inputAbs'
            inputFile = 'absCoords.txt'
            outputFile = 'calcOffsets.txt'
            offset = 0
            absolute = 1
            if remember == 1:
                absolute=0

    try:
        if help == 1:
            if helpRelative == 1:
                usageRelativeCommand()
                sys.exit()
            else:
                usage()
                sys.exit()
    except:
        print ''

    source = "".join(args)
    try:                    #try to read the input file
        filein = open(inputFile)
    except:                     #in case of failure: if the input file is the default one generate an example. If it is a user specified file name: Print 'File not found'
        if inputFile == 'offsets.txt' or inputFile=='absCoords.txt':
            exampleFiles(which)
            #filein = open(inputFile)
        else:
            print '\nFile not found\n'
            sys.exit(2)
    lines= filein.readlines()
    fileout=open(outputFile,'w')

    #Command line usage relative->absolute
    if offset == 0 and absolute == 0:
        try:
            if which == 'inputRel' and remember==1:

                centralCoord = raw_input('Central coordinates (e.g.: 01:34:11.8 +30:50:23.4): ')
                centralCoord = centralCoord.split()

                angle=input('Insert Rotation Angle (deg)            : ')

                offset = raw_input('Offset (e.g.:3 3 arcsec)                : ')

                offset=offset.split()
                offset[0] = float(offset[0])
                offset[1] = float(offset[1])
                i =offset
                try:
                    s=  i[2]
                except:
                    i = i+['arcmin']
                if i[2] in ["''","arcsec",'"']:
                    print i
                    i[0] = i[0]/60
                    i[1] = i[1]/60
                    print i
                if i[2] in ["deg"]:
                    i[0] = i[0]*60
                    i[1] = i[1]*60
                if i[2] in ["arcmin","'"]:
                    i[0] = i[0]
                    i[1] = i[1]
                if i[2] in ["hours","h"]:

                    i[0] = i[0]*15*60
                    i[1] = i[1]*15*60

                if i[2] in ["minutes,","min","m"]:
                    i[0] = i[0]*15
                    i[1] = i[1]*15
                if i[2] in ["seconds","sec","s"]:
                    i[0] = i[0]*15/60
                    i[1] = i[1]*15/60
                if i[2] in ["radian","rad"]:
                    i[0] = i[0]/2*math.pi*360*60
                    i[1] = i[1]/2*math.pi*360*60
                offset =i


                print '\nThe corresponding coordinates are:\n'+                  relative2Absolute(centralCoord,offset,angle)



        #Command line usage absolute->relative
            if which == 'inputAbs':

                centralCoord = raw_input('Central coordinates (e.g.: 01:34:11.8, +30:50:23.4)   :')
                centralCoord = centralCoord.split(',')

                angle=input('Insert Rotation Angle (deg)                :')
                absoluteCoord = raw_input('Second coordinates (e.g.: 01:34:11.8, +30:50:23.4)   :')
                absoluteCoord = absoluteCoord.split(',')
                print '\nThe corresponding offsets are:\n'+absolute2Relative(centralCoord,absoluteCoord,angle)+' in arcmins'
        except:
                if remember == 1:
                    print '\nCommand Line usage.\nPlease define the conversion type:\n \nrelative->absolute: -r/--relative or \nabsolute->relative: -a/--absolute.\n'
                sys.exit(2)

    #conversion absolute to Offsets



    if (offset == 0 and absolute == 1) or distanceCalc == 1:
        print 'here'
        centralCoord = [str(lines[0].split()[0].strip()),str(lines[0].split()[1].split()[0])] #read first row
        print centralCoord
        angle = float(lines[1].split()[0].strip())  #read secondt row
        print angle

        list=[]

        for i in lines[2:]:         #read the remaining rows if they have correct structure

            i = i.split()

            if len(i)==2:

                i[1]=i[1].split()
                print i[1]
                s=[]                # assure that comments input files are possible with #
                for j in i[1]:


                    if '#' in j:
                        break
                    else:
                        s += [j]
                i[1] = s                # end comments part

            if len(i) < 2:
                continue




            #i[1]=i[1].split()[0]

            # to print out only  the results use (for e.g. further processing in other programs):
            #fileout.write(str(absolute2Relative(centralCoord,[i[0],i[1]],angle))+'\n')
            a,b = absolute2Relative(centralCoord,[i[0],i[1][0]],angle)

            if absolute==1:
                fileout.write(str(i[0])+' '+str(i[1])+' -> '+"%1.3f"%a+' '+"%1.3f"%b+'\n')
            if distanceCalc ==1:
                cDistance= calcDistance(a,b,Distance)
                print cDistance
                list+=["%1.3f"%cDistance]
                fileout.write(str(i[0])+' '+str(i[1])+' -> '+"%1.3f"%float(cDistance)+'\n')

        print list
        print 'Output written to '+str(outputFile)

    #conversion Offsets to absolute
    if offset == 1 and absolute == 0:

        # read the first line of the input file which is the central Coord


        centralCoord = [str(lines[0].split(',')[0].strip()),str(lines[0].split(',')[1].split()[0])]
        #read the second line of the input file which is the angle of the rotated coordinate system. If no Rotation is present the angle is zero.
        angle = float(lines[1].split()[0].strip())




        for i in lines[2:]:


            i = i.split()

            s=[]                # assure that comments input files are possible with #
            for j in i:


                if '#' in j:
                    break
                else:
                    s += [j]
            i = s               # end comments part

            if len(i) < 2:
                continue

            #assure that the Offsets are treatet as floats
            i[0] = float(i[0])
            i[1] = float(i[1])


            # the default offsets are in arcmin but others can be specified:
            # This part calculates the corresponding angle in arcmin
            if len(i) == 3:

                angleType = i[2]

                if i[2] in ["''","arcsec",'"']:
                    i[0] = i[0]/60
                    i[1] = i[1]/60
                if i[2] in ["deg"]:
                    i[0] = i[0]*60
                    i[1] = i[1]*60
                if i[2] in ["arcmin","'"]:
                    continue
                if i[2] in ["hours","h"]:

                    i[0] = i[0]*15*60
                    i[1] = i[1]*15*60

                if i[2] in ["minutes,""min","m"]:
                    i[0] = i[0]*15
                    i[1] = i[1]*15
                if i[2] in ["seconds","sec","s"]:
                    i[0] = i[0]*15/60
                    i[1] = i[1]*15/60
                if i[2] in ["radian","rad"]:
                    i[0] = i[0]/2*math.pi*360*60
                    i[1] = i[1]/2*math.pi*360*60

            #print i
            #print only the results:
            #fileout.write(str(relative2Absolute(centralCoord,[float(i[0]),float(i[1])],angle))+'\n')
            fileout.write(str("%1.3f"%(i[0]))+' '+str("%1.3f"%(i[1]))+' -> '+str(relative2Absolute(centralCoord,[float(i[0]),float(i[1])],angle))+'\n')
        print 'Output written to '+str(outputFile)

# Conversion between equatorial System and Degrees

def equatorial2Deg(coord):


    CoordsplitRA = coord[0].strip().split(':')
    CoordsplitDec = coord[1].strip().split(':')
    sign = CoordsplitDec[0][0]
    if sign != '+' and sign != '-':
        sign='+'
    print CoordsplitRA
    print CoordsplitDec

    if sign == '+':
        deg=[(float(CoordsplitRA[0])*(360./24) +
              float(CoordsplitRA[1])*(360./24/60) +
              float(CoordsplitRA[2])*(360./24/60/60)),
             (float(CoordsplitDec[0]) +
              float(CoordsplitDec[1])*(1./60) +
              float(CoordsplitDec[2])*1./60/60)]
    if sign == '-':
        deg=[(float(CoordsplitRA[0])*(360./24)+float(CoordsplitRA[1])*(360./24/60)+float(CoordsplitRA[2])*(360./24/60/60)),(float(CoordsplitDec[0])-float(CoordsplitDec[1])*(1./60)-float(CoordsplitDec[2])*1./60/60)]


    return deg


#Conversion between Degrees and equatorial Coordinates

def deg2Equatorial(Deg):


    Coord=[]
    Coord+=[str(int(Deg[0]/15))+':'+str(int(((Deg[0]/15)-int(Deg[0]/15))*60))+':'+"%1.2f"%(float(str((((Deg[0]/15-int(Deg[0]/15))*60)-int((Deg[0]/15-int(Deg[0]/15))*60))*60)))]

    Coord+=[(str(int(Deg[1]))+':'+str(int(math.fabs(int((float(Deg[1])-int(Deg[1]))*60))))+':'+"%1.2f"%(math.fabs(float(str(float(((float(Deg[1])-int(Deg[1]))*60)-int((float(Deg[1])-int(Deg[1]))*60))*60)))))]

    return str(Coord[0])+' '+str(Coord[1])



#Calculate the absolute coordinates. It needs the central Coordinates: centralCoord, the Offset in  [OffsetRA[arcmin],OffsetDec[arcmin]] and, if presetnt , the Rotation angel between the equatorial system and the one in which the offsets are calculated.

def relative2Absolute(centralCoord,Offset,angle):
    angle =  angle*2*math.pi/360.

    Offset[0] = Offset[0]/60.
    Offset[1] = Offset[1]/60.

    #Rotation
    xAbs = math.cos(angle)*Offset[0]+math.sin(angle)*Offset[1]
    yAbs = -1*math.sin(angle)*Offset[0]+math.cos(angle)*Offset[1]
    degCoords = equatorial2Deg(centralCoord)
    degCoords[1]=float(degCoords[1])+float(yAbs)

    #correction for Declination
    degCoords[0]=float(degCoords[0])+(float(xAbs)/math.cos(math.radians(degCoords[1])))


    return deg2Equatorial(degCoords)


#calculate the Relative Coordinates (Offsets) from the absolute Coordinates and a central position

def absolute2Relative(centralCoord,absoluteCoord,angle):
    angle =  math.radians(angle)

    print absoluteCoord
    centralDeg=equatorial2Deg(centralCoord)
    absoluteDeg=equatorial2Deg(absoluteCoord)

    diff=[absoluteDeg[0]-centralDeg[0],absoluteDeg[1]-centralDeg[1]]
    diff=[diff[0]*math.cos(math.radians(absoluteDeg[1])),diff[1]]  #correction for declination


    # inverse Rotaion Matrix
    xDiffRot = math.cos(angle)*diff[0]-math.sin(angle)*diff[1]
    yDiffRot = math.sin(angle)*diff[0]+math.cos(angle)*diff[1]
    xDiffRotMinutes = xDiffRot*60
    yDiffRotMinutes = yDiffRot*60
    #return [xDiffRot,yDiffRot]
    return float(xDiffRotMinutes),float(yDiffRotMinutes)




def calcDistance(a,b,sourceDistance,output='kpc'):
    '''
    This function caluclates Distances between two coordinates.
      Input:
        a       :  Offset in Az
        b       :  Offset in Dec
        output      :  choose the unit of the Distance (pc, kpc or Mpc)
        sourceDistance  :  Distance to the source in pc
      Output:
        Distance in kpc
    '''

    Distance = math.sqrt(a**2+b**2)*60*sourceDistance*4.848e-6
    if output=='kpc':
        Distance = Distance/1e3
    elif output == 'pc':
        pass
    elif output=='Mpc':
        Distance =Distance/1e6
    return Distance


#start the main program

if __name__ == "__main__":
    main(sys.argv[1:])
