import matplotlib.pyplot as plt
import astrolyze.functions.astro_functions as astFunc
import astrolyze.functions.constants as const
import numpy as np
from copy import deepcopy as copy

def plot_sed(p2, nu_or_lambda='nu', color='black', kappa='easy', 
            linewidth=0.5, xRange='normal'):
    '''
    Plot a multi component greybody model.

    :param p2:  formatted: p2 = [[t1,t2,t3,...], [N1,N2,N3,...], [beta]] 
    :param nu_or_lambda: plot against frequency ``'nu'`` or wavelenght 
        ``'lambda'``
    :param kappa: The kappa to use. ``'easy'`` or ``'Kruegel'``. Please refer
       to :py:func:`functions.astroFunctions.greyBody` for more information.
    :param xRange: PLEASE ADD DESCRIPTION
    :param linewidth: The linewidth of the plotted lines. Default to 0.5.
    :param linewidth: float
    :param color: the color of the plotted lines. Default to ``'black'``
    :type color: matplotlib conform color.
    '''

    if xRange == 'LTIR':
    # Plot the SED in the range of the determination
        # of the L_TIR: 3-1100 micron 
        xmin =  3e-6# micron
        xmax =  4000e-6 # micron
        # conversion to frequency in GHz
        xmin = const.c/xmax/1e9
        xmax = const.c/xmin/1e9
        step = 0.1

    if xRange == 'normal':
        # arbitrary range definition 
        xmin = 1e-2
        xmax = 3e5	
        step = 0.5
    if type(xRange) == list:
        xmin = xRange[0]
        xmax = xRange[1]
        if len(xRange) < 3:
            step = 0.1
        else:
            step = xRange[2]
    x = np.arange(xmin,xmax,step)
    # multiComponent gives the summed 'model' and the components 'grey'. 
    # 'grey' is a List 
    if nu_or_lambda == 'nu':
        model,grey = astFunc.multi_component_grey_body(p2,x,'nu',kappa)
    if nu_or_lambda=='lambda':
        model,grey = astFunc.multi_component_grey_body(p2,x,'nu',kappa)
        y=copy(x)
        modelLambda =copy(model)
        greyLambda = [] 
        for i in range(len(grey)):
            greyLambda += [copy(grey[i])]
        for i in range(len(x)):
            y[i]=(const.c/(x[len(x)-i-1]*1e9))/1e-6
        for i in range(len(model)):
            modelLambda[i]=model[len(model)-i-1]
            for j in range(len(greyLambda)):
                greyLambda[j][i] = grey[j][len(grey[j])-i-1]
        x=y
        model =modelLambda
        grey = greyLambda
    plt.loglog(x,model,ls='-',color=color,label='_nolegend_',lw=0.5,marker='')
    linestyles = [':','-.','-']
    j=0
    for i in grey:
        plt.loglog(x,i,color=color,ls=linestyles[j],lw=0.5,marker='')
        j+=1

def create_figure(pList, data, parted=[1], plotLegend='no', label=['Fluxes'], 
                 color=['black'], marker=['x'], markerSize=6, titleString='', 
                 xLabel='',yLabel='', textString1=None, nu_or_lambda='nu', 
                 fontdict=None, printFitResult=True, fitResultLoc=[50,10],  
                 lineWidth=0.5, kappa='easy', chisq='', xRange='normal', 
                 plotXlabel=None, plotYlabel=None, noXtics=False, 
                 noYtics=False, lineColor='black', ylim=[1e-3,3e1], 
                 xlim=[500,20000]):
    '''
    Plots the total SED of M33. Resamples the older gildas output. 
    input: 
    pList: Multi-Component GreyBody parameters. pList = [[t1,t2,t3,...],
    [N1,N2,N3,...],[beta]] 
    '''
    # ensures the current figure is an empty page
    # generates the Text that shows the fit results for this Plot
    # if pList corresponds to a 2 component grey Body it prints the N1/N2 value
    textString = ''
    for i in range(len(pList[0])):
        textString += (r'T'+str(i+1)+'='+str('%1.1f'%pList[0][i])+' K\nM'+
                      str(i+1)+'='+str("%1.2e"%pList[1][i])+' M$_{\odot}$\n')
    if len(pList[0])==2:
        textString += r'N1/N2 = '+str('%i'%(pList[1][0]/pList[1][1]))+'\n'
    textString += (r'beta = '+str("%1.2f"%pList[2][0])+'\n$\\chi^2$ ='+
                   str("%1.2f"%chisq)+'\n')
    print textString
    # sets the limits of the plot Page
    plotSize = 0.9 # percentace of plotpage larger than the plotted values?
    if nu_or_lambda=='nu':
        xLimNu=[min(data[0])-min(data[0])*plotSize,max(data[0])+max(data[0])]
    if nu_or_lambda == 'lambda':
        newData = []
        for i in data[0]:
            #print i
            newData+=[const.c / (i * 1e9) / 1e-6]
        data[0]=newData
        xLimNu=[min(data[0]) - min(data[0]) * plotSize*2,
                max(data[0]) + max(data[0]) * plotSize]
    # Steering the xtick-Labels if more than one plot is to be connected.	
    axNow= plt.gca()
    plt.setp( axNow.get_xticklabels(), visible=True)
    if noXtics == True:	
        plt.setp( axNow.get_xticklabels(), visible=False)
    if noYtics == True:
        plt.setp( axNow.get_yticklabels(), visible=False) 	
    plt.xlim(xlim[0],xlim[1])
    plt.ylim(ylim[0],ylim[1])
    # reads the Data for Plotting	
    # PLots the model given in pList
    plot_sed(pList,nu_or_lambda,kappa=kappa,xRange=xRange,color=lineColor)
    markersize = 6
    #Plotting the data points
    if len(parted)==1:	
        plt.errorbar(data[0], data[1], yerr=data[2], fmt='o', marker='p', 
                     mfc='None', mew=0.5, mec='#00ffff', ms=markersize, 
                     color='black',lw=lineWidth)
    else:
        for i in range(len(parted)):
            if i == 0:
                plt.errorbar(data[0][0:parted[i]], data[1][0:parted[i]], 
                             yerr=data[2][0:parted[i]], fmt=marker[i], 
                             marker=marker[i], mfc='None', label=label[i], 
                             mew=0.5, mec=color[i], ms=markersize, 
                             color=color[i], lw=lineWidth)
            else:
                plt.errorbar(data[0][parted[i-1]:parted[i]], 
                             data[1][parted[i-1]:parted[i]], 
                             yerr=data[2][parted[i-1]:parted[i]], fmt=marker[i],
                             marker=marker[i], mfc='None', label=label[i],
                             mew=0.5, mec=color[i], ms=markersize, 
                             color=color[i], lw=lineWidth)
    # setting up legend,title, xlabel.
    if plotLegend == True:	
        fontdict={'size':'13'}
        plt.legend(loc='upper right', numpoints=1, fancybox=False, 
                   prop=fontdict, markerscale=1)	
    fontdict={'size':'22'}
    if printFitResult==True:
        fontdict={'size':'12'}
        plt.text(fitResultLoc[0], fitResultLoc[1], s=textString, 
                 fontdict=fontdict, alpha=0.4)
        fontdict={'size':'22'}
    fontdict={'size':'10'}
    if textString1 != None:
        plt.text(5,10,s=textString1, fontdict=fontdict)
    plt.title(titleString)
    if plotXlabel==True:
        plt.xlabel(xLabel)
    if plotYlabel==True:
        plt.ylabel(yLabel)
    plt.axis([xlim[0],xlim[1],ylim[0],ylim[1]])
    axNow= plt.gca()
    plt.setp( axNow.get_xticklabels(), visible=True)
