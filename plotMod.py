from pylab import *
from numpy import *
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('legend',**{'fancybox':'True'})
rc('legend',**{'fancybox':'True'})
rc('text', usetex=True)



def plotMod(fname,numberOfReactants,plotX = 0,saveFig = 0,index = [8,9], figname = "1.png",ind=1,norm=True):
    """
    Convert the data file from the Cython file to an numpy array, and plot it
    input:
        fname: data file name
        numberOfReactants
        plotX (Default: No). plot the data?

    output
        cNT: Array with the data
    """


    f = open(fname)
    cNT2 = np.fromfile(f, dtype=([('time','f8'),('mols',str(numberOfReactants)+'i4')]), count=-1, sep='')
    close(fname)
    cNT = np.zeros([len(cNT2),numberOfReactants+1])
    cNT[:,0] = cNT2['time']
    cNT[:,1:] = cNT2['mols']


    if plotX:
        hold(True)
        clf()
        fig = figure(1)
        fig.set_size_inches(6*0.9,3.7*0.9*0.906)
        fig.set_facecolor('white')
        ax = fig.add_subplot(1,1,ind)

        for ind in index:
            A = cNT[:,ind]
            if norm: A = (A-mean(A))/std(A)
            ax.plot(cNT[:,0],A,)
        xlim([50,3000])
        show()


    if saveFig:
        savefig(''.join([fname,".pdf"]), bbox_inches=0 ,dpi=100)



    return cNT

