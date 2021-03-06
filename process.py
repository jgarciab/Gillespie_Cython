from numpy import *
import time as tm
from scipy.stats import mstats, skew, probplot
import numpy as np

import init
import loopGillespie
import loopInit
import plotMod
import plotSpectrum


# Use imports from another folder
import sys
sys.path.append( '/home/j/Dropbox/Python/marModel8/Peaks')


def interp(cNT,timeReaction):
    dataInterp = zeros([timeReaction,shape(cNT)[1]-1])
    realTime = cNT[:,0]
    time = linspace(np.min(realTime),np.max(realTime),timeReaction)

    for i in range(shape(cNT)[1]-1):
        dataInterp[:,i] = np.interp(time, cNT[:,0], cNT[:,(i+1)])


    return [time, dataInterp]

def statData(time, dataInterp, numberOfReactants,numberOfReactions,printOp = True):
    quantilesArray = zeros([10,numberOfReactants])
    z = shape(dataInterp)[0]
    z = 2*int(z/3)

    #Quantiles, mean, stdv, skew
    for i in range(numberOfReactants):
        quantilesArray[0:4,i] = mstats.mquantiles(dataInterp[z:,i],prob=[0.01 , 0.5 , 0.75, 0.99])
        quantilesArray[4,i] = mean(dataInterp[z:,i])
        quantilesArray[5,i] = std(dataInterp[z:,i])
        quantilesArray[6,i] = skew(dataInterp[z:,i])
        quantilesArray[7,i] = quantilesArray[5,i]/quantilesArray[4,i] # Noise: Unambiguous
        # Noise strength: Get trends that may be obscured by the 1/sqrt(N) scaling of noise due to finite effects
        quantilesArray[8,i] = quantilesArray[5,i]**2/quantilesArray[4,i]


    print max(dataInterp[z:,-2])
    quantilesArray[-1,:] = range(numberOfReactants)

    quantilesArrayLog = zeros([10,numberOfReactants])

    if printOp == True:
        set_printoptions(precision=3,suppress=True,linewidth=150)

        print
        print("       5%        50%        75%        95%        Mean       STD        Skewness   Noise      Noise Str  Index")
        print(transpose(quantilesArray))



    return [quantilesArray,quantilesArrayLog]


def progress(timeReaction,normBool=True,plotBool=True,printOp = True):

    tic = tm.clock()
    # Get the mathematical model, initial conditions, etc
    print("Step 1: InitInfo",(tm.clock()-tic))
    [M, N, input, output, numberOfReactions,numberOfReactants] = init.initAndProcess()


    print("Step 2: InitArrays",(tm.clock()-tic))
    [updateNmatrix, vecMol2, vecN] = loopInit.loopInfo(numberOfReactions,input,output)

    # Gillespie algorithm
    print("Step 3: Loop",(tm.clock()-tic))

    #[N, sumMar5,std5] = loopGillespie.loopProgress(numberOfReactions, numberOfReactants, N, vecN,updateNmatrix, vecMol2,M,20*60)
    [N, sumMar5,std5] = loopGillespie.loopProgress(numberOfReactions, numberOfReactants, N, vecN,updateNmatrix, vecMol2,M,timeReaction)
    fname = "data2ssBP.dat"

    figName = ''

    #Plot
    print("Step 4: Plots",(tm.clock()-tic))
    cNT = plotMod.plotMod(fname,numberOfReactants,plotX = plotBool,saveFig = 0,index = [4,8],figname = figName,norm=normBool)

    #Statistics
    print("Step 5: Interpolation",(tm.clock()-tic))
    [time, dataInterp] = interp(cNT,timeReaction)

    print("Step 6: Some Statistics",(tm.clock()-tic))
    [quantilesArray,quantilesArrayLog] = statData(time, dataInterp, numberOfReactants,numberOfReactions,printOp = printOp)

    print("Step 7: Fourier transform and autocorr",(tm.clock()-tic))
    #plotSpectrum.plotSpectrum(dataInterp,time,Fourier = 0, Corr = 1,index = [-1,-2])

    return [quantilesArray,dataInterp]




def main():
    progress(300000,normBool=False,plotBool=True,printOp=True)


if __name__=="__main__": main()
