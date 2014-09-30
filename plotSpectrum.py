from pylab import plot, show, title, xlabel, ylabel, subplot, figure, hold, legend,xlim,savefig,clf,ylim
from scipy import fft
from numpy import  nonzero, shape, correlate, mean, std,linspace,savetxt, copy,array,arange,floor

def plotSpectrum(data,timeR,Fourier=1,Corr=1,index = array([5,6]),figname='a.pdf',ind=1):

    [ind1,ind2] = index
    k = min([60*5,len(data[:,0])/3])

    y = data[k:,ind1]
    x = data[k:,ind2]
    time = timeR[k:]

    if Fourier:
        figure(1)
        subplot(2,1,1)
        hold(True)
        plot(time,y)
        plot(time,x)
        xlabel('Time')
        ylabel('Amplitude')
        subplot(2,1,2)
        Y = fft(y)
        Y = Y[1:]
        n= len(Y)

        #Crazy math in here
        powerY = abs(Y[arange(0,floor(n/2),dtype='int')])**2
        nyquist = 1./2.
        freq = arange(0,n/2,dtype='float')/(n/2)*nyquist
        period=1./freq
        # a is the number of points per time step (Depends on the interpolation in process - interp)
        a = 1
        plot(period/a,powerY)
        ylabel('Power')
        xlabel('Period (Min/oscillation)')
        #savefig(''.join([figname,".png"]), bbox_inches=0 ,dpi=100)


    if Corr:
        modeC = "same"
        x = (x - mean(x))/std(x)
        y =  (y - mean(y))/std(y)

        timeInt = time[1]-time[0]
        numPoints = len(time)

        fig2 = figure(2)
        fig2.set_size_inches(6*0.9,3.7*0.9*0.906)
        fig2.set_facecolor('white')

        ax2 = fig2.add_subplot(1,1,ind)
        n= correlate(y,y,modeC)
        d= correlate(y,x,modeC)
        lagR = (nonzero(n == max(n))[0] - nonzero(d == max(d))[0])*timeR[1]
        title("CC Act-DS %f Max corr: %f" %(lagR, max(d)/numPoints))
        ax2.plot(linspace(len(x)/2*timeInt,-len(x)/2*timeInt,len(x)),d/numPoints)
        xlim([-180,180])
        #ylim([-0.25,0.75])

        #savefig(''.join([figname,".png"]), bbox_inches=0 ,dpi=100)

    show()