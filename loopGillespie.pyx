# cython: boundscheck=False
# cython: wraparound=False

import numpy as np
cimport numpy as np
import time

cdef extern from "math.h":
    double log(float theta)
    double sin(double x)
    int abs (int x)
cdef extern from "stdlib.h":
    double rand()
    void srand(unsigned int seed )
    double RAND_MAX
    int round(double x)

from libc.stdio cimport *
from libc.limits cimport *

DTYPE = np.int
ctypedef np.int_t DTYPE_t

DTYPE2 = np.float
ctypedef np.float_t DTYPE2_t

#@profile
def loopProgress(int numberOfReactions, int numberOfReactants, np.ndarray[DTYPE_t, ndim=1] N, list vecN, np.ndarray[DTYPE_t, ndim=2] updateNmatrix, list vecMol2, np.ndarray[DTYPE2_t, ndim=1] M2,int timeReaction):
    cdef double h[1000], a[1000], M[1000]
    cdef double Lsum, Rsum
    cdef int reactionN, sumA, Na ,N_cython[1000], updNm[1000][1000], vecMol2_cython[1000][1000], vecMol2_index[1000],vecN_cython[1000][2], vecN_index[1000]
    cdef double rnum1, rnum2, a0, timeR, t
    cdef double stdA, partialA,sumMarA
    cdef tuple numberReactInvolved
    cdef int i,j,k,ui
    cdef long count


    cdef char fname[80], fline[80]
    cdef char  fline2[80]

    cdef FILE *fp
    cdef FILE *fp2
    sprintf(fname, "data2ssBP.dat")
    fp = fopen(fname, "w+")

    t = 0
    count = 0
    a0 = 0

    stdA = 0
    partialA = 0
    sumMarA = 0

    for  i from 0 <= i < numberOfReactions:
        M[i] = M2[i]

    srand(int(np.random.rand()*1000))

    for  i from 0 <= i < numberOfReactions:
        for j from 0 <= j < np.size(vecMol2[i]):
            vecMol2_cython[i][j] = vecMol2[i][0][j]
        for k from 0 <= k < np.size(vecN[i]):
            vecN_cython[i][k]=vecN[i][0][k]

        vecMol2_index[i] = np.size(vecMol2[i])
        vecN_index[i] = np.size(vecN[i])

    for  i from 0 <= i < numberOfReactants:
        N_cython[i] = N[i]
        for  j from 0 <= j < numberOfReactions:
            updNm[j][i] = updateNmatrix[j,i]

    for  i from 0 <= i < numberOfReactions:
        Na = 1
        for  uj from 0 <= uj < vecN_index[i]:
            j = vecN_cython[i][uj]
            Na = Na * N_cython[j]

        h[i] = Na

        a[i] = h[i]*M[i]
        a0 += a[i]


    # Loop
    while(t<timeReaction):
        rnum1 = (rand()+1)/(RAND_MAX + 1)
        rnum2 = (rand()+1)/(RAND_MAX + 1)


         #Time to the next reaction
        timeR = -log((rnum1))/a0
        t = t + timeR
        if(count%100 == 0):
            fwrite(&t,sizeof(double),1,fp)
            fwrite(N_cython,sizeof(int),numberOfReactants,fp)



        reactionN = 0
        Rsum = a[0]
        while(rnum2*a0 > Rsum ):
            reactionN += 1
            Rsum += a[reactionN]


        # Update N_cython
        for   i from 0 <= i < numberOfReactants:
            N_cython[i] -= updNm[reactionN][i]
            if N_cython[i] < 0:
                print reactionN, rnum1,rnum2, a0, N,N_cython[i],i,updNm[reactionN][i]
                return [np.NAN, np.NAN, np.NAN]


        #Reactions that changes (i) when reactionN occurs
        for  ui from 0 <= ui < vecMol2_index[reactionN]:
            i = vecMol2_cython[reactionN][ui]

            #Indexes of the molecules involved
            Na = 1
            for  uj from 0 <= uj < vecN_index[i]:
                j = vecN_cython[i][uj]
                Na = Na * N_cython[j]

            h[i] = Na
            a0 -= a[i]
            a[i] = h[i]*M[i]
            a0 += a[i]


        count += 1

        #This calculate the mean and std
        #sumMarA = sumMarA  + N_cython[3]*timeR
        #if t >1200:
        #    partialA = sumMarA/t
        #    stdA = stdA +  ((N_cython[3]-partialA)*(N_cython[3]-partialA))*timeR


    fclose(fp)


    for i in range(numberOfReactants):
        N[i] = N_cython[i]

    #sumMarA = partialA
    #stdA = stdA/(t-1200)

    return [N, sumMarA,stdA]
