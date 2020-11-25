/* mex_create_array.h */
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include "mex.h"

double** set_mxDoubleMatrix(const mxArray *pmi){
    mwSize M,N;
    mwIndex j;
    double **pm;
    M = mxGetM(pmi); N = mxGetN(pmi);
    pm = (double **) mxMalloc(N*sizeof(double*));
    pm[0] = mxGetDoubles(pmi);
    for(j=1;j<N;j++){
        pm[j] = pm[j-1]+M;
    }
    return pm;
}

float** set_mxSingleMatrix(const mxArray *pmi){
    mwSize M,N;
    mwIndex j;
    float **pm;
    M = mxGetM(pmi); N = mxGetN(pmi);
    pm = (float **) mxMalloc(N*sizeof(float*));
    pm[0] = mxGetSingles(pmi);
    for(j=1;j<N;j++){
        pm[j] = pm[j-1]+M;
    }
    return pm;
}

bool** set_mxLogicalMatrix(const mxArray *pmi){
    mwSize M,N;
    mwIndex j;
    bool **pm;
    M = mxGetM(pmi); N = mxGetN(pmi);
    pm = (bool **) mxMalloc(N*sizeof(bool*));
    pm[0] = mxGetLogicals(pmi);
    for(j=1;j<N;j++){
        pm[j] = pm[j-1]+M;
    }
    return pm;
}

int8_T** set_mxInt8Matrix(const mxArray *pmi){
    mwSize M,N;
    mwIndex j;
    int8_T **pm;
    M = mxGetM(pmi); N = mxGetN(pmi);
    pm = (int8_T **) mxMalloc(N*sizeof(int8_T*));
    pm[0] = mxGetInt8s(pmi);
    for(j=1;j<N;j++){
        pm[j] = pm[j-1]+M;
    }
    return pm;
}

int16_T** set_mxInt16Matrix(const mxArray *pmi){
    mwSize M,N;
    mwIndex j;
    int16_T **pm;
    M = mxGetM(pmi); N = mxGetN(pmi);
    pm = (int16_T **) mxMalloc(N*sizeof(int32_T*));
    pm[0] = mxGetInt16s(pmi);
    for(j=1;j<N;j++){
        pm[j] = pm[j-1]+M;
    }
    return pm;
}

int32_T** set_mxInt32Matrix(const mxArray *pmi){
    mwSize M,N;
    mwIndex j;
    int32_T **pm;
    M = mxGetM(pmi); N = mxGetN(pmi);
    pm = (int32_T **) mxMalloc(N*sizeof(int32_T*));
    pm[0] = mxGetInt32s(pmi);
    for(j=1;j<N;j++){
        pm[j] = pm[j-1]+M;
    }
    return pm;
}