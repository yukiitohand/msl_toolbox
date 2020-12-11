/* mex_create_array.h */
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include "mex.h"

#ifndef MEX_CREATE_ARRAY_H
#define MEX_CREATE_ARRAY_H

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

/* create a column oriented MxN matrix accessed by ar2d[n][m] */
void createDoubleMatrix(double ***ar2d, double **ar_base, size_t N, size_t M)
{
    size_t ni;
    
    *ar2d = (double**) malloc(sizeof(double*) * N);
    *ar_base = (double*) malloc(sizeof(double) * N * M);
    (*ar2d)[0] = *ar_base;
    for(ni=1;ni<N;ni++){
        (*ar2d)[ni] = (*ar2d)[ni-1] + M;
    }
}

void createSingleMatrix(float ***ar2d, float **ar_base, size_t N, size_t M)
{
    size_t ni;
    
    *ar2d = (float**) malloc(sizeof(float*) * N);
    *ar_base = (float*) malloc(sizeof(float) * N * M);
    (*ar2d)[0] = *ar_base;
    for(ni=1;ni<N;ni++){
        (*ar2d)[ni] = (*ar2d)[ni-1] + M;
    }
}

/* create a column oriented MxN matrix accessed by ar2d[n][m] */
void createInt32Matrix(int32_T ***ar2d, int32_T **ar_base, size_t N, size_t M)
{
    size_t ni;
    
    *ar2d = (int32_T**) malloc(sizeof(int32_T*) * N);
    *ar_base = (int32_T*) malloc(sizeof(int32_T) * N * M);
    (*ar2d)[0] = &(*ar_base)[0];
    for(ni=1;ni<N;ni++){
        (*ar2d)[ni] = (*ar2d)[ni-1] + M;
    }
}

/* Int32 Pointer matrix */
void createInt32PMatrix(int32_T ****ar2d, int32_T ***ar_base, size_t N, size_t M)
{
    size_t ni;
    
    *ar2d = (int32_T***) malloc(sizeof(int32_T**) * N);
    *ar_base = (int32_T**) malloc(sizeof(int32_T*) * N * M);
    (*ar2d)[0] = &(*ar_base)[0];
    for(ni=1;ni<N;ni++){
        (*ar2d)[ni] = (*ar2d)[ni-1] + M;
    }
}

void createDoublePMatrix(double ****ar2d, double ***ar_base, size_t N, size_t M)
{
    size_t ni;
    
    *ar2d = (double***) malloc(sizeof(double**) * N);
    *ar_base = (double**) malloc(sizeof(double*) * N * M);
    (*ar2d)[0] = &(*ar_base)[0];
    for(ni=1;ni<N;ni++){
        (*ar2d)[ni] = (*ar2d)[ni-1] + M;
    }
}

#endif