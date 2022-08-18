/* =====================================================================
 * msl_create_mapping_mastcam2msldemc_mex_v2.c
 * Create a mapping cell array from mastcam image to DEM.
 * 
 * INPUTS:
 * 0 msldemc_imUFOVnnx         [Ldem x Sdem] int16 matrix 
 * 1 msldemc_imUFOVnny         [Ldem x Sdem] int16 matrix 
 * 2 mastcam_msldemcUFOV_nnx   [L_im x S_im] int32 matrix
 * 3 mastcam_msldemcUFOV_nny   [L_im x S_im] int32 matrix 
 * 
 * (The index [0,0] here corresponds to the index [1,1] in MATLAB.)
 * 
 * OUTPUTS:
 * 0 map_mastcam2msldemc [Lim x Sim] int32 cell array, whose elements 
 *    stores [ * x 2] array of corresponding pixel indicies.
 * 
 *
 *
 * This is a MEX file for MATLAB.
 * ===================================================================== */
#include <stdint.h>
#include "io64.h"
#include "mex.h"
#include "math.h"
#include "matrix.h"
#include <string.h>
#include <stdio.h>

#include <stdlib.h>
#include "envi.h"
#include "mex_create_array.h"

/* main computation routine */
void create_mapper_mastcam2msldemc(
        int16_T **msldemc_imUFOVnnx, int16_T **msldemc_imUFOVnny,
        int32_T **mastcam_nnx, int32_T **mastcam_nny,
        int32_T L_dem, int32_T S_dem,
        int32_T L_im, int32_T S_im, mxArray *map_mastcam2msldemc)
{
    
    
    mwSize *mapper_countar_base;
    mwSize **mapper_countar;
    // int32_T **mapper_countar_base_c, ***mapper_countar_base_c;
    int32_T **mapar_base;
    int32_T ***mapar;
    int32_T si,li;
    int16_T imx,imy;
    mwSize mwsi,mwli,n2;
    // mwSize imLen;
    mwSize mwL_dem, mwS_dem;
    mwSize mwL_im, mwS_im;
    mwIndex imidx[2];
    mwIndex imidx1d;
    
    
    mwL_dem = (mwSize) L_dem; mwS_dem = (mwSize) S_dem;
    mwL_im = (mwSize) L_im; mwS_im = (mwSize) S_im;
    
    // imLen = (size_t) L_im * (size_t) S_im;
    
    // printf("%d\n",S_im);
    
    /* First counting the number of elements to be stored in each  */
    /* initialize countar array */
    mapper_countar_base = (mwSize*) mxMalloc(sizeof(mwSize) * (size_t) mwL_im * (size_t) mwS_im);
    mapper_countar = (mwSize**) mxMalloc(sizeof(mwSize*) * (size_t) mwS_im);
    mapper_countar[0] = &mapper_countar_base[0];
    for(mwsi=1;mwsi<mwS_im;mwsi++)
        mapper_countar[mwsi] = mapper_countar[mwsi-1] + mwL_im;
    for(mwsi=0;mwsi<mwS_im;mwsi++){
        for(mwli=0;mwli<mwL_im;mwli++){
            if(mastcam_nnx[mwsi][mwli] > -1){
                mapper_countar[mwsi][mwli] = 1;
            } else {
                mapper_countar[mwsi][mwli] = 0;
            }
        }
    }
    // printf("%d\n",S_im);
    
    /* count the number of pixels stored for each mastcam image pixels */
    for(si=0;si<S_dem;si++){
        for(li=0;li<L_dem;li++){
            if(msldemc_imUFOVnnx[si][li]>-1){
                imx = msldemc_imUFOVnnx[si][li];
                imy = msldemc_imUFOVnny[si][li];
                if((mastcam_nnx[imx][imy] != si) || (mastcam_nny[imx][imy] != li)){
                    mapper_countar[imx][imy]++;
                }
            }
        }
    }
    // printf("%d\n",S_im);
    
    /* */
    //mapper_countar_base_c = (mwSize*) mxMalloc(sizeof(mwSize) * (size_t) L_im * (size_t) S_im);
    //mapper_countar_c = (mwSize**) mxMalloc(sizeof(mwSize*) * (size_t) S_im);
    //mapper_countar_c[0] = &mapper_countar_base_c[0];
    //for(si=1;si<S_im;si++)
    //    mapper_countar_c[si] = mapper_countar_c[si-1] + L_im;
    //for(si=0;si<S_im;si++){
    //    for(li=0;li<L_im;li++){
    //        mapper_countar_c[si][li] = 0;
    //    }
    //}
    
    
    /* Prepare a dummy matrix for easy access to the output cell array */
    mapar_base = (int32_T**) mxMalloc(sizeof(int32_T*) * (size_t) L_im * (size_t) S_im);
    mapar = (int32_T***) mxMalloc(sizeof(int32_T**) * (size_t) S_im);
    mapar[0] = &mapar_base[0];
    for(si=1;si<S_im;si++)
        mapar[si] = mapar[si-1] + L_im;
    for(si=0;si<S_im;si++){
        for(li=0;li<L_im;li++){
            mapar[si][li] = NULL;
        }
    }
    // printf("%d\n",S_im);
    
    for(mwsi=0;mwsi<mwS_im;mwsi++){
        imidx[1] = mwsi;
        for(mwli=0;mwli<mwL_im;mwli++){
            if(mapper_countar[mwsi][mwli]>0){
                imidx[0] = mwli;
                imidx1d = mxCalcSingleSubscript(map_mastcam2msldemc,2,imidx);
                // printf("%d,%d,%d,%d\n",imidx[0],imidx[1],imidx1d,mapper_countar[mwsi][mwli]);
                mxSetCell(map_mastcam2msldemc,imidx1d,
                        mxCreateNumericMatrix(2,mapper_countar[mwsi][mwli],mxINT32_CLASS,mxREAL));
                mapar[mwsi][mwli] = mxGetInt32s(mxGetCell(map_mastcam2msldemc,imidx1d));
            }
        }
    }
    // printf("%d\n",S_im);
    
    /* initialize again the counter array for counting current indicies */
    for(mwsi=0;mwsi<mwS_im;mwsi++){
        for(mwli=0;mwli<mwL_im;mwli++){
            mapper_countar[mwsi][mwli] = 0;
        }
    }
    // printf("%d\n",S_im);
    // printf("%d\n",mapar[1231][408][17]);
    
    /* first fill out the associated pixel ids in mastcam_nn */
    for(si=0;si<S_im;si++){
        for(li=0;li<L_im;li++){
            if(mastcam_nnx[si][li] > -1){
                n2 = mapper_countar[si][li]*2;
                mapar[si][li][n2]   = mastcam_nnx[si][li]+1; // for matlab indexing
                mapar[si][li][n2+1] = mastcam_nny[si][li]+1; // for matlab indexing
                mapper_countar[si][li]++;
            }
        }
    }
    
    
    for(si=0;si<S_dem;si++){
        // printf("si=%d\n",si);
        for(li=0;li<L_dem;li++){
            // printf("si=%d,li=%d\n",si,li);
            imx = msldemc_imUFOVnnx[si][li];
            imy = msldemc_imUFOVnny[si][li];
            if(imx>-1){
                if((mastcam_nnx[imx][imy] != si) || (mastcam_nny[imx][imy] != li)){
                    n2 = mapper_countar[imx][imy]*2;
                    mapar[imx][imy][n2]   = si+1; // for matlab indexing
                    mapar[imx][imy][n2+1] = li+1; // for matlab indexing
                    mapper_countar[imx][imy]++;
                }
            }
        }
    }
    
    // printf("%d\n",S_im);
    mxFree(mapar);
    mxFree(mapar_base);
    mxFree(mapper_countar);
    mxFree(mapper_countar_base);
    
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int16_T **msldemc_imUFOVnnx, **msldemc_imUFOVnny;
    int32_T **mastcam_nnx, **mastcam_nny;
    mwSize L_dem,S_dem;
    mwSize L_im,S_im;
    mxArray *map_mastcam2msldemc;
    
    mwIndex si,li;

    /* -----------------------------------------------------------------
     * CHECK PROPER NUMBER OF INPUTS AND OUTPUTS
     * ----------------------------------------------------------------- */
    /* if(nrhs!=11) {
        mexErrMsgIdAndTxt("proj_mastcam2MSLDEM_v4_mex:nrhs","Eleven inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("proj_mastcam2MSLDEM_v4_mex:nlhs","Five outputs required.");
    }
    */
    /* make sure the first input argument is scalar */
    /*
    if( !mxIsChar(prhs[0]) ) {
        mexErrMsgIdAndTxt("proj_mastcam2MSLDEM_v4_mex:notChar","Input 0 needs to be a character vector.");
    }
    */
    /* -----------------------------------------------------------------
     * I/O SETUPs
     * ----------------------------------------------------------------- */
    
    /* INPUT 0,1 msldemc_imUFOVnnx,msldemc_imUFOVnny */
    msldemc_imUFOVnnx = set_mxInt16Matrix(prhs[0]);
    msldemc_imUFOVnny = set_mxInt16Matrix(prhs[1]);
    L_dem = mxGetM(prhs[0]);
    S_dem = mxGetN(prhs[1]);
    
    /* INPUT 2,3 msldemc projection parameters */
    mastcam_nnx = set_mxInt32Matrix(prhs[2]);
    mastcam_nny = set_mxInt32Matrix(prhs[3]);
    L_im   = mxGetM(prhs[2]);
    S_im   = mxGetN(prhs[2]);
    
    

    /* OUTPUT 0 crism_imUFOVmask */
    map_mastcam2msldemc = mxCreateCellMatrix(L_im,S_im);
    plhs[0] = map_mastcam2msldemc;
    
    // Initialize matrices
    
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    create_mapper_mastcam2msldemc(msldemc_imUFOVnnx,msldemc_imUFOVnny,
            mastcam_nnx,mastcam_nny,
            (int32_T) L_dem, (int32_T) S_dem, 
            (int32_T) L_im, (int32_T) S_im, map_mastcam2msldemc);
    
    /* free memories */
    mxFree(msldemc_imUFOVnnx);
    mxFree(msldemc_imUFOVnny);
    mxFree(mastcam_nnx);
    mxFree(mastcam_nny);
    
}
