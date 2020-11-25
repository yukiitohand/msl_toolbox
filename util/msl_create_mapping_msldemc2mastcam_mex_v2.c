/* =====================================================================
 * msl_create_mapping_msldemc2mastcam_mex_v2.c
 * Create a mapping cell array from DEM to mastcam image
 * 
 * INPUTS:
 * 0 msldemc_imUFOVnnx         [Ldem x Sdem] int16 matrix 
 * 1 msldemc_imUFOVnny         [Ldem x Sdem] int16 matrix 
 * 2 mastcam_nnx               [L_im x S_im] int32 matrix
 * 3 mastcam_nny               [L_im x S_im] int32 matrix 
 * 
 * (The index [0,0] here corresponds to the index [1,1] in MATLAB.)
 * 
 * OUTPUTS:
 * 0 mapidx_msldemc2mastcam [Ldem x Sdem] int32 array, whose elements 
 *    stores for indicies associated to those of "mapcell_mastcam2msldemc"
 *
 * 1 mapcell_msldemc2mastcam [* x 1] cell array whose elements stores 
 * [ * x 2] array of corresponding pixel indicies.
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
void create_mapper_msldemc2mastcam(
        int16_T **msldemc_imUFOVnnx, int16_T **msldemc_imUFOVnny,
        int32_T **mastcam_nnx, int32_T **mastcam_nny,
        int32_T L_dem, int32_T S_dem, int32_T L_im, int32_T S_im,
        int32_T **mapidx_msldemc2mastcam, mxArray *mapcell_msldemc2mastcam)
{
    
    
    mwSize *mapper_countar_base;
    mwSize **mapper_countar;
    // int32_T **mapper_countar_base_c, ***mapper_countar_base_c;
    int16_T **mapar_base;
    int16_T ***mapar;
    int32_T si,li;
    int32_T imx,imy;
    mwSize mwsi,mwli,n2;
    // mwSize imLen;
    mwSize mwL_dem, mwS_dem;
    mwSize mwL_im, mwS_im;
    mwIndex imidx[2];
    mwIndex idx1d;
    
    
    mwL_dem = (mwSize) L_dem; mwS_dem = (mwSize) S_dem;
    mwL_im = (mwSize) L_im; mwS_im = (mwSize) S_im;
    
    // imLen = (size_t) L_im * (size_t) S_im;
    
    // printf("%d\n",S_im);
    
    /* First counting the number of elements to be stored in each  */
    /* initialize countar array */
    mapper_countar_base = (mwSize*) mxMalloc(sizeof(mwSize) * (size_t) mwL_dem * (size_t) mwS_dem);
    mapper_countar = (mwSize**) mxMalloc(sizeof(mwSize*) * (size_t) mwS_dem);
    mapper_countar[0] = &mapper_countar_base[0];
    for(mwsi=1;mwsi<mwS_dem;mwsi++)
        mapper_countar[mwsi] = mapper_countar[mwsi-1] + mwL_dem;
    for(mwsi=0;mwsi<mwS_dem;mwsi++){
        for(mwli=0;mwli<mwL_dem;mwli++){
            if(msldemc_imUFOVnnx[mwsi][mwli] > -1){
                mapper_countar[mwsi][mwli] = 1;
            } else {
                mapper_countar[mwsi][mwli] = 0;
            }
        }
    }
    // printf("%d\n",S_im);
    
    /* count the number of pixels stored for each msldemc image pixels */
    for(si=0;si<S_im;si++){
        for(li=0;li<L_im;li++){
            if(mastcam_nnx[si][li]>-1){
                imx = mastcam_nnx[si][li];
                imy = mastcam_nny[si][li];
                if((msldemc_imUFOVnnx[imx][imy] != si) || (msldemc_imUFOVnny[imx][imy] != li)){
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
    mapar_base = (int16_T**) mxMalloc(sizeof(int16_T*) * (size_t) L_dem * (size_t) S_dem);
    mapar = (int16_T***) mxMalloc(sizeof(int16_T**) * (size_t) S_dem);
    mapar[0] = &mapar_base[0];
    for(si=1;si<S_dem;si++)
        mapar[si] = mapar[si-1] + L_dem;
    for(si=0;si<S_dem;si++){
        for(li=0;li<L_dem;li++){
            mapar[si][li] = NULL;
        }
    }
    // printf("%d\n",S_im);
    
    for(mwsi=0;mwsi<mwS_dem;mwsi++){
        for(mwli=0;mwli<mwL_dem;mwli++){
            if(mapper_countar[mwsi][mwli]>0){
                idx1d = (mwIndex) mapidx_msldemc2mastcam[mwsi][mwli];
                // printf("%d,%d\n",idx1d,mapper_countar[mwsi][mwli]);
                mxSetCell(mapcell_msldemc2mastcam, idx1d,
                        mxCreateNumericMatrix(2,mapper_countar[mwsi][mwli],mxINT16_CLASS,mxREAL));
                mapar[mwsi][mwli] = mxGetInt16s(mxGetCell(mapcell_msldemc2mastcam,idx1d));
            }
        }
    }
    // printf("%d\n",S_im);
    
    /* initialize again the counter array for counting current indicies */
    for(mwsi=0;mwsi<mwS_dem;mwsi++){
        for(mwli=0;mwli<mwL_dem;mwli++){
            mapper_countar[mwsi][mwli] = 0;
        }
    }
    // printf("%d\n",S_im);
    // printf("%d\n",mapar[1231][408][17]);
    
    /* first fill out the associated pixel ids in msldem_nn */
    for(si=0;si<S_dem;si++){
        for(li=0;li<L_dem;li++){
            if(msldemc_imUFOVnnx[si][li] > -1){
                n2 = mapper_countar[si][li]*2;
                mapar[si][li][n2]   = msldemc_imUFOVnnx[si][li]+1; // for matlab indexing
                mapar[si][li][n2+1] = msldemc_imUFOVnny[si][li]+1; // for matlab indexing
                mapper_countar[si][li]++;
            }
        }
    }
    
    
    for(si=0;si<S_im;si++){
        // printf("si=%d\n",si);
        for(li=0;li<L_im;li++){
            // printf("si=%d,li=%d\n",si,li);
            imx = mastcam_nnx[si][li];
            imy = mastcam_nny[si][li];
            if(imx>-1){
                if((msldemc_imUFOVnnx[imx][imy] != si) || (msldemc_imUFOVnny[imx][imy] != li)){
                    n2 = mapper_countar[imx][imy]*2;
                    mapar[imx][imy][n2]   = (int16_T) si+1; // for matlab indexing
                    mapar[imx][imy][n2+1] = (int16_T) li+1; // for matlab indexing
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
    int32_T **mapidx_msldemc2mastcam;
    mwSize L_dem,S_dem;
    mwSize L_im,S_im;
    mxArray *mapcell_msldemc2mastcam;
    mwSize cell_len[1];
    
    mwIndex si,li;
    int32_T cur_idx;

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
    plhs[0] = mxCreateNumericMatrix(L_dem,S_dem,mxINT32_CLASS,mxREAL);
    mapidx_msldemc2mastcam = set_mxInt32Matrix(plhs[0]);
    
    // Initialize matrices
    cur_idx = 0;
    for(si=0;si<S_dem;si++){
        for(li=0;li<L_dem;li++){
            if(msldemc_imUFOVnnx[si][li]>-1){
                mapidx_msldemc2mastcam[si][li] = cur_idx;
                cur_idx++;
            } else {
                mapidx_msldemc2mastcam[si][li] = -1;
            }
        }
    }
    
    cell_len[0] = (mwSize) cur_idx;
    mapcell_msldemc2mastcam = mxCreateCellArray(1,cell_len);
    plhs[1] = mapcell_msldemc2mastcam;
    
    
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    create_mapper_msldemc2mastcam(msldemc_imUFOVnnx,msldemc_imUFOVnny,
            mastcam_nnx,mastcam_nny,
            (int32_T) L_dem, (int32_T) S_dem, (int32_T) L_im, (int32_T) S_im, 
            mapidx_msldemc2mastcam,mapcell_msldemc2mastcam);
    
    /* free memories */
    mxFree(msldemc_imUFOVnnx);
    mxFree(msldemc_imUFOVnny);
    mxFree(mastcam_nnx);
    mxFree(mastcam_nny);
    mxFree(mapidx_msldemc2mastcam);
    
}
