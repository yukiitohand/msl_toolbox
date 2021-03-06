/* =====================================================================
 * get_msldemtUFOVmask_wmsldemc_mex.c
 * Evaluate any pixels in MSLDEM image whether or not they exist in the 
 * 
 * INPUTS:
 * 0 msldemc_img           Double array [L_demc x S_demc]
 * 1 msldemc_northing      Double array [L_demc]
 * 2 msldemc_easting       Double array [S_demc]
 * 3 msldemc_imFOVmask     int8_T [L_demc x S_demc]
 * 4 S_im                  int
 * 5 L_im                  int
 * 6 cahv_mdl              CAHV_MODEL
 * 7 msldemt_img           Double array [L_demt x S_demt]
 * 8 msldemt_northing      Double array [L_demt]
 * 9 msldemt_easting       Double array [S_demt]
 * 10 msldemt_imFOVmask    int8_T [L_demt x S_demt]
 *
 * The origin of msldemc_img, msldemc_northing, and msldemc_easting is 
 * cam_C.
 * 
 * 
 * OUTPUTS:
 * 0  msldemt_inImage     [L_demc x S_demc]  Boolean
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
#include "cahvor.h"
#include "lib_proj_mastcamMSLDEM.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double **msldemc_img;
    double *msldemc_northing;
    double *msldemc_easting;
    double **msldemt_img;
    double *msldemt_northing;
    double *msldemt_easting;
    CAHV_MODEL cahv_mdl;
    int8_T **msldemc_imFOVmask;
    mwSize S_im,L_im;
    
    int8_T **msldemt_imFOVmask;
    
    mwIndex si,li;

    
    mwSize S_demc, L_demc, S_demt, L_demt;

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
    
    /* INPUT 0 msldemc_img */
    msldemc_img = set_mxDoubleMatrix(prhs[0]);
    L_demc = mxGetM(prhs[0]);
    S_demc = mxGetN(prhs[0]);
    
    /* INPUT 1/2 msldem northing easting */
    msldemc_northing = mxGetDoubles(prhs[1]);
    msldemc_easting = mxGetDoubles(prhs[2]);
    
    /* INPUT 3 msldemc imFOV */
    msldemc_imFOVmask = set_mxInt8Matrix(prhs[3]);
    
    /* INPUT 4/5 image S_im, L_im */
    S_im = (mwSize) mxGetScalar(prhs[4]);
    L_im = (mwSize) mxGetScalar(prhs[5]);
    
    /* INPUT 6 CAHV model */
    cahv_mdl = mxGet_CAHV_MODEL(prhs[6]);
    
    /* INPUT 7 msldemt_img */
    msldemt_img = set_mxDoubleMatrix(prhs[7]);
    L_demt = mxGetM(prhs[7]);
    S_demt = mxGetN(prhs[7]);
    
    /* INPUT 8/9 msldem northing easting */
    msldemt_northing = mxGetDoubles(prhs[8]);
    msldemt_easting = mxGetDoubles(prhs[9]);
    
    /* INPUT 10 msldemt imFOV */
    plhs[0] = mxDuplicateArray(prhs[10]);
    msldemt_imFOVmask = set_mxInt8Matrix(plhs[0]);
    
    
    /* OUTPUT 0 Unobstructed FOV mask */
    // plhs[0] = mxCreateNumericMatrix(L_demt,S_demt,mxINT8_CLASS,mxREAL);
    // msldemt_inImage = set_mxInt8Matrix(plhs[0]);
    
    // Initialize matrices
    // for(si=0;si<S_demt;si++){
    //    for(li=0;li<L_demt;li++){
    //        msldemt_inImage[si][li] = msldemc_imFOVmask[si][li];
    //    }
    // }
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    mask_obstructed_pts_in_msldemt_using_msldemc(
            (int32_T) S_demc, (int32_T) L_demc, 
            msldemc_northing, msldemc_easting, msldemc_img, msldemc_imFOVmask,
            (int32_T) S_demt, (int32_T) L_demt, 
            msldemt_northing, msldemt_easting, msldemt_img, msldemt_imFOVmask,
            (int32_T) S_im, (int32_T) L_im, cahv_mdl);
    
    /* free memories */
    mxFree(msldemc_img);
    mxFree(msldemt_img);
    mxFree(msldemc_imFOVmask);
    mxFree(msldemt_imFOVmask);
    
}