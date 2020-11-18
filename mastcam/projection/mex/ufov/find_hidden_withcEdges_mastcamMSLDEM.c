/* =====================================================================
 * find_hidden_withcEdges_mastcamMSLDEM.c
 * Evaluate any pixels in MSLDEM image whether or not they exist in the 
 * 
 * INPUTS:
 * 0 msldemcedge_inImage   int8 array [(L_demc-1) x (S_demc-1)]
 * 1 msldemc_imFOVmask     int8 array [L_demc x S_demc] 
 *
 * OUTPUTS:
 * 0  msldemc_inImage     [L_demc x S_demc]  Boolean
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
void find_hidden(int32_T msldemc_samples, int32_T msldemc_lines,
        int8_T **msldemcedge_inImage, int8_T **msldemc_imFOVmask, 
        int8_T **msldemc_inImage)
{
    int32_T c,l;
    int32_T S_demcm1;
    
    //L_demcm1 = msldemc_lines - 1;
    S_demcm1 = msldemc_samples - 1;
    //L_demcm1_dbl = (double) L_demcm1_dbl;
    //S_demcm1_dbl = (double) S_demcm1_dbl;
    // printf("%d\n",msldemc_samples);
    
    for(l=0;l<msldemc_lines;l++){
        //printf("l=%d\n",l);
        // decide the first and last indexes to be assessed.
        for(c=0;c<S_demcm1;c++){
            if(msldemcedge_inImage[c][l]>0){
                if(msldemc_imFOVmask[c][l]>0)
                    msldemc_inImage[c][l] = 1;
                if(msldemc_imFOVmask[c+1][l]>0)
                    msldemc_inImage[c+1][l] = 1;
            }
        }
    }
    
    
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{  
    int8_T **msldemcedge_inImage;
    int8_T **msldemc_imFOVmask;
    int8_T **msldemc_inImage;
    

    
    mwSize S_demc, L_demc;
    mwSize si,li;

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
    
    /* INPUT 0 msldem_imgpath */
    msldemcedge_inImage = set_mxInt8Matrix(prhs[0]);
    
    /* INPUT 1 msldem_imgpath */
    msldemc_imFOVmask = set_mxInt8Matrix(prhs[1]);
    
    L_demc = mxGetM(prhs[1]);
    S_demc = mxGetN(prhs[1]);
    
    /* OUTPUT 0/1/2 north-east-elevation */
    plhs[0] = mxCreateNumericMatrix(L_demc,S_demc,mxINT8_CLASS,mxREAL);
    msldemc_inImage = set_mxInt8Matrix(plhs[0]);
    
    // Initialize matrices
    for(si=0;si<S_demc;si++){
        for(li=0;li<L_demc;li++){
            msldemc_inImage[si][li] = 0;
        }
    }
    // printf("%d\n",S_im);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    find_hidden((int32_T) S_demc, (int32_T) L_demc, msldemcedge_inImage,
            msldemc_imFOVmask, msldemc_inImage);
    
    /* free memories */
    mxFree(msldemc_imFOVmask);
    mxFree(msldemcedge_inImage);
    mxFree(msldemc_inImage);
    
}
