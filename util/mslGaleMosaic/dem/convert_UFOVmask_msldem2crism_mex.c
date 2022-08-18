/* =====================================================================
 * convert_UFOVmask_msldem2crism_mex.c
 * Convert the UFOV mask image in the MSDEM resolution to crism resolution.
 * 
 * INPUTS:
 * 0 msldemc_imUFOVmask [Ldem x Sdem] int8 matrix (basically boolean)
 * 1 msldemc_northing  [Ldem] double northing array of msldemc
 * 2 msldemc_easting   [Sdem] double easting array of msldemc
 * 3 crism_n0    Scalar northing at the index [0,0] of crism image
 * 4 crism_e0    Scalar easting at the index [0,0] of crism image
 * 5 crism_dn    Scalar northing step size of crism image
 * 6 crism_de    Scalar easting step size of crism image
 * 7 crism_nsz   Scalar, integer size of crism (northing, y)
 * 8 crism_esz   Scalar, integer size of crism (easting, x)
 * 
 * (The index [0,0] here corresponds to the index [1,1] in MATLAB.)
 * 
 * OUTPUTS:
 * 0 crism_imUFOVmask [Lcrism x Scrism] int16 matrix
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
void convert_imUFOVmask_msl2crism(int8_T **msldemc_imUFOVmask, int16_T **crism_imUFOVmask,
        int32_T Ldem, int32_T Sdem, double *msldemc_northing, double *msldemc_easting,
        int32_T Lcrism, int32_T Scrism, double crism_n0, double crism_e0, double crism_dn, double crism_de)
{
    size_t si,li;
    int32_T *dem_north_cridx;
    int32_T *dem_east_cridx;
    double tmp;
    
    dem_north_cridx = (int32_T*) malloc(sizeof(int32_T) * (size_t) Ldem);
    for(li=0;li<Ldem;li++){
        tmp = round((crism_n0-msldemc_northing[li])/crism_dn);
        if(tmp>-1 || tmp<Lcrism){
            dem_north_cridx[li] = (int32_T) tmp;
        }else{
            dem_north_cridx[li] = -2;
        }
    }
    
    dem_east_cridx = (int32_T*) malloc(sizeof(int32_T) * (size_t) Sdem); 
    for(si=0;si<Sdem;si++){
        tmp = round((msldemc_easting[si] - crism_e0)/crism_de);
        if(tmp > -1 && tmp < Scrism){
            dem_east_cridx[si] = (int32_T) tmp;
        } else {
            dem_east_cridx[si] = -2;
        }
    }
    
    for(si=0;si<Sdem;si++){
        if(dem_east_cridx[si]>-1){
            for(li=0;li<Ldem;li++){
                if(dem_north_cridx[li]>-1){
                    if(msldemc_imUFOVmask[si][li]>0){
                        crism_imUFOVmask[dem_east_cridx[si]][dem_north_cridx[li]]++;
                    }
                }
            }
        }
    }
    free(dem_north_cridx);
    free(dem_east_cridx);
    
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int8_T **msldemc_imUFOVmask;
    int16_T **crism_imUFOVmask;
    mwSize Ldem,Sdem;
    double *msldemc_northing,*msldemc_easting;
    double crism_n0,crism_e0,crism_dn,crism_de;
        mwSize Lcrism,Scrism;
    
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
    
    /* INPUT 0 msldemc_imUFOVmask */
    msldemc_imUFOVmask = set_mxInt8Matrix(prhs[0]);
    Ldem = (mwSize) mxGetM(prhs[0]);
    Sdem = (mwSize) mxGetN(prhs[0]);
    
    /* INPUT 1-4 msldemc projection parameters */
    msldemc_northing = mxGetDoubles(prhs[1]);
    msldemc_easting  = mxGetDoubles(prhs[2]);
    
    /* INPUT 5-10 crism projection parameters */
    crism_n0 = mxGetScalar(prhs[3]);
    crism_e0 = mxGetScalar(prhs[4]);
    crism_dn = mxGetScalar(prhs[5]);
    crism_de = mxGetScalar(prhs[6]);
    Lcrism   = (mwSize) mxGetScalar(prhs[7]);
    Scrism   = (mwSize) mxGetScalar(prhs[8]);

    /* OUTPUT 0 crism_imUFOVmask */
    plhs[0] = mxCreateNumericMatrix(Lcrism,Scrism,mxINT16_CLASS,mxREAL);
    crism_imUFOVmask = set_mxInt16Matrix(plhs[0]);
    
    // Initialize matrices
    for(si=0;si<Scrism;si++){
        for(li=0;li<Lcrism;li++){
            crism_imUFOVmask[si][li] = 0;
        }
    }
    
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    convert_imUFOVmask_msl2crism(msldemc_imUFOVmask,crism_imUFOVmask,
            (int32_T) Ldem, (int32_T) Sdem,
            msldemc_northing, msldemc_easting,
            (int32_T) Lcrism, (int32_T) Scrism,
            crism_n0, crism_e0, crism_dn, crism_de);
    
    /* free memories */
    mxFree(msldemc_imUFOVmask);
    mxFree(crism_imUFOVmask);
    
}
