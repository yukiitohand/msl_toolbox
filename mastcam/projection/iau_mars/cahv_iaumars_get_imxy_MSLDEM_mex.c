/* =====================================================================
 * cahv_iaumars_get_imxy_MSLDEM_mex.c
 * Read MSLDEM image data
 * Perform projection of mastcam pixels onto MSLDEM data
 * 
 * INPUTS:
 * 0 msldem_imgpath        char*
 * 1 msldem_header         struct
 * 2 mslrad_offset         Double Scalar
 * 3 msldemc_header        struct
 * 4 msldemc_latitude      Double array [L_demc]
 * 5 msldemc_longitude     Double array [S_demc]
 * 6 msldemc_imFOVmaskd    int8 [L_demc x S_demc]
 * 7 cammdl               CAHVOR model class
 * 
 * 
 * OUTPUTS:
 * 0 msldemc_imx           Doubles [L_demc x S_demc]
 * 1 msldemc_imy           Doubles [L_demc x S_demc]
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
#include "cahvor.h"
#include "mex_create_array.h"

/* main computation routine */
void get_imxy_MSLDEM(char *msldem_imgpath, EnviHeader msldem_hdr, double mslrad_offset,
        int32_T msldemc_imxy_sample_offset, int32_T msldemc_imxy_line_offset,
        int32_T msldemc_samples, int32_T msldemc_lines,
        double *msldemc_latitude, double *msldemc_longitude,
        int8_T **msldemc_imFOVmaskd, 
        CAHV_MODEL cahv_mdl,
        double **msldemc_imx, double **msldemc_imy)
{
    int32_T c,l;
    int16_T sxy;
    long skip_pri;
    long skip_l, skip_r;
    float *elevl;
    size_t ncpy;
    size_t sz=sizeof(float);
    FILE *fid;
    double radius_cl;
    double pmcx,pmcy,pmcz;
    double apmc,hpmc,vpmc;
    double *cam_C, *cam_A, *cam_H, *cam_V;
    
    double *cos_lon, *sin_lon;
    double cos_latl, sin_latl;
    double x_iaumars,y_iaumars,z_iaumars;
    
    cam_C = cahv_mdl.C; cam_A = cahv_mdl.A; cam_H = cahv_mdl.H; cam_V = cahv_mdl.V;
    
    cos_lon = (double*) malloc(sizeof(double) * (size_t) msldemc_samples);
    sin_lon = (double*) malloc(sizeof(double) * (size_t) msldemc_samples);
    for(c=0;c<msldemc_samples;c++){
        cos_lon[c] = cos(msldemc_longitude[c]);
        sin_lon[c] = sin(msldemc_longitude[c]);
    }
    
    fid = fopen(msldem_imgpath,"rb");
    
    /* skip lines */
    skip_pri = (long) msldem_hdr.samples * (long) msldemc_imxy_line_offset * (long) sz;
    // printf("%d*%d*%d=%ld\n",msldem_header.samples,msldemc_imxy_line_offset,s,skip_pri);
    fseek(fid,skip_pri,SEEK_CUR);
    
    /* read the data */
    ncpy = sz * (size_t) msldemc_samples;
    elevl = (float*) malloc(ncpy);
    skip_l = (long) sz * (long) msldemc_imxy_sample_offset;
    skip_r = ((long) msldem_hdr.samples - (long) msldemc_samples)* (long) sz - skip_l;
    
    
    // printf("%d,%d,%d\n",skip_l,msldemc_samples*s,skip_r);
    for(l=0;l<msldemc_lines;l++){
        fseek(fid,skip_l,SEEK_CUR);
        fread(elevl,sz,msldemc_samples,fid);
        fseek(fid,skip_r,SEEK_CUR);
        
        cos_latl = cos(msldemc_latitude[l]);
        sin_latl = sin(msldemc_latitude[l]);
        for(c=0;c<msldemc_samples;c++){
            if(msldemc_imFOVmaskd[c][l]>1){
                radius_cl = (double) elevl[c] + mslrad_offset;
                x_iaumars = radius_cl * cos_latl * cos_lon[c];
                y_iaumars = radius_cl * cos_latl * sin_lon[c];
                z_iaumars = radius_cl * sin_latl;
                
                pmcx = x_iaumars - cam_C[0];
                pmcy = y_iaumars - cam_C[1];
                pmcz = z_iaumars - cam_C[2];
                
                apmc = cam_A[0] * pmcx + cam_A[1] * pmcy + cam_A[2] * pmcz;
                hpmc = cam_H[0] * pmcx + cam_H[1] * pmcy + cam_H[2] * pmcz;
                vpmc = cam_V[0] * pmcx + cam_V[1] * pmcy + cam_V[2] * pmcz;
                
                msldemc_imx[c][l] = hpmc / apmc;
                msldemc_imy[c][l] = vpmc / apmc;
            }

        }

    }
    free(elevl);
    free(cos_lon);
    free(sin_lon);
    fclose(fid);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *msldem_imgpath;
    EnviHeader msldem_hdr;
    double mslrad_offset;
    CAHV_MODEL cahv_mdl;
    mwSize msldemc_imxy_sample_offset,msldemc_imxy_line_offset;
    double *msldemc_latitude;
    double *msldemc_longitude;
    int8_T **msldemc_imFOVmaskd;    
    double **msldemc_imx;
    double **msldemc_imy;
    

    
    mwSize si,li;
    mwSize msldemc_samples, msldemc_lines;

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
    msldem_imgpath = mxArrayToString(prhs[0]);
    
    /* INPUT 1 msldem_header */
    msldem_hdr = mxGetEnviHeader(prhs[1]);
    
    mslrad_offset = mxGetScalar(prhs[2]);
    
    /* INPUT 2 msldemc_sheader*/
    msldemc_imxy_sample_offset = (mwSize) mxGetScalar(mxGetField(prhs[3],0,"sample_offset"));
    msldemc_imxy_line_offset = (mwSize) mxGetScalar(mxGetField(prhs[3],0,"line_offset"));
    msldemc_samples = (mwSize) mxGetScalar(mxGetField(prhs[3],0,"samples"));
    msldemc_lines = (mwSize) mxGetScalar(mxGetField(prhs[3],0,"lines"));
    
    /* INPUT 3/4 msldem northing easting */
    msldemc_latitude  = mxGetDoubles(prhs[4]);
    msldemc_longitude = mxGetDoubles(prhs[5]);
    

    /* INPUT 6 msldem imFOV */
    msldemc_imFOVmaskd = set_mxInt8Matrix(prhs[6]);
    
    /* INPUT 7 camera model */
    cahv_mdl = mxGet_CAHV_MODEL(prhs[7]);
    
    /* OUTPUT 0/1 msldem imxy */
    plhs[0] = mxCreateDoubleMatrix(msldemc_lines,msldemc_samples,mxREAL);
    msldemc_imx = set_mxDoubleMatrix(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(msldemc_lines,msldemc_samples,mxREAL);
    msldemc_imy = set_mxDoubleMatrix(plhs[1]);  
    

    
    // Initialize matrices
    for(si=0;si<msldemc_samples;si++){
        for(li=0;li<msldemc_lines;li++){
            msldemc_imx[si][li] = NAN;
            msldemc_imy[si][li] = NAN;            
        }
    }
    
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    get_imxy_MSLDEM(msldem_imgpath, msldem_hdr, mslrad_offset,
        (int32_T) msldemc_imxy_sample_offset, (int32_T) msldemc_imxy_line_offset,
        (int32_T) msldemc_samples, (int32_T) msldemc_lines,
        msldemc_latitude, msldemc_longitude,
        msldemc_imFOVmaskd, 
        cahv_mdl,
        msldemc_imx, msldemc_imy);
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldemc_imx);
    mxFree(msldemc_imy);
    mxFree(msldemc_imFOVmaskd);
    
}
