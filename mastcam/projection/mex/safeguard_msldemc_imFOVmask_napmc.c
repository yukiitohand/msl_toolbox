/* =====================================================================
 * safeguard_msldemc_imFOVmask_napmc.c
 * Detect potential pixels within FOV the pmc of which havs negative dot 
 * products with the camera axis vector.
 * 
 * INPUTS:
 * 0 msldem_imgpath        char*
 * 1 msldem_header         struct
 * 2 msldemc_header        struct
 * 3 msldemc_northing      Double array [L_demc]
 * 4 msldemc_easting       Double array [S_demc]
 * 5 msldem_imFOVmaskd     int8 [L_demc x S_demc]
 * 8 cmmdl                 CAHV model class
 * 
 * 
 * OUTPUTS:
 * 0 msldemc_safeguard     int8 [L_demc * S_demc]
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

/* main computation routine */
void safeguard_imFOVmask(char *msldem_imgpath, EnviHeader msldem_hdr, 
        int32_T msldemc_imxy_sample_offset, int32_T msldemc_imxy_line_offset,
        int32_T msldemc_samples, int32_T msldemc_lines,
        double *msldemc_northing, double *msldemc_easting,
        int8_T **msldemc_imFOVmaskd, 
        int32_T S_im, int32_T L_im,
        CAHV_MODEL cahv_mdl,
        int8_T **msldemc_safeguard_mask)
{
    int32_T c,l;
    int32_T msldemc_samplesm1,msldemc_linesm1;
    long skip_pri;
    long skip_l, skip_r;
    float *elevl;
    size_t ncpy;
    size_t sz=sizeof(float);
    FILE *fid;
    double dem_cl;
    double pmcx,pmcy,pmcz;
    double apmcx,apmcy,apmcz;
    double *APmCys;
    double apmc;
    double *cam_C, double *cam_A;
    cam_C = cahv_mdl.C; cam_A = cahv_mdl.A;
    
    msldemc_samplesm1 = msldemc_samples-1;
    msldemc_linesm1 = msldemc_lines - 1;
    
    APmCys = (double*) malloc((size_t) msldemc_samples * sizeof(double));
    for(c=0;c<msldemc_samples;c++){
        pmcy  = msldemc_easting[c] - cam_C[1];
        APmCys[c] = cam_A[1] * pmcy;
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
        pmcx  = msldemc_northing[l] - cam_C[0];
        apmcx = cam_A[0] * pmcx;
        for(c=0;c<msldemc_samples;c++){
            dem_cl = (double) elevl[c];
            apmcz = cam_A[2] * pmcz;
            apmcy = APmCys[c];
            apmc = apmcx + apmcy + apmcz;
            if(apmc<0){
                if(l==0){
                    printf("a\n");
                    if(c==0){
                        if(msldemc_imFOVmaskd[c+1][l]>2 || msldemc_imFOVmaskd[c][l+1]>2 || msldemc_imFOVmaskd[c+1][l+1]>2)
                            msldemc_safeguard_mask[c][l] = 1;
                    } else if(c==msldemc_samplesm1){
                        if(msldemc_imFOVmaskd[c-1][l]>2 || msldemc_imFOVmaskd[c-1][l+1]>2 || msldemc_imFOVmaskd[c][l+1]>2)
                            msldemc_safeguard_mask[c][l] = 1;
                    } else{
                        if(msldemc_imFOVmaskd[c-1][l]>2 || msldemc_imFOVmaskd[c+1][l]>2
                                || msldemc_imFOVmaskd[c-1][l+1]>2 || msldemc_imFOVmaskd[c][l+1]>2 || msldemc_imFOVmaskd[c+1][l+1]>2)
                            msldemc_safeguard_mask[c][l] = 1;
                    }
                } else if(l==msldemc_linesm1){
                    if(c==0){
                        if(msldemc_imFOVmaskd[c][l-1]>2 || msldemc_imFOVmaskd[c+1][l-1]>2 || msldemc_imFOVmaskd[c+1][l]>2)
                            msldemc_safeguard_mask[c][l] = 1;
                    } else if(c==msldemc_samplesm1){
                        if(msldemc_imFOVmaskd[c-1][l-1]>2 || msldemc_imFOVmaskd[c][l-1]>2 || msldemc_imFOVmaskd[c-1][l]>2)
                            msldemc_safeguard_mask[c][l] = 1;
                    } else{
                        if(msldemc_imFOVmaskd[c-1][l-1]>2 || msldemc_imFOVmaskd[c][l-1]>2 || msldemc_imFOVmaskd[c+1][l-1]>2
                                || msldemc_imFOVmaskd[c-1][l]>2 || msldemc_imFOVmaskd[c+1][l]>2)
                            msldemc_safeguard_mask[c][l] = 1;
                    }
                } else {
                    if(c==0){
                        if(msldemc_imFOVmaskd[c][l-1]>2 || msldemc_imFOVmaskd[c+1][l-1]>2 || msldemc_imFOVmaskd[c+1][l]>2 || msldemc_imFOVmaskd[c][l+1]>2 || msldemc_imFOVmaskd[c+1][l+1]>2)
                            msldemc_safeguard_mask[c][l] = 1;
                    } else if(c==msldemc_samplesm1){
                        if(msldemc_imFOVmaskd[c-1][l-1]>2 || msldemc_imFOVmaskd[c][l-1]>2 || msldemc_imFOVmaskd[c-1][l]>2 || msldemc_imFOVmaskd[c-1][l+1]>2 || msldemc_imFOVmaskd[c][l+1]>2)
                            msldemc_safeguard_mask[c][l] = 1;
                    } else{
                        if(msldemc_imFOVmaskd[c-1][l-1]>2 || msldemc_imFOVmaskd[c][l-1]>2 || msldemc_imFOVmaskd[c+1][l-1]>2 
                                || msldemc_imFOVmaskd[c-1][l]>2 || msldemc_imFOVmaskd[c+1][l]>2
                                || msldemc_imFOVmaskd[c-1][l+1]>2 || msldemc_imFOVmaskd[c][l+1]>2 || msldemc_imFOVmaskd[c+1][l+1]>2)
                            msldemc_safeguard_mask[c][l] = 1;
                    }
                }
            }
        }

    }
    free(elevl);
    free(APmCys);
    fclose(fid);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *msldem_imgpath;
    EnviHeader msldem_hdr;
    mwSize msldemc_imxy_sample_offset,msldemc_imxy_line_offset;
    double *msldemc_northing;
    double *msldemc_easting;
    int8_T **msldemc_imFOVmaskd;
    
    mwSize S_im,L_im;
    CAHV_MODEL cahv_mdl;
    
    int8_T **msldemc_safeguard_mask;
    

    
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
    
    /* INPUT 2 msldemc_sheader*/
    msldemc_imxy_sample_offset = (mwSize) mxGetScalar(mxGetField(prhs[2],0,"sample_offset"));
    msldemc_imxy_line_offset = (mwSize) mxGetScalar(mxGetField(prhs[2],0,"line_offset"));
    msldemc_samples = (mwSize) mxGetScalar(mxGetField(prhs[2],0,"samples"));
    msldemc_lines = (mwSize) mxGetScalar(mxGetField(prhs[2],0,"lines"));
    
    /* INPUT 3/4 msldem northing easting */
    msldemc_northing = mxGetDoubles(prhs[3]);
    msldemc_easting = mxGetDoubles(prhs[4]);
    

    /* INPUT 5 msldem imFOV */
    msldemc_imFOVmaskd = set_mxInt8Matrix(prhs[5]);
    
    
    /* INPUT 7 camera model */
    cahv_mdl = mxGet_CAHV_MODEL(prhs[6]);
    
    
    /* OUTPUT 0/1 msldem imxy */
    plhs[0] = mxCreateNumericMatrix(msldemc_lines,msldemc_samples,mxINT8_CLASS,mxREAL);
    msldemc_safeguard_mask = set_mxInt8Matrix(plhs[0]);
    
    // Initialize matrices
    for(si=0;si<msldemc_samples;si++){
        for(li=0;li<msldemc_lines;li++){
            msldemc_safeguard_mask[si][li] = 0;           
        }
    }
    
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    safeguard_imFOVmask(msldem_imgpath, msldem_hdr, 
        (int32_T) msldemc_imxy_sample_offset, (int32_T) msldemc_imxy_line_offset,
        (int32_T) msldemc_samples, (int32_T) msldemc_lines,
        msldemc_northing, msldemc_easting,
        msldemc_imFOVmaskd, 
        (int32_T) S_im, (int32_T) L_im,
        cahv_mdl,
        msldemc_safeguard_mask);
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldemc_imFOVmaskd);
    mxFree(msldemc_safeguard_mask);
    
}
