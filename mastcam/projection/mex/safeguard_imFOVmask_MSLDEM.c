/* =====================================================================
 * get_imxy_MSLDEM_mex.c
 * Read MSLDEM image data
 * Perform projection of mastcam pixels onto MSLDEM data
 * 
 * INPUTS:
 * 0 msldem_imgpath        char*
 * 1 msldem_header         struct
 * 2 msldemc_header        struct
 * 3 msldemc_northing      Double array [L_demc]
 * 4 msldemc_easting       Double array [S_demc]
 * 5 msldem_imFOVmaskd    int8 [L_demc x S_demc]
 * 8 cammdl               CAHVOR model class
 * 
 * 
 * OUTPUTS:
 * 0 msldem_safeguard      int8 [L_demc * S_demc]
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

/* main computation routine */
void safeguard_imFOVmask(char *msldem_imgpath, EnviHeader msldem_hdr, 
        int32_T msldemc_imxy_sample_offset, int32_T msldemc_imxy_line_offset,
        int32_T msldemc_samples, int32_T msldemc_lines,
        double *msldemc_northing, double *msldemc_easting,
        int8_T **msldemc_imFOVmaskd, 
        int32_T S_im, int32_T L_im,
        double *cam_C, double *cam_A,
        int8_T **msldemc_safeguard_mask)
{
    int32_T c,l;
    int32_T msldemc_samplesm1,msldemc_linesm1;
    int16_T sxy;
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
                    if(c==0){
                        if(msldemc_imFOVmaskd[c+1][l] || msldemc_imFOVmaskd[c][l+1] || msldemc_imFOVmaskd[c+1][l+1])
                            msldemc_safeguard_mask[c][l] = 1;
                    } else if(c==msldemc_samplesm1){
                        if(msldemc_imFOVmaskd[c-1][l] || msldemc_imFOVmaskd[c-1][l+1] || msldemc_imFOVmaskd[c][l+1])
                            msldemc_safeguard_mask[c][l] = 1;
                    } else{
                        if(msldemc_imFOVmaskd[c-1][l] || msldemc_imFOVmaskd[c+1][l] 
                                || msldemc_imFOVmaskd[c-1][l+1] || msldemc_imFOVmaskd[c][l+1] || msldemc_imFOVmaskd[c+1][l+1])
                            msldemc_safeguard_mask[c][l] = 1;
                    }
                } else if(l==msldemc_linesm1){
                    if(c==0){
                        if(msldemc_imFOVmaskd[c][l-1] || msldemc_imFOVmaskd[c+1][l-1] || msldemc_imFOVmaskd[c+1][l])
                            msldemc_safeguard_mask[c][l] = 1;
                    } else if(c==msldemc_samplesm1){
                        if(msldemc_imFOVmaskd[c-1][l-1] || msldemc_imFOVmaskd[c][l-1] || msldemc_imFOVmaskd[c-1][l])
                            msldemc_safeguard_mask[c][l] = 1;
                    } else{
                        if(msldemc_imFOVmaskd[c-1][l-1] || msldemc_imFOVmaskd[c][l-1] || msldemc_imFOVmaskd[c+1][l-1] 
                                || msldemc_imFOVmaskd[c-1][l] || msldemc_imFOVmaskd[c+1][l])
                            msldemc_safeguard_mask[c][l] = 1;
                    }
                } else {
                    if(c==0){
                        if(msldemc_imFOVmaskd[c][l-1] || msldemc_imFOVmaskd[c+1][l-1] || msldemc_imFOVmaskd[c+1][l] || msldemc_imFOVmaskd[c][l+1] || msldemc_imFOVmaskd[c+1][l+1])
                            msldemc_safeguard_mask[c][l] = 1;
                    } else if(c==msldemc_samplesm1){
                        if(msldemc_imFOVmaskd[c-1][l-1] || msldemc_imFOVmaskd[c][l-1] || msldemc_imFOVmaskd[c-1][l] || msldemc_imFOVmaskd[c-1][l+1] || msldemc_imFOVmaskd[c][l+1])
                            msldemc_safeguard_mask[c][l] = 1;
                    } else{
                        if(msldemc_imFOVmaskd[c-1][l-1] || msldemc_imFOVmaskd[c][l-1] || msldemc_imFOVmaskd[c+1][l-1] 
                                || msldemc_imFOVmaskd[c-1][l] || msldemc_imFOVmaskd[c+1][l]
                                || msldemc_imFOVmaskd[c-1][l+1] || msldemc_imFOVmaskd[c][l+1] || msldemc_imFOVmaskd[c+1][l+1])
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
    mxArray *cam_C_mxar, *cam_C_mxard;
    mxArray *cam_A_mxar, *cam_A_mxard;
    double *cam_C, *cam_A;
    
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
    msldem_hdr.samples = (int32_T) mxGetScalar(mxGetField(prhs[1],0,"samples"));
    msldem_hdr.lines = (int32_T) mxGetScalar(mxGetField(prhs[1],0,"lines"));
    msldem_hdr.bands = (int32_T) mxGetScalar(mxGetField(prhs[1],0,"bands"));
    msldem_hdr.data_type = (int32_T) mxGetScalar(mxGetField(prhs[1],0,"data_type"));
    msldem_hdr.byte_order = (int32_T) mxGetScalar(mxGetField(prhs[1],0,"byte_order"));
    msldem_hdr.data_ignore_value = mxGetScalar(mxGetField(prhs[1],0,"data_ignore_value"));
    
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
    cam_C_mxar = mxGetProperty(prhs[8],0,"C");
    cam_C_mxard = mxDuplicateArray(cam_C_mxar);
    cam_C = mxGetDoubles(cam_C_mxard);
    
    cam_A_mxar = mxGetProperty(prhs[8],0,"A");
    cam_A_mxard = mxDuplicateArray(cam_A_mxar);
    cam_A = mxGetDoubles(cam_A_mxard);
    
    
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
        cam_C, cam_A,
        msldemc_safeguard_mask);
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldemc_imFOVmaskd);
    mxFree(msldemc_safeguard_mask);
    
}
