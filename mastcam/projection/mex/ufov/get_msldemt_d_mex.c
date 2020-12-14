/* =====================================================================
 * get_msldemt_d_mex.c
 * GET the image coordinate of vertices of each pixel 
 * in the DEM image.
 * 
 * INPUTS:
 * 0 msldem_imgpath        char*
 * 1 msldem_header         struct
 * 2 msldemc_header        struct
 * 3 msldemc_northing      Double array [L_demc]
 * 4 msldemc_easting       Double array [S_demc]
 * 5 msldemc_imFOVmaskd    int8 [L_demc x S_demc]
 * 6 S_im                  int
 * 7 L_im                  int
 * 8 cammdl               CAHVOR model class
 * 
 * 
 * OUTPUTS:
 * 0 msldemc_imxm           Doubles [(L_demc-1) x (S_demc-1)]
 * 1 msldemc_imym           Doubles [(L_demc-1) x (S_demc-1)]
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
void get_imxyclm_MSLDEM(char *msldem_imgpath, EnviHeader msldem_hdr, 
        int32_T msldemc_imxy_sample_offset, int32_T msldemc_imxy_line_offset,
        int32_T msldemc_samples, int32_T msldemc_lines,
        double *msldemc_northing, double *msldemc_easting,
        int8_T **msldemc_imFOVmaskd, int8_T **msldemc_imUFOVmask,
        CAHV_MODEL cahv_mdl,
        double* msldemt_northing, double* msldemt_easting, double** msldemt_img, int8_T **msldemt_imFOVmask)
{
    int32_T c,l;
    int16_T sxy;
    long skip_pri;
    long skip_l, skip_r;
    float *elevl, *elevlp1;
    size_t ncpy;
    size_t sz=sizeof(float);
    FILE *fid;
    double dem_cl;
    double pmcx,pmcy,pmcz;
    double apmcx,apmcy,apmcz;
    double *APmCys;
    double apmc;
    double *cam_C, *cam_A;
    int32_T msldemc_samplesm1, msldemc_linesm1;
    
    msldemc_samplesm1 = msldemc_samples - 1;
    msldemc_linesm1 = msldemc_lines - 1;
    
    cam_C = cahv_mdl.C; cam_A = cahv_mdl.A;
    
    APmCys = (double*) malloc((size_t) msldemc_samplesm1 * sizeof(double));
    for(c=0;c<msldemc_samplesm1;c++){
        msldemt_easting[c] = 0.5*(msldemc_easting[c]+msldemc_easting[c+1]);
        pmcy  = msldemt_easting[c] - cam_C[1];
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
    elevlp1 = (float*) malloc(sz*msldemc_samples);
    skip_l = (long) sz * (long) msldemc_imxy_sample_offset;
    skip_r = ((long) msldem_hdr.samples - (long) msldemc_samples)* (long) sz - skip_l;
    
    fseek(fid,skip_l,SEEK_CUR);
    fread(elevlp1,sz,msldemc_samples,fid);
    fseek(fid,skip_r,SEEK_CUR);
    
    // printf("%d,%d,%d\n",skip_l,msldemc_samples*s,skip_r);
    for(l=0;l<msldemc_linesm1;l++){
        memcpy(elevl,elevlp1,ncpy);
        fseek(fid,skip_l,SEEK_CUR);
        fread(elevlp1,sz,msldemc_samples,fid);
        fseek(fid,skip_r,SEEK_CUR);
        
        msldemt_northing[l] = 0.5*(msldemc_northing[l] + msldemc_northing[l+1]);
        pmcx  = msldemt_northing[l] - cam_C[0];
        apmcx = cam_A[0]*pmcx;
        for(c=0;c<msldemc_samplesm1;c++){
            if((msldemc_imFOVmaskd[c][l]>0) || (msldemc_imFOVmaskd[c+1][l]>0) || (msldemc_imFOVmaskd[c][l+1]>0) || (msldemc_imFOVmaskd[c+1][l+1]>0)){
                if( (msldemc_imUFOVmask[c][l]==0) || (msldemc_imUFOVmask[c+1][l]==0) || (msldemc_imUFOVmask[c][l+1]==0) || (msldemc_imUFOVmask[c+1][l+1]==0) ){
                    
                    dem_cl = 0.5*( (double) elevl[c+1] + (double) elevlp1[c] );
                    pmcz  = -dem_cl-cam_C[2];
                    apmcz = cam_A[2] * pmcz;
                    apmcy = APmCys[c];
                    apmc = apmcx + apmcy + apmcz;

                    msldemt_img[c][l] = dem_cl;

                    if(apmc>0){
                        msldemt_imFOVmask[c][l] = 2;
                    } else {
                        msldemt_imFOVmask[c][l] = 1;
                    }
                }
                    
            }

        }

    }
    free(elevl);
    free(elevlp1);
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
    int8_T **msldemc_imUFOVmask;
    CAHV_MODEL cahv_mdl;
    
    double **msldemt_img;
    double *msldemt_northing,*msldemt_easting;
    int8_T **msldemt_imFOVmask;
    

    
    mwSize si,li;
    mwSize msldemc_samples, msldemc_lines;
    mwSize msldemc_samplesm1,msldemc_linesm1;

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
    
    /* INPUT 5 msldem imFOV */
    msldemc_imUFOVmask = set_mxInt8Matrix(prhs[6]);
    
    
    /* INPUT 7 camera model */
    cahv_mdl = mxGet_CAHV_MODEL(prhs[7]);
    
    /* OUTPUT 0 msldemt */
    msldemc_samplesm1 = msldemc_samples-1;
    msldemc_linesm1 = msldemc_lines-1;
    plhs[0] = mxCreateDoubleMatrix(msldemc_linesm1,msldemc_samplesm1,mxREAL);
    msldemt_img = set_mxDoubleMatrix(plhs[0]);
    
    plhs[1] = mxCreateNumericArray(1,&msldemc_linesm1,mxDOUBLE_CLASS,mxREAL);
    msldemt_northing = mxGetDoubles(plhs[1]);
    
    plhs[2] = mxCreateNumericArray(1,&msldemc_samplesm1,mxDOUBLE_CLASS,mxREAL);
    msldemt_easting = mxGetDoubles(plhs[2]);
    
    plhs[3] = mxCreateNumericMatrix(msldemc_linesm1,msldemc_samplesm1,mxINT8_CLASS,mxREAL);
    msldemt_imFOVmask = set_mxInt8Matrix(plhs[3]);
    
    
    // Initialize matrices
    for(li=0;li<msldemc_linesm1;li++){
        msldemt_northing[li] = NAN;         
    }
    
    for(si=0;si<msldemc_samplesm1;si++){
        msldemt_easting[si] = NAN;
        for(li=0;li<msldemc_linesm1;li++){
            msldemt_img[si][li] = NAN;
            msldemt_imFOVmask[si][li] = 0;
        }
    }
    
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    get_imxyclm_MSLDEM(msldem_imgpath, msldem_hdr, 
        (int32_T) msldemc_imxy_sample_offset, (int32_T) msldemc_imxy_line_offset,
        (int32_T) msldemc_samples, (int32_T) msldemc_lines,
        msldemc_northing, msldemc_easting,
        msldemc_imFOVmaskd,msldemc_imUFOVmask,
        cahv_mdl,
        msldemt_northing,msldemt_easting, msldemt_img, msldemt_imFOVmask);
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldemc_imFOVmaskd);
    mxFree(msldemc_imUFOVmask);
    mxFree(msldemt_img);
    mxFree(msldemt_imFOVmask);
    
    
}
