/* =====================================================================
 * cahv_get_imxy_MSLDEM_mex.c
 * Read MSLDEM image data
 * Perform projection of mastcam pixels onto MSLDEM data
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
void get_imxy_MSLDEM(char *msldem_imgpath, EnviHeader msldem_hdr, 
        int32_T msldemc_imxy_sample_offset, int32_T msldemc_imxy_line_offset,
        int32_T msldemc_samples, int32_T msldemc_lines,
        double *msldemc_northing, double *msldemc_easting,
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
    double dem_cl;
    double pmcx,pmcy,pmcz;
    double apmcx,apmcy,apmcz;
    double hpmcx,hpmcy,hpmcz;
    double vpmcx,vpmcy,vpmcz;
    double *APmCys,*HPmCys,*VPmCys;
    double apmc,hpmc,vpmc;
    double *cam_C, *cam_A, *cam_H, *cam_V;
    
    cam_C = cahv_mdl.C; cam_A = cahv_mdl.A; cam_H = cahv_mdl.H; cam_V = cahv_mdl.V;
    
    
    
    APmCys = (double*) malloc((size_t) msldemc_samples * sizeof(double));
    HPmCys = (double*) malloc((size_t) msldemc_samples * sizeof(double));
    VPmCys = (double*) malloc((size_t) msldemc_samples * sizeof(double));
    for(c=0;c<msldemc_samples;c++){
        pmcy  = msldemc_easting[c] - cam_C[1];
        APmCys[c] = cam_A[1] * pmcy;
        HPmCys[c] = cam_H[1] * pmcy;
        VPmCys[c] = cam_V[1] * pmcy;
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
        hpmcx = cam_H[0] * pmcx;
        vpmcx = cam_V[0] * pmcx;
        for(c=0;c<msldemc_samples;c++){
            if(msldemc_imFOVmaskd[c][l]>1){
                dem_cl = (double) elevl[c];
                pmcz  = -dem_cl-cam_C[2];
                apmcz = cam_A[2] * pmcz;
                hpmcz = cam_H[2] * pmcz;
                vpmcz = cam_V[2] * pmcz;

                apmcy = APmCys[c];
                hpmcy = HPmCys[c];
                vpmcy = VPmCys[c];

                apmc = apmcx + apmcy + apmcz;
                hpmc = hpmcx + hpmcy + hpmcz;
                vpmc = vpmcx + vpmcy + vpmcz;
                
                msldemc_imx[c][l] = hpmc / apmc;
                msldemc_imy[c][l] = vpmc / apmc;
            }

        }

    }
    free(elevl);
    free(APmCys);
    free(HPmCys);
    free(VPmCys);
    fclose(fid);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *msldem_imgpath;
    EnviHeader msldem_hdr;
    CAHV_MODEL cahv_mdl;
    mwSize msldemc_imxy_sample_offset,msldemc_imxy_line_offset;
    double *msldemc_northing;
    double *msldemc_easting;
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
    //cahv_mdl = mxGet_CAHV_MODEL(prhs[6]);
    // printf("%f,%f,%f\n",cahv_mdl.C[0],cahv_mdl.C[1],cahv_mdl.C[2]);
    //printf("%f\n",cahv_mdl.Hdash[0]);
    
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
    get_imxy_MSLDEM(msldem_imgpath, msldem_hdr, 
        (int32_T) msldemc_imxy_sample_offset, (int32_T) msldemc_imxy_line_offset,
        (int32_T) msldemc_samples, (int32_T) msldemc_lines,
        msldemc_northing, msldemc_easting,
        msldemc_imFOVmaskd, 
        cahv_mdl,
        msldemc_imx, msldemc_imy);
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldemc_imx);
    mxFree(msldemc_imy);
    mxFree(msldemc_imFOVmaskd);
    
}
