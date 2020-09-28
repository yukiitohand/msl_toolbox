/* =====================================================================
 * get_imxycm_MSLDEM_mex.c
 * GET coordinate values in the image coordinate of column edge of each pixel 
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
 * 0 msldemc_imxm           Doubles [(L_demc) x (S_demc-1)]
 * 1 msldemc_imym           Doubles [(L_demc) x (S_demc-1)]
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
void get_imxyclm_MSLDEM(char *msldem_imgpath, EnviHeader msldem_hdr, 
        int32_T msldemc_imxy_sample_offset, int32_T msldemc_imxy_line_offset,
        int32_T msldemc_samples, int32_T msldemc_lines,
        double *msldemc_northing, double *msldemc_easting,
        int8_T **msldemc_imFOVmaskd, 
        int32_T S_im, int32_T L_im,
        double *cam_C, double *cam_A, double *cam_H, double *cam_V,
        double **msldemc_imxm, double **msldemc_imym)
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
    int32_T msldemc_samplesm1;
    
    msldemc_samplesm1 = msldemc_samples - 1;
    
    APmCys = (double*) malloc((size_t) msldemc_samplesm1 * sizeof(double));
    HPmCys = (double*) malloc((size_t) msldemc_samplesm1 * sizeof(double));
    VPmCys = (double*) malloc((size_t) msldemc_samplesm1 * sizeof(double));
    for(c=0;c<msldemc_samplesm1;c++){
        pmcy  = 0.5*(msldemc_easting[c]+msldemc_easting[c+1]) - cam_C[1];
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
        
        pmcx  = 0.5*(msldemc_northing[l] + msldemc_northing[l+1]) - cam_C[0];
        apmcx = cam_A[0] * pmcx;
        hpmcx = cam_H[0] * pmcx;
        vpmcx = cam_V[0] * pmcx;
        for(c=0;c<msldemc_samplesm1;c++){
            if((msldemc_imFOVmaskd[c][l]==2) && (msldemc_imFOVmaskd[c+1][l]==2)){
                dem_cl = 0.5*((double) (elevl[c] + elevl[c+1]));
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
                
                msldemc_imxm[c][l] = hpmc / apmc;
                msldemc_imym[c][l] = vpmc / apmc;
            }

        }

    }
    free(elevl);
    free(APmCys);
    free(HPmCys);
    free(VPmCys);
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
    mxArray *cam_H_mxar, *cam_H_mxard;
    mxArray *cam_V_mxar, *cam_V_mxard;
    double *cam_C, *cam_A, *cam_H, *cam_V;
    
    double **msldemc_imxm;
    double **msldemc_imym;
    

    
    mwSize si,li;
    mwSize msldemc_samples, msldemc_lines;
    mwSize msldemc_samplesm1;

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
    
    /* INPUT 6/7 image S_im, L_im */
    S_im = (mwSize) mxGetScalar(prhs[6]);
    L_im = (mwSize) mxGetScalar(prhs[7]);
    
    /* INPUT 7 camera model */
    cam_C_mxar = mxGetProperty(prhs[8],0,"C");
    cam_C_mxard = mxDuplicateArray(cam_C_mxar);
    cam_C = mxGetDoubles(cam_C_mxard);
    
    cam_A_mxar = mxGetProperty(prhs[8],0,"A");
    cam_A_mxard = mxDuplicateArray(cam_A_mxar);
    cam_A = mxGetDoubles(cam_A_mxard);
    
    cam_H_mxar = mxGetProperty(prhs[8],0,"H");
    cam_H_mxard = mxDuplicateArray(cam_H_mxar);
    cam_H = mxGetDoubles(cam_H_mxard);
    
    cam_V_mxar = mxGetProperty(prhs[8],0,"V");
    cam_V_mxard = mxDuplicateArray(cam_V_mxar);
    cam_V = mxGetDoubles(cam_V_mxard);
    
    
    /* OUTPUT 0/1 msldem imxy */
    msldemc_samplesm1 = msldemc_samples-1;
    plhs[0] = mxCreateDoubleMatrix(msldemc_lines,msldemc_samplesm1,mxREAL);
    msldemc_imxm = set_mxDoubleMatrix(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(msldemc_lines,msldemc_samplesm1,mxREAL);
    msldemc_imym = set_mxDoubleMatrix(plhs[1]);
    

    
    // Initialize matrices
    for(si=0;si<msldemc_samplesm1;si++){
        for(li=0;li<msldemc_lines;li++){
            msldemc_imxm[si][li] = NAN;
            msldemc_imym[si][li] = NAN;            
        }
    }
    
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    get_imxyclm_MSLDEM(msldem_imgpath, msldem_hdr, 
        (int32_T) msldemc_imxy_sample_offset, (int32_T) msldemc_imxy_line_offset,
        (int32_T) msldemc_samples, (int32_T) msldemc_lines,
        msldemc_northing, msldemc_easting,
        msldemc_imFOVmaskd, 
        (int32_T) S_im, (int32_T) L_im,
        cam_C, cam_A, cam_H, cam_V,
        msldemc_imxm, msldemc_imym);
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldemc_imxm);
    mxFree(msldemc_imym);
    mxFree(msldemc_imFOVmaskd);
    
}
