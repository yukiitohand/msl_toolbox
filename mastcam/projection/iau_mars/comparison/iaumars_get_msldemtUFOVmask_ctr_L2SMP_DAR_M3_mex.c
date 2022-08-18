/* =====================================================================
 * iaumars_get_msldemtUFOVmask_ctr_L2SMP_DAR_M3_mex.c
 * L2  : Library type 2 msldemc will be read from a file not an input. 
 * SMP : SiMPle iteration without prior binning 
 * DAR : Dynamic numeric ARray (c,l,radius)
 * M3  : 3x3 matrix inversion object image coordinate.
 * 
 * INPUTS:
 * 0 msldem_imgpath        char* path to the image
 * 1 msldem_hdr            EnviHeader
 * 2 mslrad_offset         Double Scalar
 * 3 msldemc_imFOVhdr      Struct
 * 4 msldemc_latitude      Double array [L_demc]
 * 5 msldemc_longitude     Double array [S_demc]
 * 6 msldemc_imFOVmask     int8_T [L_demc x S_demc]
 * 7 S_im                  int
 * 8 L_im                  int
 * 9 cahv_mdl              CAHV_MODEL
 * 10 dym                  int8_t dynamic array masking
 * 
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
#include <time.h>

#include <stdlib.h>
#include "envi.h"
#include "mex_create_array.h"
#include "cahvor.h"
#include "lib_proj_mastcamMSLDEM_IAUMars_L2SMP_DAR_M3.h"

void create_SMPDAR_iaumars_ctr(double S_im, double L_im, CAHV_MODEL cahv_mdl,
        char *msldem_imgpath, EnviHeader msldem_hdr, double mslrad_offset,
        int32_T msldemc_imxy_sample_offset, int32_T msldemc_imxy_line_offset,
        int32_T msldemc_samples, int32_T msldemc_lines,
        double *msldemc_latitude, double *msldemc_longitude, 
        int8_T **msldemc_imFOVmask, 
        int64_t *count_papmc, int32_T **c_papmc, int32_T **l_papmc,
        double **rad_papmc,
        int32_T *count_napmc, int32_T **c_napmc, int32_T **l_napmc,
        double **rad_napmc)
{
    int32_T c,l;
    double pmcx,pmcy,pmcz,apmc,ppvx,ppvy;
    double *cam_C, *cam_A, *cam_H, *cam_V;
    
    long skip_pri;
    long skip_l, skip_r;
    float *elevl;
    size_t ncpy;
    size_t sz=sizeof(float);
    FILE *fid;
    float data_ignore_value_float;
    // double dem_cl;
    
    double *cos_lon, *sin_lon;
    double radius_tmp;
    double cos_latl, sin_latl;
    double x_iaumars, y_iaumars, z_iaumars;
    
    // S_imm1 = S_im - 1;
    // L_imm1 = L_im - 1;
    
    
    cam_C = cahv_mdl.C; cam_A = cahv_mdl.A; cam_H = cahv_mdl.H; cam_V = cahv_mdl.V;
    
    *count_papmc = 0;
    *count_napmc = 0;
    
    cos_lon = (double*) malloc(sizeof(double) * (size_t) msldemc_samples);
    sin_lon = (double*) malloc(sizeof(double) * (size_t) msldemc_samples);
    for(c=0;c<msldemc_samples;c++){
        cos_lon[c] = cos(msldemc_longitude[c]);
        sin_lon[c] = sin(msldemc_longitude[c]);
    }
    
    /* Initialize the counting matrix */
    
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
    data_ignore_value_float = (float) msldem_hdr.data_ignore_value + 1.0;
    
    // printf("%d\n",msldemc_samples);
    for(l=0;l<msldemc_lines;l++){
        /* read a line */ 
        fseek(fid,skip_l,SEEK_CUR);
        fread(elevl,sz,msldemc_samples,fid);
        fseek(fid,skip_r,SEEK_CUR);
        for(c=0;c<msldemc_samples;c++){
            if(elevl[c]<data_ignore_value_float)
                elevl[c] = NAN;
        }
        cos_latl   = cos(msldemc_latitude[l]);
        sin_latl   = sin(msldemc_latitude[l]);
        for(c=0;c<msldemc_samples;c++){
            if((msldemc_imFOVmask[c][l]>1)){
                (*count_papmc)++;
                
            } else if(msldemc_imFOVmask[c][l]==1) {
                /*  */
                (*count_napmc)++;
            }
        }
    }
    
    
    /* */
     if(*count_papmc>0){
        *c_papmc   = (int32_T*) malloc(sizeof(int32_T) * (size_t) *count_papmc);
        *l_papmc   = (int32_T*) malloc(sizeof(int32_T) * (size_t) *count_papmc);
        *rad_papmc = (double*)  malloc(sizeof(double)  * (size_t) *count_papmc);
    } else {
        *c_papmc   = NULL;
        *l_papmc   = NULL;
        *rad_papmc = NULL;
    }
    *count_papmc = 0;
    if(*count_napmc>0){
        *c_napmc   = (int32_T*) malloc(sizeof(int32_T) * (size_t) *count_napmc);
        *l_napmc   = (int32_T*) malloc(sizeof(int32_T) * (size_t) *count_napmc);
        *rad_napmc = (double*)  malloc(sizeof(double)  * (size_t) *count_napmc);
    } else {
        *c_napmc   = NULL;
        *l_napmc   = NULL;
        *rad_napmc = NULL;
    }
    *count_napmc = 0;
    
    
    printf("a\n");
    
    fseek(fid,skip_pri,SEEK_SET);
    
    for(l=0;l<msldemc_lines;l++){
        fseek(fid,skip_l,SEEK_CUR);
        fread(elevl,sz,msldemc_samples,fid);
        fseek(fid,skip_r,SEEK_CUR);
        for(c=0;c<msldemc_samples;c++){
            if(elevl[c]<data_ignore_value_float)
                elevl[c] = NAN;
        }
        cos_latl   = cos(msldemc_latitude[l]);
        sin_latl   = sin(msldemc_latitude[l]);
        for(c=0;c<msldemc_samples;c++){
            if(msldemc_imFOVmask[c][l]>1){
                radius_tmp = (double) elevl[c] + mslrad_offset;
                x_iaumars  = radius_tmp * cos_latl * cos_lon[c];
                y_iaumars  = radius_tmp * cos_latl * sin_lon[c];
                z_iaumars  = radius_tmp * sin_latl;
                pmcx = x_iaumars - cam_C[0];
                pmcy = y_iaumars - cam_C[1];
                pmcz = z_iaumars - cam_C[2];
                
                apmc   = pmcx*cam_A[0] + pmcy*cam_A[1] + pmcz*cam_A[2];
                
                //ppvx = (pmcx*cam_H[0] + pmcy*cam_H[1] + pmcz*cam_H[2])/apmc;
                //ppvy = (pmcx*cam_V[0] + pmcy*cam_V[1] + pmcz*cam_V[2])/apmc;

                (*c_papmc)[*count_papmc] = c;
                (*l_papmc)[*count_papmc] = l;
                (*rad_papmc)[*count_papmc] = radius_tmp;
                (*count_papmc)++;

            } else if(msldemc_imFOVmask[c][l]==1){
                radius_tmp = (double) elevl[c] + mslrad_offset;
                (*c_napmc)[*count_napmc] = c;
                (*l_napmc)[*count_napmc] = l;
                (*rad_napmc)[*count_napmc] = radius_tmp;
                (*count_napmc)++;
            }
        }
    }
    
    fclose(fid);
    free(elevl);
    free(cos_lon);
    free(sin_lon);
}
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *msldem_imgpath;
    EnviHeader msldem_header;
    mwSize msldemc_imxy_sample_offset,msldemc_imxy_line_offset;
    mwSize msldemc_samples,msldemc_lines;
    double *msldemc_latitude;
    double *msldemc_longitude;
    double mslrad_offset;
    CAHV_MODEL cahv_mdl;
    int8_T **msldemc_imFOVmask;
    double S_im,L_im;
    int8_t dym;
    
    int8_T **msldemt_imUFOVmask;
    
    mwIndex si,li;
    
    
    int32_t *c_papmc, *l_papmc;
    int64_t count_papmc;
    double *rad_papmc;
    
    int32_T *c_napmc, *l_napmc;
    int32_T count_napmc;
    double *rad_napmc;
    
    clock_t strt_time, end_time;
    double cpu_time_used;
    
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
    
    /* INPUT 1 msldem_header and 2 offset */
    msldem_header = mxGetEnviHeader(prhs[1]);
    
    mslrad_offset     = mxGetScalar(prhs[2]);
    
    /* INPUT 3 msldemc_sheader*/
    msldemc_imxy_sample_offset = (mwSize) mxGetScalar(mxGetField(prhs[3],0,"sample_offset"));
    msldemc_imxy_line_offset = (mwSize) mxGetScalar(mxGetField(prhs[3],0,"line_offset"));
    msldemc_samples = (mwSize) mxGetScalar(mxGetField(prhs[3],0,"samples"));
    msldemc_lines = (mwSize) mxGetScalar(mxGetField(prhs[3],0,"lines"));
    //msldemc_img = set_mxDoubleMatrix(prhs[0]);
    //L_demc = mxGetM(prhs[0]);
    //S_demc = mxGetN(prhs[0]);
    
    /* INPUT 4/5 msldem northing easting */
    msldemc_latitude  = mxGetDoubles(prhs[4]);
    msldemc_longitude = mxGetDoubles(prhs[5]);
    
    
    /* INPUT 6 msldemc imFOV */
    msldemc_imFOVmask = set_mxInt8Matrix(prhs[6]);
    
    /* INPUT 7/8 image S_im, L_im */
    S_im = mxGetScalar(prhs[7]);
    L_im = mxGetScalar(prhs[8]);
    
    /* INPUT 9 CAHV model */
    cahv_mdl = mxGet_CAHV_MODEL(prhs[9]);
    
    /* INPUT 10 Dynamic Array Masking FLAG */
    dym = (int8_t) mxGetScalar(prhs[10]);
    
    /* OUTPUT 0 msldemc imFOV */
    plhs[0] = mxCreateNumericMatrix(msldemc_lines,msldemc_samples,mxINT8_CLASS,mxREAL);
    msldemt_imUFOVmask = set_mxInt8Matrix(plhs[0]);
    
    // plhs[1] = mxCreateNumericMatrix(L_im,S_im,mxINT32_CLASS,mxREAL);
    // bin_count_im = set_mxInt32Matrix(plhs[1]);
    
    // Initialize matrices
    for(si=0;si<msldemc_samples;si++){
        for(li=0;li<msldemc_lines;li++){
            if(msldemc_imFOVmask[si][li]>1){
                msldemt_imUFOVmask[si][li] = 2;
            } else {
                msldemt_imUFOVmask[si][li] = msldemc_imFOVmask[si][li];
            }
        }
    }
    
    /* -----------------------------------------------------------------
     * Create Simple dynamic array for points to be tested.
     * ----------------------------------------------------------------- */
    strt_time = clock();
    create_SMPDAR_iaumars_ctr(S_im, L_im, cahv_mdl,
            msldem_imgpath, msldem_header, mslrad_offset,
            msldemc_imxy_sample_offset, msldemc_imxy_line_offset,
            msldemc_samples, msldemc_lines, 
            msldemc_latitude, msldemc_longitude, msldemc_imFOVmask,
            &count_papmc,&c_papmc,&l_papmc,&rad_papmc,
            &count_napmc,&c_napmc,&l_napmc,&rad_napmc);
    
    end_time = clock();
    cpu_time_used = ((double) (end_time - strt_time)) / CLOCKS_PER_SEC;
    printf("Preprocessing took %f [sec].\n",cpu_time_used);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    strt_time = clock();
    if(dym==1){
        mask_obstructed_pts_in_msldemt_using_msldemc_iaumars_L2SMP_DARDYM_M3(
                msldem_imgpath, msldem_header, mslrad_offset,
            (int32_T) msldemc_imxy_sample_offset, (int32_T) msldemc_imxy_line_offset,
            (int32_T) msldemc_samples, (int32_T) msldemc_lines,
            msldemc_latitude, msldemc_longitude, msldemc_imFOVmask,
            (int32_T) msldemc_samples, (int32_T) msldemc_lines,
            msldemc_latitude, msldemc_longitude, msldemt_imUFOVmask,
            count_papmc, c_papmc, l_papmc, rad_papmc,
            count_napmc, c_napmc, l_napmc, rad_napmc,
            S_im, L_im, cahv_mdl);
    } else if(dym==0){
        printf("Not printed for DYM_FLAG=%d\n",dym);
    } else {
        printf("Undefined DYM_FLAG=%d\n",dym);
    }
    
    end_time = clock();
    cpu_time_used = ((double) (end_time - strt_time)) / CLOCKS_PER_SEC;
    printf("Main operation took %f [sec].\n",cpu_time_used);
        
    
    /* Freeing memories */
    
    if(count_papmc>0){
        free(c_papmc);
        free(l_papmc);
        free(rad_papmc);
    }
    
    if(count_napmc>0){
        free(c_napmc);
        free(l_napmc);
        free(rad_napmc);
    }
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldemc_imFOVmask);
    mxFree(msldemt_imUFOVmask);
    
}