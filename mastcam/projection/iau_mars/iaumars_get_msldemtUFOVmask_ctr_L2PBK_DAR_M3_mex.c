/* =====================================================================
 * iaumars_get_msldemtUFOVmask_ctr_L2PBK_DAR_M3_mex.c
 * L2  : msldemc will be read from a file not an input.
 * PBK : Prior Binning into bins with the auxiliary size defined by two parameters K_L and K_S
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
 * 10 K_L                  double reciprocal of the length of the bin in the image line direction.
 * 11 K_S                  double reciprocal of the length of the bin in the image sample direction.
 * 12 dyu                  (int8_t) dynamic array masking flag
 *
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

#include <stdlib.h>
#include "envi.h"
#include "mex_create_array.h"
#include "cahvor.h"
#include "lib_proj_mastcamMSLDEM_IAUMars_L2PBK_DAR_M3.h"

void bin_msldemt_ctr_iaumars_L2PBK_DAR(double S_im, double L_im, CAHV_MODEL cahv_mdl,
        char *msldem_imgpath, EnviHeader msldem_hdr, double mslrad_offset,
        int32_T msldemc_imxy_sample_offset, int32_T msldemc_imxy_line_offset,
        int32_T msldemc_samples, int32_T msldemc_lines,
        double *msldemc_latitude, double *msldemc_longitude, 
        int8_T **msldemc_imFOVmask, 
        int32_T **bin_count_im,int32_T ***bin_im_c, int32_T ***bin_im_l,
        double ***bin_rad,
        double K_L, double K_S,
        int32_T *count_napmc, int32_T **c_napmc, int32_T **l_napmc,
        double **rad_napmc)
{
    int32_T xi,yi,c,l;
    int32_t binL,binS;
    int32_t binLm1,binSm1;
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
    
    binL = (int32_t) (K_L * L_im); binS = (int32_t) (K_S * S_im);
    binLm1 = binL-1; binSm1 = binS-1;
    
    cam_C = cahv_mdl.C; cam_A = cahv_mdl.A; cam_H = cahv_mdl.H; cam_V = cahv_mdl.V;
    
    *count_napmc = 0;
    
    cos_lon = (double*) malloc(sizeof(double) * (size_t) msldemc_samples);
    sin_lon = (double*) malloc(sizeof(double) * (size_t) msldemc_samples);
    for(c=0;c<msldemc_samples;c++){
        cos_lon[c] = cos(msldemc_longitude[c]);
        sin_lon[c] = sin(msldemc_longitude[c]);
    }
    
    /* Initialize the counting matrix */
    for(xi=0;xi<binS;xi++){
        for(yi=0;yi<binL;yi++){
            bin_count_im[xi][yi] = 0;
        }
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
                radius_tmp = (double) elevl[c] + mslrad_offset;
                x_iaumars  = radius_tmp * cos_latl * cos_lon[c];
                y_iaumars  = radius_tmp * cos_latl * sin_lon[c];
                z_iaumars  = radius_tmp * sin_latl;
                pmcx = x_iaumars - cam_C[0];
                pmcy = y_iaumars - cam_C[1];
                pmcz = z_iaumars - cam_C[2];
                
                apmc =  pmcx*cam_A[0] + pmcy*cam_A[1] + pmcz*cam_A[2];
                ppvx = (pmcx*cam_H[0] + pmcy*cam_H[1] + pmcz*cam_H[2])/apmc;
                ppvy = (pmcx*cam_V[0] + pmcy*cam_V[1] + pmcz*cam_V[2])/apmc;
                xi = (int32_T) floor(K_S*(ppvx+0.5));
                yi = (int32_T) floor(K_L*(ppvy+0.5));
                //xi = xi<0?0:xi;
                //xi = xi>S_imm1?:xi;
                if(xi<0)
                    xi=0;
                else if(xi>binSm1)
                    xi = binSm1;

                if(yi<0)
                    yi=0;
                else if(yi>binLm1)
                    yi = binLm1;

                ++bin_count_im[xi][yi];
            } else if(msldemc_imFOVmask[c][l]==1) {
                /*  */
                (*count_napmc)++;
            }
        }
    }
    
    for(xi=0;xi<binS;xi++){
        for(yi=0;yi<binL;yi++){
            if(bin_count_im[xi][yi]>0){
                bin_im_c[xi][yi] = (int32_T*) malloc(bin_count_im[xi][yi]*sizeof(int32_T));
                bin_im_l[xi][yi] = (int32_T*) malloc(bin_count_im[xi][yi]*sizeof(int32_T));
                bin_rad[xi][yi] = (double*) malloc(bin_count_im[xi][yi]*sizeof(double));
            } else {
                bin_im_c[xi][yi] = NULL;
                bin_im_l[xi][yi] = NULL;
                bin_rad[xi][yi] = NULL;
            }
        }
    }
    
    /* */
    if(*count_napmc>0){
        *c_napmc = (int32_T*) malloc(sizeof(int32_T) * (size_t) *count_napmc);
        *l_napmc = (int32_T*) malloc(sizeof(int32_T) * (size_t) *count_napmc);
        *rad_napmc = (double*) malloc(sizeof(double) * (size_t) *count_napmc);
    } else {
        *c_napmc = NULL;
        *l_napmc = NULL;
        *rad_napmc = NULL;
    }
    *count_napmc = 0;
    
    for(xi=0;xi<binS;xi++){
        for(yi=0;yi<binL;yi++){
            bin_count_im[xi][yi] = 0;
        }
    }
    
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
                
                ppvx = (pmcx*cam_H[0] + pmcy*cam_H[1] + pmcz*cam_H[2])/apmc;
                ppvy = (pmcx*cam_V[0] + pmcy*cam_V[1] + pmcz*cam_V[2])/apmc;
                xi = (int32_T) floor(K_S*(ppvx+0.5));
                yi = (int32_T) floor(K_L*(ppvy+0.5));
                if(xi<0)
                    xi=0;
                else if(xi>binSm1)
                    xi = binSm1;

                if(yi<0)
                    yi=0;
                else if(yi>binLm1)
                    yi = binLm1;

                bin_im_c[xi][yi][bin_count_im[xi][yi]] = c;
                bin_im_l[xi][yi][bin_count_im[xi][yi]] = l;
                bin_rad[xi][yi][bin_count_im[xi][yi]]  = radius_tmp;
                ++bin_count_im[xi][yi];

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
    int8_t dyu;
    
    int8_T **msldemt_imUFOVmask;
    
    mwIndex si,li;
    
    int32_T **bin_count_im, *bin_count_im_base;
    int32_T ***bin_im_c, **bin_im_c_base;
    int32_T ***bin_im_l, **bin_im_l_base;
    double ***bin_rad, **bin_rad_base;
    int32_T *c_napmc, *l_napmc;
    int32_T count_napmc;
    double *rad_napmc;
    
    double K_L,K_S;
    mwSize binL,binS;

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
    
    /* INPUT 1 msldem_header and 2 radius offset */
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
    
    /* INPUT 10/11 bin size parameters */
    K_L = mxGetScalar(prhs[10]);
    K_S = mxGetScalar(prhs[11]);
    
    /* INPUT 12 dynamic masking flag */
    dyu = (int8_t) mxGetScalar(prhs[12]);
    
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
     * binning the image
     * ----------------------------------------------------------------- */
    binL = (mwSize) (K_L * L_im); binS = (mwSize) (K_S * S_im);
    createInt32Matrix(&bin_count_im,&bin_count_im_base,(size_t) binS, (size_t) binL);
    createInt32PMatrix(&bin_im_c, &bin_im_c_base, (size_t) binS, (size_t) binL);
    createInt32PMatrix(&bin_im_l, &bin_im_l_base, (size_t) binS, (size_t) binL);
    createDoublePMatrix(&bin_rad, &bin_rad_base, (size_t) binS, (size_t) binL);
    
    bin_msldemt_ctr_iaumars_L2PBK_DAR(S_im, L_im, cahv_mdl,
            msldem_imgpath, msldem_header, mslrad_offset,
            msldemc_imxy_sample_offset, msldemc_imxy_line_offset,
            msldemc_samples, msldemc_lines, 
            msldemc_latitude, msldemc_longitude, msldemc_imFOVmask,
            bin_count_im, bin_im_c, bin_im_l, bin_rad,
            K_L,K_S,
            &count_napmc,&c_napmc,&l_napmc,&rad_napmc);
    
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    if(dyu==1){
        mask_obstructed_pts_in_msldemt_using_msldemc_iaumars_L2PBK_DARDYM_M3(
                msldem_imgpath, msldem_header, mslrad_offset,
            (int32_T) msldemc_imxy_sample_offset, (int32_T) msldemc_imxy_line_offset,
            (int32_T) msldemc_samples, (int32_T) msldemc_lines,
            msldemc_latitude, msldemc_longitude, msldemc_imFOVmask,
            (int32_T) msldemc_samples, (int32_T) msldemc_lines,
            msldemc_latitude, msldemc_longitude, msldemt_imUFOVmask,
            bin_count_im, bin_im_c, bin_im_l,bin_rad,
            K_L,K_S,
            count_napmc, c_napmc, l_napmc, rad_napmc,
            S_im, L_im, cahv_mdl);
    } else if(dyu==0){
        mask_obstructed_pts_in_msldemt_using_msldemc_iaumars_L2PBK_DAR_M3(
                msldem_imgpath, msldem_header, mslrad_offset,
            (int32_T) msldemc_imxy_sample_offset, (int32_T) msldemc_imxy_line_offset,
            (int32_T) msldemc_samples, (int32_T) msldemc_lines,
            msldemc_latitude, msldemc_longitude, msldemc_imFOVmask,
            (int32_T) msldemc_samples, (int32_T) msldemc_lines,
            msldemc_latitude, msldemc_longitude, msldemt_imUFOVmask,
            bin_count_im, bin_im_c, bin_im_l,bin_rad,
            K_L,K_S,
            count_napmc, c_napmc, l_napmc, rad_napmc,
            S_im, L_im, cahv_mdl);
    } else {
        printf("Undefined DYU_FLAG=%d\n",dyu);
    }
        
    
    /* Freeing memories */
    for(si=0;si<binS;si++){
        for(li=0;li<binL;li++){
            if(bin_count_im[si][li]>0){
                free(bin_im_c[si][li]);
                free(bin_im_l[si][li]);
                free(bin_rad[si][li]);
            }
        }
    }
    
    free(bin_im_c);
    free(bin_im_l);
    free(bin_im_c_base);
    free(bin_im_l_base);
    
    free(bin_rad);
    free(bin_rad_base);
    
    
    // mxFree(bin_count_im);
    free(bin_count_im);
    free(bin_count_im_base);
    
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