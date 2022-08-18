/* =====================================================================
 * iaumars_get_msldemtUFOVmask_d_wmsldemc_L2K_mex.c
 * Evaluate any pixels in MSLDEM image whether or not they exist in the 
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
 *
 * 
 * 
 * OUTPUTS:
 * 0  msldemt_inImage     [L_demc-1 x S_demc-1]  Boolean
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
#include "lib_proj_mastcamMSLDEM_L2K_IAUMars.h"

void bin_msldemt_xyz_wAHVint_L2K_d(double S_im, double L_im, CAHV_MODEL cahv_mdl,
        char *msldem_imgpath, EnviHeader msldem_hdr,  double mslrad_offset,
        int32_T msldemc_imxy_sample_offset, int32_T msldemc_imxy_line_offset,
        int32_T msldemc_samples, int32_T msldemc_lines,
        double *msldemc_latitude, double *msldemc_longitude, int8_T **msldemc_imFOVmask,
        double *msldemt_latitude, double *msldemt_longitude, int8_T **msldemt_imFOVmask,
        int32_T **bin_count_im,int32_T ***bin_im_c, int32_T ***bin_im_l,
        double ***bin_imx, double ***bin_imy,double ***bin_rad,
        double K_L, double K_S,
        int32_T *count_napmc, int32_T **c_napmc, int32_T **l_napmc, double **rad_napmc)
{
    int32_T xi,yi,c,l;
    int32_t binL,binS;
    int32_t binLm1,binSm1;
    double pmcx,pmcy,pmcz,apmc,ppvx,ppvy;
    double apmcx,hpmcx,vpmcx;
    double *cam_C, *cam_A, *cam_H, *cam_V;
    
    long skip_pri;
    long skip_l, skip_r;
    float *elevl,*elevlp1;
    size_t ncpy;
    size_t sz=sizeof(float);
    FILE *fid;
    float data_ignore_value_float;
    double dem_cl;
    int32_T msldemc_samplesm1, msldemc_linesm1;
    
    double *cos_tlon, *sin_tlon;
    double tlonc, tlatl;
    double radius_tmp;
    double cos_tlatl, sin_tlatl;
    double x_iaumars, y_iaumars, z_iaumars;
    
    msldemc_samplesm1 = msldemc_samples - 1;
    msldemc_linesm1   = msldemc_lines - 1;
    
    binL = (int32_t) (K_L * L_im); binS = (int32_t) (K_S * S_im);
    binLm1 = binL-1; binSm1 = binS-1;
    
    cam_C = cahv_mdl.C; cam_A = cahv_mdl.A;
    cam_H = cahv_mdl.H; cam_V = cahv_mdl.V;
    
    *count_napmc = 0;
    
    /* Initialize the counting matrix */
    for(xi=0;xi<binS;xi++){
        for(yi=0;yi<binL;yi++){
            bin_count_im[xi][yi] = 0;
        }
    }
    
    cos_tlon = (double*) malloc(sizeof(double) * (size_t) msldemc_samplesm1);
    sin_tlon = (double*) malloc(sizeof(double) * (size_t) msldemc_samplesm1);
    for(c=0;c<msldemc_samplesm1;c++){
        tlonc = 0.5*(msldemc_longitude[c]+msldemc_longitude[c+1]);
        cos_tlon[c] = cos(tlonc);
        sin_tlon[c] = sin(tlonc);
        msldemt_longitude[c] = tlonc;
    }
    
    
    
    fid = fopen(msldem_imgpath,"rb");
    /* skip lines */
    skip_pri = (long) msldem_hdr.samples * (long) msldemc_imxy_line_offset * (long) sz;
    // printf("%d*%d*%d=%ld\n",msldem_header.samples,msldemc_imxy_line_offset,s,skip_pri);
    fseek(fid,skip_pri,SEEK_CUR);
    /* read the data */
    ncpy = sz * (size_t) msldemc_samples;
    elevl = (float*) malloc(ncpy);
    elevlp1 = (float*) malloc(ncpy);
    skip_l = (long) sz * (long) msldemc_imxy_sample_offset;
    skip_r = ((long) msldem_hdr.samples - (long) msldemc_samples)* (long) sz - skip_l;
    data_ignore_value_float = (float) msldem_hdr.data_ignore_value + 1.0;
    
    fseek(fid,skip_l,SEEK_CUR);
    fread(elevlp1,sz,msldemc_samples,fid);
    fseek(fid,skip_r,SEEK_CUR);
    for(c=0;c<msldemc_samples;c++){
        if(elevlp1[c]<data_ignore_value_float)
            elevlp1[c] = NAN;
    }
    
    
    // printf("%d\n",msldemc_samples);
    for(l=0;l<msldemc_linesm1;l++){
        /* read a line */ 
        memcpy(elevl,elevlp1,ncpy);
        
        fseek(fid,skip_l,SEEK_CUR);
        fread(elevlp1,sz,msldemc_samples,fid);
        fseek(fid,skip_r,SEEK_CUR);
        for(c=0;c<msldemc_samples;c++){
            if(elevlp1[c]<data_ignore_value_float)
                elevlp1[c] = NAN;
        }
        
        tlatl = 0.5*(msldemc_latitude[l]+msldemc_latitude[l+1]);
        cos_tlatl = cos(tlatl);
        sin_tlatl = sin(tlatl);
        msldemt_latitude[l] = tlatl;
        for(c=0;c<msldemc_samplesm1;c++){
            if((msldemc_imFOVmask[c][l]>0) || (msldemc_imFOVmask[c+1][l]>0) || 
                    (msldemc_imFOVmask[c][l+1]>0) || (msldemc_imFOVmask[c+1][l+1]>0)){
                radius_tmp = 0.5*( (double) elevl[c+1] + (double) elevlp1[c] ) + mslrad_offset;
                x_iaumars  = radius_tmp * cos_tlatl * cos_tlon[c];
                y_iaumars  = radius_tmp * cos_tlatl * sin_tlon[c];
                z_iaumars  = radius_tmp * sin_tlatl;
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
                
                if(apmc>0){
                    ++bin_count_im[xi][yi];
                    msldemt_imFOVmask[c][l] = 2;
                } else if(apmc<0) {
                    (*count_napmc)++;
                    msldemt_imFOVmask[c][l] = 1;
                }
            } 
        }
    }
    
    for(xi=0;xi<binS;xi++){
        for(yi=0;yi<binL;yi++){
            if(bin_count_im[xi][yi]>0){
                bin_im_c[xi][yi] = (int32_T*) malloc(bin_count_im[xi][yi]*sizeof(int32_T));
                bin_im_l[xi][yi] = (int32_T*) malloc(bin_count_im[xi][yi]*sizeof(int32_T));
                bin_imx[xi][yi] = (double*) malloc(bin_count_im[xi][yi]*sizeof(double));
                bin_imy[xi][yi] = (double*) malloc(bin_count_im[xi][yi]*sizeof(double));
                bin_rad[xi][yi] = (double*) malloc(bin_count_im[xi][yi]*sizeof(double));
            } else {
                bin_im_c[xi][yi] = NULL;
                bin_im_l[xi][yi] = NULL;
                bin_imx[xi][yi]  = NULL;
                bin_imy[xi][yi]  = NULL;
                bin_rad[xi][yi]  = NULL;
            }
        }
    }
    
    /* */
    if(*count_napmc>0){
        *c_napmc = (int32_T*) malloc(sizeof(int32_T) * (size_t) *count_napmc);
        *l_napmc = (int32_T*) malloc(sizeof(int32_T) * (size_t) *count_napmc);
        *rad_napmc = (double*) malloc(sizeof(double) * (size_t) *count_napmc);
    } else {
        *c_napmc   = NULL;
        *l_napmc   = NULL;
        *rad_napmc = NULL;
    }
    *count_napmc = 0;
    
    for(xi=0;xi<S_im;xi++){
        for(yi=0;yi<L_im;yi++){
            bin_count_im[xi][yi] = 0;
        }
    }
    
    fseek(fid,skip_pri,SEEK_SET);
    fseek(fid,skip_l,SEEK_CUR);
    fread(elevlp1,sz,msldemc_samples,fid);
    fseek(fid,skip_r,SEEK_CUR);
    for(c=0;c<msldemc_samples;c++){
        if(elevlp1[c]<data_ignore_value_float)
            elevlp1[c] = NAN;
    }
    
    for(l=0;l<msldemc_linesm1;l++){
        memcpy(elevl,elevlp1,ncpy);
        
        fseek(fid,skip_l,SEEK_CUR);
        fread(elevlp1,sz,msldemc_samples,fid);
        fseek(fid,skip_r,SEEK_CUR);
        for(c=0;c<msldemc_samples;c++){
            if(elevlp1[c]<data_ignore_value_float)
                elevlp1[c] = NAN;
        }
        tlatl = 0.5*(msldemc_latitude[l]+msldemc_latitude[l+1]);
        cos_tlatl = cos(tlatl);
        sin_tlatl = sin(tlatl);
        for(c=0;c<msldemc_samplesm1;c++){
            if(msldemt_imFOVmask[c][l]>1){
                radius_tmp = 0.5*( (double) elevl[c+1] + (double) elevlp1[c] ) + mslrad_offset;
                x_iaumars  = radius_tmp * cos_tlatl * cos_tlon[c];
                y_iaumars  = radius_tmp * cos_tlatl * sin_tlon[c];
                z_iaumars  = radius_tmp * sin_tlatl;
                pmcx = x_iaumars - cam_C[0];
                pmcy = y_iaumars - cam_C[1];
                pmcz = z_iaumars - cam_C[2];
                apmc =  pmcx*cam_A[0] + pmcy*cam_A[1] + pmcz*cam_A[2];
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
                bin_imx[xi][yi][bin_count_im[xi][yi]]  = ppvx;
                bin_imy[xi][yi][bin_count_im[xi][yi]]  = ppvy;
                bin_rad[xi][yi][bin_count_im[xi][yi]]  = radius_tmp;
                ++bin_count_im[xi][yi];

            } else if(msldemt_imFOVmask[c][l]==1){
                radius_tmp = 0.5*( (double) elevl[c+1] + (double) elevlp1[c] ) + mslrad_offset;
                (*c_napmc)[*count_napmc] = c;
                (*l_napmc)[*count_napmc] = l;
                (*rad_napmc)[*count_napmc] = radius_tmp;
                (*count_napmc)++;
            }
        }
    }
    
    fclose(fid);
    free(elevl);
    free(elevlp1);
    free(cos_tlon);
    free(sin_tlon);
}
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *msldem_imgpath;
    EnviHeader msldem_header;
    double mslrad_offset;
    mwSize msldemc_imxy_sample_offset,msldemc_imxy_line_offset;
    mwSize msldemc_samples,msldemc_lines;
    double *msldemc_latitude;
    double *msldemc_longitude;
    CAHV_MODEL cahv_mdl;
    int8_T **msldemc_imFOVmask;
    double S_im,L_im;
    
    mwSize msldemt_samples, msldemt_lines;
    double *msldemt_latitude;
    double *msldemt_longitude;
    int8_T **msldemt_imUFOVmask;
    
    mwIndex si,li;
    
    int32_T **bin_count_im, *bin_count_im_base;
    int32_T ***bin_im_c, **bin_im_c_base;
    int32_T ***bin_im_l, **bin_im_l_base;
    double ***bin_imx, **bin_imx_base;
    double ***bin_imy, **bin_imy_base;
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
    
    /* INPUT 1 msldem_header */
    msldem_header = mxGetEnviHeader(prhs[1]);
    
    mslrad_offset = mxGetScalar(prhs[2]);
    
    /* INPUT 2 msldemc_sheader*/
    msldemc_imxy_sample_offset = (mwSize) mxGetScalar(mxGetField(prhs[3],0,"sample_offset"));
    msldemc_imxy_line_offset = (mwSize) mxGetScalar(mxGetField(prhs[3],0,"line_offset"));
    msldemc_samples = (mwSize) mxGetScalar(mxGetField(prhs[3],0,"samples"));
    msldemc_lines = (mwSize) mxGetScalar(mxGetField(prhs[3],0,"lines"));
    //msldemc_img = set_mxDoubleMatrix(prhs[0]);
    //L_demc = mxGetM(prhs[0]);
    //S_demc = mxGetN(prhs[0]);
    
    /* INPUT 1/2 msldem northing easting */
    msldemc_latitude  = mxGetDoubles(prhs[4]);
    msldemc_longitude = mxGetDoubles(prhs[5]);
    
    /* INPUT 3 msldemc imFOV */
    msldemc_imFOVmask = set_mxInt8Matrix(prhs[6]);
    
    /* INPUT 4/5 image S_im, L_im */
    S_im = mxGetScalar(prhs[7]);
    L_im = mxGetScalar(prhs[8]);
    
    /* INPUT 6 CAHV model */
    cahv_mdl = mxGet_CAHV_MODEL(prhs[9]);
    
    K_L = mxGetScalar(prhs[10]);
    K_S = mxGetScalar(prhs[11]);
    
    
    /* OUTPUT 0 msldemc imFOV */
    msldemt_samples = msldemc_samples-1;
    msldemt_lines = msldemc_lines-1;
    plhs[0] = mxCreateNumericMatrix(msldemt_lines,msldemt_samples,
            mxINT8_CLASS,mxREAL);
    msldemt_imUFOVmask = set_mxInt8Matrix(plhs[0]);
    
    // Initialize matrices
    for(si=0;si<msldemt_samples;si++){
        for(li=0;li<msldemt_lines;li++){
            msldemt_imUFOVmask[si][li] = 0;
        }
    }
    
    /* -----------------------------------------------------------------
     * binning the image
     * ----------------------------------------------------------------- */
    binL = (mwSize) (K_L * L_im); binS = (mwSize) (K_S * S_im);
    createInt32Matrix(&bin_count_im,&bin_count_im_base,(size_t) binS, (size_t) binL);
    createInt32PMatrix(&bin_im_c, &bin_im_c_base, (size_t) binS, (size_t) binL);
    createInt32PMatrix(&bin_im_l, &bin_im_l_base, (size_t) binS, (size_t) binL);
    createDoublePMatrix(&bin_imx, &bin_imx_base, (size_t) binS, (size_t) binL);
    createDoublePMatrix(&bin_imy, &bin_imy_base, (size_t) binS, (size_t) binL);
    createDoublePMatrix(&bin_rad, &bin_rad_base, (size_t) binS, (size_t) binL);
    
    msldemt_latitude  = (double*) malloc(sizeof(double)*(size_t) msldemt_lines);
    msldemt_longitude = (double*) malloc(sizeof(double)*(size_t) msldemt_samples);
    
    bin_msldemt_xyz_wAHVint_L2K_d(S_im, L_im, cahv_mdl,
            msldem_imgpath, msldem_header, mslrad_offset,
            msldemc_imxy_sample_offset, msldemc_imxy_line_offset,
            msldemc_samples, msldemc_lines, 
            msldemc_latitude, msldemc_longitude, msldemc_imFOVmask,
            msldemt_latitude, msldemt_longitude, msldemt_imUFOVmask,
            bin_count_im, bin_im_c, bin_im_l, bin_imx, bin_imy, bin_rad,
            K_L,K_S,
            &count_napmc,&c_napmc,&l_napmc,&rad_napmc);
    
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */    
    mask_obstructed_pts_in_msldemt_using_msldemc_L2K_iaumars(
        msldem_imgpath, msldem_header, mslrad_offset,
        (int32_T) msldemc_imxy_sample_offset, (int32_T) msldemc_imxy_line_offset,
        (int32_T) msldemc_samples, (int32_T) msldemc_lines,
        msldemc_latitude, msldemc_longitude, msldemc_imFOVmask,
        (int32_t) msldemt_samples, (int32_t) msldemt_lines,
        msldemt_latitude, msldemt_longitude, msldemt_imUFOVmask,
        bin_count_im, bin_im_c, bin_im_l,
        bin_imx, bin_imy, bin_rad,
        K_L,K_S,
        count_napmc, c_napmc, l_napmc, rad_napmc,
        S_im, L_im, cahv_mdl);
        
    
    /* Freeing memories */
    for(si=0;si<binS;si++){
        for(li=0;li<binL;li++){
            if(bin_count_im[si][li]>0){
                free(bin_im_c[si][li]);
                free(bin_im_l[si][li]);
                free(bin_imx[si][li]);
                free(bin_imy[si][li]);
                free(bin_rad[si][li]);
            }
        }
    }
    
    free(bin_im_c);
    free(bin_im_l);
    free(bin_im_c_base);
    free(bin_im_l_base);
    
    free(bin_imx);
    free(bin_imy);
    free(bin_imx_base);
    free(bin_imy_base);
    free(bin_rad);
    free(bin_rad_base);
    
    free(bin_count_im);
    free(bin_count_im_base);
    
    if(count_napmc>0){
        free(c_napmc);
        free(l_napmc);
        free(rad_napmc);
    }
    
    free(msldemt_latitude);
    free(msldemt_longitude);
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldemc_imFOVmask);
    mxFree(msldemt_imUFOVmask);
    
}