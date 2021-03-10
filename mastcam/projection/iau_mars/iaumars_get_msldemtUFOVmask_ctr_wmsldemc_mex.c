/* =====================================================================
 * iaumars_get_msldemtUFOVmask_ctr_wmsldemc_mex.c
 * Evaluate any pixels in MSLDEM image whether or not they exist in the 
 * 
 * INPUTS:
 * 0 msldem_img            char* path to the image
 * 4 msldemc_latitude      Double array [L_demc]
 * 5 msldemc_longitude     Double array [S_demc]
 * 6 msldemc_imFOVmask     int8_T [L_demc x S_demc]
 * 7 S_im                  int
 * 8 L_im                  int
 * 9 cahv_mdl              CAHV_MODEL
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
#include "lib_proj_mastcamMSLDEM_IAUMars.h"

void bin_msldemt_xyz_wAHVint_iaumars_ctr(int32_T S_im, int32_T L_im, CAHV_MODEL cahv_mdl,
        double **msldemc_img,
        int32_T msldemc_samples, int32_T msldemc_lines,
        double *msldemc_latitude, double *msldemc_longitude, 
        int8_T **msldemc_imFOVmask, 
        int32_T **bin_count_im,int32_T ***bin_im_c, int32_T ***bin_im_l,
        double ***bin_imx, double ***bin_imy,double ***bin_rad,
        int32_T *count_napmc, int32_T **c_napmc, int32_T **l_napmc,
        double **rad_napmc)
{
    int32_T xi,yi,c,l;
    int32_T S_imm1,L_imm1;
    double pmcx,pmcy,pmcz,apmc,ppvx,ppvy;
    double *cam_C, *cam_A, *cam_H, *cam_V;
    
    double *cos_lon, *sin_lon;
    double radius_tmp;
    double cos_latl, sin_latl;
    double x_iaumars, y_iaumars, z_iaumars;
    
    S_imm1 = S_im - 1;
    L_imm1 = L_im - 1;
    
    cam_C = cahv_mdl.C; cam_A = cahv_mdl.A; cam_H = cahv_mdl.H; cam_V = cahv_mdl.V;
    
    *count_napmc = 0;
    
    cos_lon = (double*) malloc(sizeof(double) * (size_t) msldemc_samples);
    sin_lon = (double*) malloc(sizeof(double) * (size_t) msldemc_samples);
    for(c=0;c<msldemc_samples;c++){
        cos_lon[c] = cos(msldemc_longitude[c]);
        sin_lon[c] = sin(msldemc_longitude[c]);
    }
    
    /* Initialize the counting matrix */
    for(xi=0;xi<S_im;xi++){
        for(yi=0;yi<L_im;yi++){
            bin_count_im[xi][yi] = 0;
        }
    }
    
    // Count the number of points that fall into each bin.
    for(l=0;l<msldemc_lines;l++){
        cos_latl   = cos(msldemc_latitude[l]);
        sin_latl   = sin(msldemc_latitude[l]);
        for(c=0;c<msldemc_samples;c++){
            if((msldemc_imFOVmask[c][l]>1)){
                radius_tmp = msldemc_img[c][l];
                x_iaumars  = radius_tmp * cos_latl * cos_lon[c];
                y_iaumars  = radius_tmp * cos_latl * sin_lon[c];
                z_iaumars  = radius_tmp * sin_latl;
                pmcx = x_iaumars - cam_C[0];
                pmcy = y_iaumars - cam_C[1];
                pmcz = z_iaumars - cam_C[2];
                
                apmc =  pmcx*cam_A[0] + pmcy*cam_A[1] + pmcz*cam_A[2];
                ppvx = (pmcx*cam_H[0] + pmcy*cam_H[1] + pmcz*cam_H[2])/apmc;
                ppvy = (pmcx*cam_V[0] + pmcy*cam_V[1] + pmcz*cam_V[2])/apmc;
                xi = (int32_T) floor(ppvx);
                yi = (int32_T) floor(ppvy);
                //xi = xi<0?0:xi;
                //xi = xi>S_imm1?:xi;
                if(xi<0)
                    xi=0;
                else if(xi>S_imm1)
                    xi = S_imm1;

                if(yi<0)
                    yi=0;
                else if(yi>L_imm1)
                    yi = L_imm1;

                ++bin_count_im[xi][yi];
            } else if(msldemc_imFOVmask[c][l]==1) {
                /*  */
                (*count_napmc)++;
            }
        }
    }
    
    for(xi=0;xi<S_im;xi++){
        for(yi=0;yi<L_im;yi++){
            if(bin_count_im[xi][yi]>0){
                bin_im_c[xi][yi] = (int32_T*) malloc(bin_count_im[xi][yi]*sizeof(int32_T));
                bin_im_l[xi][yi] = (int32_T*) malloc(bin_count_im[xi][yi]*sizeof(int32_T));
                bin_imx[xi][yi] = (double*) malloc(bin_count_im[xi][yi]*sizeof(double));
                bin_imy[xi][yi] = (double*) malloc(bin_count_im[xi][yi]*sizeof(double));
                bin_rad[xi][yi] = (double*) malloc(bin_count_im[xi][yi]*sizeof(double));
            } else {
                bin_im_c[xi][yi] = NULL;
                bin_im_l[xi][yi] = NULL;
                bin_imx[xi][yi] = NULL;
                bin_imy[xi][yi] = NULL;
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
    
    for(xi=0;xi<S_im;xi++){
        for(yi=0;yi<L_im;yi++){
            bin_count_im[xi][yi] = 0;
        }
    }
    
    
    for(l=0;l<msldemc_lines;l++){
        cos_latl   = cos(msldemc_latitude[l]);
        sin_latl   = sin(msldemc_latitude[l]);
        for(c=0;c<msldemc_samples;c++){
            if(msldemc_imFOVmask[c][l]>1){
                radius_tmp = msldemc_img[c][l];
                x_iaumars  = radius_tmp * cos_latl * cos_lon[c];
                y_iaumars  = radius_tmp * cos_latl * sin_lon[c];
                z_iaumars  = radius_tmp * sin_latl;
                pmcx = x_iaumars - cam_C[0];
                pmcy = y_iaumars - cam_C[1];
                pmcz = z_iaumars - cam_C[2];
                
                apmc   = pmcx*cam_A[0] + pmcy*cam_A[1] + pmcz*cam_A[2];
                
                ppvx = (pmcx*cam_H[0] + pmcy*cam_H[1] + pmcz*cam_H[2])/apmc;
                ppvy = (pmcx*cam_V[0] + pmcy*cam_V[1] + pmcz*cam_V[2])/apmc;
                xi = (int32_T) floor(ppvx);
                yi = (int32_T) floor(ppvy);
                if(xi<0)
                    xi=0;
                else if(xi>S_imm1)
                    xi = S_imm1;

                if(yi<0)
                    yi=0;
                else if(yi>L_imm1)
                    yi = L_imm1;

                bin_im_c[xi][yi][bin_count_im[xi][yi]] = c;
                bin_im_l[xi][yi][bin_count_im[xi][yi]] = l;
                bin_imx[xi][yi][bin_count_im[xi][yi]]  = ppvx;
                bin_imy[xi][yi][bin_count_im[xi][yi]]  = ppvy;
                bin_rad[xi][yi][bin_count_im[xi][yi]]  = radius_tmp;
                ++bin_count_im[xi][yi];

            } else if(msldemc_imFOVmask[c][l]==1){
                radius_tmp = msldemc_img[c][l];
                (*c_napmc)[*count_napmc] = c;
                (*l_napmc)[*count_napmc] = l;
                (*rad_napmc)[*count_napmc] = radius_tmp;
                (*count_napmc)++;
            }
        }
    }
    
    free(cos_lon);
    free(sin_lon);
}
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double **msldemc_img;
    mwSize msldemc_samples,msldemc_lines;
    double *msldemc_latitude;
    double *msldemc_longitude;
    CAHV_MODEL cahv_mdl;
    int8_T **msldemc_imFOVmask;
    mwSize S_im,L_im;
    
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
    
    /* INPUT 0 msldemc_img */
    msldemc_img     = set_mxDoubleMatrix(prhs[0]);
    msldemc_lines   = mxGetM(prhs[0]);
    msldemc_samples = mxGetN(prhs[0]);
    
    
    /* INPUT 1/2 msldemc latitude & longitude */
    msldemc_latitude  = mxGetDoubles(prhs[1]);
    msldemc_longitude = mxGetDoubles(prhs[1]);
    
    
    /* INPUT 3 msldemc imFOV */
    msldemc_imFOVmask = set_mxInt8Matrix(prhs[2]);
    
    /* INPUT 4/5 image S_im, L_im */
    S_im = (mwSize) mxGetScalar(prhs[3]);
    L_im = (mwSize) mxGetScalar(prhs[4]);
    
    /* INPUT 6 CAHV model */
    cahv_mdl = mxGet_CAHV_MODEL(prhs[5]);
    
    
    /* OUTPUT 0 msldemc imFOV */
    plhs[0] = mxCreateNumericMatrix(msldemc_lines,msldemc_samples,mxINT8_CLASS,mxREAL);
    msldemt_imUFOVmask = set_mxInt8Matrix(plhs[0]);
    
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
    createInt32Matrix(&bin_count_im,&bin_count_im_base,(size_t) S_im, (size_t) L_im);
    createInt32PMatrix(&bin_im_c, &bin_im_c_base, (size_t) S_im, (size_t) L_im);
    createInt32PMatrix(&bin_im_l, &bin_im_l_base, (size_t) S_im, (size_t) L_im);
    createDoublePMatrix(&bin_imx, &bin_imx_base, (size_t) S_im, (size_t) L_im);
    createDoublePMatrix(&bin_imy, &bin_imy_base, (size_t) S_im, (size_t) L_im);
    createDoublePMatrix(&bin_rad, &bin_rad_base, (size_t) S_im, (size_t) L_im);
    
    bin_msldemt_xyz_wAHVint_iaumars_ctr((int32_T) S_im, (int32_T) L_im, cahv_mdl,
            msldemc_img, msldemc_samples, msldemc_lines, 
            msldemc_latitude, msldemc_longitude, msldemc_imFOVmask,
            bin_count_im, bin_im_c, bin_im_l, bin_imx, bin_imy, bin_rad,
            &count_napmc,&c_napmc,&l_napmc,&rad_napmc);
    
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */    
    mask_obstructed_pts_in_msldemt_using_msldemc_iaumars(
            msldemc_img,
        (int32_T) msldemc_samples, (int32_T) msldemc_lines,
        msldemc_latitude, msldemc_longitude, msldemc_imFOVmask,
        (int32_T) msldemc_samples, (int32_T) msldemc_lines,
        msldemc_latitude, msldemc_longitude, msldemt_imUFOVmask,
        bin_count_im, bin_im_c, bin_im_l,
        bin_imx, bin_imy, bin_rad,
        count_napmc, c_napmc, l_napmc, rad_napmc,
        S_im, L_im, cahv_mdl);
        
    
    /* Freeing memories */
    for(si=0;si<S_im;si++){
        for(li=0;li<L_im;li++){
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
    
    
    // mxFree(bin_count_im);
    free(bin_count_im);
    free(bin_count_im_base);
    
    if(count_napmc>0){
        free(c_napmc);
        free(l_napmc);
        free(rad_napmc);
    }
    
    /* free memories */
    mxFree(msldemc_img);
    mxFree(msldemc_imFOVmask);
    mxFree(msldemt_imUFOVmask);
    
}