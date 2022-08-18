/* =====================================================================
 * iaumars_get_msldemtUFOVmask_ctr_L1PBK_LL0_M3_mex.c
 * L1  : msldemc is an input 2d array
 * PBK : Prior Binning into bins with the auxiliary size defined by two parameters K_L and K_S
 * LL0 : Linked List with least basic information (c,l,radius), this is most memory efficient
 * M3  : 3x3 matrix inversion object image coordinate.
 * 
 * INPUTS:
 * 0 msldemc_img           Double array [L_demc x S_demc]
 * 1 msldemc_latitude      Double array [L_demc]
 * 2 msldemc_longitude     Double array [S_demc]
 * 3 msldemc_imFOVmask     int8_T [L_demc x S_demc]
 * 4 S_im                  int
 * 5 L_im                  int
 * 6 cahv_mdl              CAHV_MODEL
 * 7 K_L                   double reciprocal of the length of the bin in the image line direction.
 * 8 K_S                   double reciprocal of the length of the bin in the image sample direction.
 * 9 dyu                   int8_t flag for dynamic linked list update.
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
#include "lib_proj_mastcamMSLDEM_IAUMars_L1PBK_LL0_M3.h"

#include <time.h>



void bin_msldemt_ctr_iaumars_L1PBK_LL0(double S_im, double L_im, CAHV_MODEL cahv_mdl,
        double **msldemc_img,
        int32_T msldemc_samples, int32_T msldemc_lines,
        double *msldemc_latitude, double *msldemc_longitude, 
        int8_T **msldemc_imFOVmask, 
        double K_L, double K_S,
        struct MSLDEMmask_LinkedList ***ll_papmc,
        struct MSLDEMmask_LinkedList **ll_napmc)
{
    int32_T xi,yi,c,l;
    int32_t binL,binS;
    int32_t binLm1,binSm1;
    double pmcx,pmcy,pmcz,apmc,ppvx,ppvy;
    double *cam_C, *cam_A, *cam_H, *cam_V;
    
    // double dem_cl;
    
    double *cos_lon, *sin_lon;
    double radius_tmp;
    double cos_latl, sin_latl;
    double x_iaumars, y_iaumars, z_iaumars;
    struct MSLDEMmask_LinkedList *ll_papmc_next;
    struct MSLDEMmask_LinkedList *ll_napmc_next;
    
    ll_napmc_next = (*ll_napmc);
    
    // S_imm1 = S_im - 1;
    // L_imm1 = L_im - 1;
    
    binL = (int32_t) (K_L * L_im); binS = (int32_t) (K_S * S_im);
    binLm1 = binL-1; binSm1 = binS-1;
    
    cam_C = cahv_mdl.C; cam_A = cahv_mdl.A; cam_H = cahv_mdl.H; cam_V = cahv_mdl.V;
    
    
    cos_lon = (double*) malloc(sizeof(double) * (size_t) msldemc_samples);
    sin_lon = (double*) malloc(sizeof(double) * (size_t) msldemc_samples);
    for(c=0;c<msldemc_samples;c++){
        cos_lon[c] = cos(msldemc_longitude[c]);
        sin_lon[c] = sin(msldemc_longitude[c]);
    }

    
    
    // printf("%d\n",msldemc_samples);
    for(l=0;l<msldemc_lines;l++){
        /* read a line */ 
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
                
                if(ll_papmc[xi][yi]!=NULL){
                    ll_papmc_next = malloc(sizeof(struct MSLDEMmask_LinkedList));
                    ll_papmc_next->c = c;
                    ll_papmc_next->l = l;
                    ll_papmc_next->radius = radius_tmp;
                    ll_papmc_next->next = ll_papmc[xi][yi];
                    ll_papmc[xi][yi]->prev = ll_papmc_next;
                    ll_papmc[xi][yi] = ll_papmc_next;
                    ll_papmc_next->prev = NULL;
                } else {
                    ll_papmc_next = malloc(sizeof(struct MSLDEMmask_LinkedList));
                    ll_papmc_next->c = c;
                    ll_papmc_next->l = l;
                    ll_papmc_next->radius = radius_tmp;
                    ll_papmc_next->next = NULL;
                    ll_papmc[xi][yi] = ll_papmc_next;
                    ll_papmc_next->prev = NULL;
                }

            } else if(msldemc_imFOVmask[c][l]==1) {
                /*  */
                if(ll_napmc_next!=NULL){
                    ll_napmc_next->next = malloc(sizeof(struct MSLDEMmask_LinkedList));
                    ll_napmc_next->next->prev = ll_napmc_next;
                    ll_napmc_next = ll_napmc_next->next;
                    ll_napmc_next->c = c;
                    ll_napmc_next->l = l;
                    ll_napmc_next->radius = radius_tmp;
                    ll_napmc_next->next = NULL;
                } else {
                    (*ll_napmc) = malloc(sizeof(struct MSLDEMmask_LinkedList));
                    ll_napmc_next = (*ll_napmc);
                    ll_napmc_next->c = c;
                    ll_napmc_next->l = l;
                    ll_napmc_next->radius = radius_tmp;
                    ll_napmc_next->next = NULL;
                    ll_napmc_next->prev = NULL;
                }
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
    mwSize msldemc_samples, msldemc_lines;
    double *msldemc_latitude;
    double *msldemc_longitude;
    double mslrad_offset;
    CAHV_MODEL cahv_mdl;
    int8_T **msldemc_imFOVmask;
    double S_im,L_im;
    int8_t dyu;
    
    int8_T **msldemt_imUFOVmask;
    
    mwIndex si,li;
    
    struct MSLDEMmask_LinkedList ***ll_papmc_bin, **ll_papmc_bin_base;
    struct MSLDEMmask_LinkedList *ll_napmc;
    struct MSLDEMmask_LinkedList *ll_tmp;
    
    
    double K_L,K_S;
    mwSize binL,binS;
    
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
    msldemc_img = set_mxDoubleMatrix(prhs[0]);
    msldemc_lines = mxGetM(prhs[0]);
    msldemc_samples = mxGetN(prhs[0]);
    
    /* INPUT 1/2 msldem northing easting */
    msldemc_latitude  = mxGetDoubles(prhs[1]);
    msldemc_longitude = mxGetDoubles(prhs[2]);
    
    
    /* INPUT 3 msldemc imFOV */
    msldemc_imFOVmask = set_mxInt8Matrix(prhs[3]);
    
    /* INPUT 4/5 image S_im, L_im */
    S_im = mxGetScalar(prhs[4]);
    L_im = mxGetScalar(prhs[5]);
    
    /* INPUT 6 CAHV model */
    cahv_mdl = mxGet_CAHV_MODEL(prhs[6]);
    
    /* INPUT 7/8 bin size parameters*/
    K_L = mxGetScalar(prhs[7]);
    K_S = mxGetScalar(prhs[8]);
    
    /* INPUT 9 Dynamic Linked List Update FLAG */
    dyu = (int8_t) mxGetScalar(prhs[9]);
    
    
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
    ll_napmc = NULL;
    ll_tmp   = NULL;
    strt_time = clock();
    binL = (mwSize) (K_L * L_im); binS = (mwSize) (K_S * S_im);
    createMSLDEMmask_LLMatrix(&ll_papmc_bin, &ll_papmc_bin_base, (size_t) binS, (size_t) binL);
    for(si=0;si<binS;si++){
        for(li=0;li<binL;li++){
            ll_papmc_bin[si][li] = NULL;
        }
    }
    bin_msldemt_ctr_iaumars_L1PBK_LL0(S_im, L_im, cahv_mdl,
            msldemc_img, msldemc_samples, msldemc_lines, 
            msldemc_latitude, msldemc_longitude, msldemc_imFOVmask,
            K_L,K_S,
            ll_papmc_bin, &ll_napmc);
    end_time = clock();
    cpu_time_used = ((double) (end_time - strt_time)) / CLOCKS_PER_SEC;
    printf("Preprocessing took %f [sec].\n",cpu_time_used);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    strt_time = clock();
    if(dyu==1){
        mask_obstructed_pts_in_msldemt_using_msldemc_iaumars_L1PBK_LL0DYU_M3(
            msldemc_img,
            (int32_T) msldemc_samples, (int32_T) msldemc_lines,
            msldemc_latitude, msldemc_longitude, msldemc_imFOVmask,
            (int32_T) msldemc_samples, (int32_T) msldemc_lines,
            msldemc_latitude, msldemc_longitude, msldemt_imUFOVmask,
            K_L,K_S,
            ll_papmc_bin, ll_napmc,
            S_im, L_im, cahv_mdl);
    } else if(dyu==0){
//         mask_obstructed_pts_in_msldemt_using_msldemc_iaumars_L1PBK_LL0_M3(
//             msldemc_img,
//             (int32_T) msldemc_samples, (int32_T) msldemc_lines,
//             msldemc_latitude, msldemc_longitude, msldemc_imFOVmask,
//             (int32_T) msldemc_samples, (int32_T) msldemc_lines,
//             msldemc_latitude, msldemc_longitude, msldemt_imUFOVmask,
//             K_L,K_S,
//             ll_papmc_bin, ll_napmc,
//             S_im, L_im, cahv_mdl);
    } else {
        printf("Undefined DYU_FLAG=%d\n",dyu);
    }
        
    end_time = clock();
    cpu_time_used = ((double) (end_time - strt_time)) / CLOCKS_PER_SEC;
    printf("Main operation took %f [sec].\n",cpu_time_used);
    
    /* Freeing memories */
    for(si=0;si<binS;si++){
        for(li=0;li<binL;li++){
            while(ll_papmc_bin[si][li] != NULL){
                ll_tmp   = ll_papmc_bin[si][li];
                ll_papmc_bin[si][li] = ll_papmc_bin[si][li]->next;
                free(ll_tmp);
            }
        }
    }
    free(ll_papmc_bin);
    free(ll_papmc_bin_base);
    
    while(ll_napmc != NULL){
        ll_tmp   = ll_napmc;
        ll_napmc = ll_napmc->next;
        free(ll_tmp);
    }
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldemc_imFOVmask);
    mxFree(msldemt_imUFOVmask);
    
}