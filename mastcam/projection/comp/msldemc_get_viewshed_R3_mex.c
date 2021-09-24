/* =====================================================================
 * viewshed_get_msldemtUFOVmask_ctr_R3_mex.c
 * Get viewshed within the FOV using XDraw algorithm
 * 
 * INPUTS:
 * 0 msldemc_img        float image
 * 1 msldemc_imFOVmask     int8_t [L_demc x S_demc]
 * 2 cahv_mdl              CAHV_MODEL
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
#include "msldem_util.h"
#include <time.h>


void get_viewshed_R3(double **msldemc_img, int32_t msldemc_samples,
        int32_t msldemc_lines, CAHV_MODEL cahv_mdl, 
        int8_t **msldemc_imFOVmask, int8_t **msldemc_imUFOVmask)
{
    int32_t s,l,c;
    int32_t msldemc_samplesm1, msldemc_linesm1;
    double *cam_C;
    double d[3];
    double plane_const;
    double s_dbl,l_dbl;
    double elevcl_dbl;
    double incr;
    int32_t ss,ll;
    double ll_dbl,s_incpt_ll,s_ref_1_dbl,s_ref_2_dbl,elev_incpt,elev_incpt_los;
    int32_t s_ref_1,s_ref_2;
    double ss_dbl,l_incpt_ss,l_ref_1_dbl,l_ref_2_dbl;
    int32_t l_ref_1,l_ref_2;
    int32_t ll_strt, ss_strt;
    int flg;
    
    cam_C = cahv_mdl.C;
    msldemc_samplesm1 = msldemc_samples-1;
    msldemc_linesm1   = msldemc_lines-1;
    
    for(l=0;l<6000;l++){
        // printf("l=%d\n",l);
        for(s=0;s<msldemc_samples;s++){
            if(msldemc_imFOVmask[s][l]>0){
                elevcl_dbl = (double) msldemc_img[s][l];
                /* Get the angle between C  */
                l_dbl = (double) l; s_dbl = (double) s;
                
                /* line parameter */
                /* d[0]*y -d[1]*x = plane_const */
                /* y-axis: sample, x-axis: line */
                d[0] = l_dbl - cam_C[0];
                d[1] = s_dbl - cam_C[1];
                d[2] = elevcl_dbl - cam_C[2];
                plane_const = d[0]*cam_C[1]-d[1]*cam_C[0];
                
                /* parse intercept with x=ll */
                flg = 1;
                if(d[0]>0){
                    ll_strt = (int32_t) ceil(cam_C[0]);
                    if(ll_strt<0)
                        ll_strt = 0;
                    if(ll_strt>msldemc_linesm1)
                        ll_strt = msldemc_linesm1;
                    for(ll=ll_strt;ll<l;ll++){
                        ll_dbl = (double) ll;
                        s_incpt_ll = (d[1]*ll_dbl+plane_const)/d[0];/* intercept at x=ll */
                        if(s_incpt_ll<0 || s_incpt_ll>msldemc_samplesm1){
                            continue;
                        }
                        s_ref_1_dbl = floor(s_incpt_ll);
                        s_ref_2_dbl = ceil(s_incpt_ll);
                        s_ref_1 = (int32_t) s_ref_1_dbl;
                        s_ref_2 = (int32_t) s_ref_2_dbl;
                        elev_incpt = msldemc_img[s_ref_1][ll] * (s_ref_2_dbl-s_incpt_ll) + msldemc_img[s_ref_2][ll] * (s_incpt_ll-s_ref_1_dbl);
                        elev_incpt_los = d[2]*(ll_dbl-cam_C[0])/d[0] + cam_C[2];
                        if(elev_incpt>elev_incpt_los){
                            msldemc_imUFOVmask[s][l] = 0;
                            break;
                            flg=0;
                        }
                    }
                } else if(d[0]<0){
                    
                }
                
                if(flg){
                    if(d[1]>0){
                        
                    } else if(d[1]<0){
                        ss_strt = (int32_t) ceil(cam_C[1]);
                        if(ss_strt<0)
                            ss_strt = 0;
                        if(ss_strt>msldemc_samplesm1)
                            ss_strt = msldemc_samplesm1;
                        for(ss=ss_strt;ss>s;ss--){
                            ss_dbl = (double) ss;
                            l_incpt_ss = (d[0]*ss_dbl-plane_const)/d[1];/* intercept at y=ss */
                            if(l_incpt_ss<0 || l_incpt_ss>msldemc_linesm1){
                                continue;
                            }
                            l_ref_1_dbl = floor(l_incpt_ss);
                            l_ref_2_dbl = ceil(l_incpt_ss);
                            l_ref_1 = (int32_t) l_ref_1_dbl;
                            l_ref_2 = (int32_t) l_ref_2_dbl;
                            elev_incpt = msldemc_img[ss][l_ref_1] * (l_ref_2_dbl-l_incpt_ss) + msldemc_img[ss][l_ref_2] * (l_incpt_ss-l_ref_1_dbl);
                            elev_incpt_los = d[2]*(ss_dbl-cam_C[1])/d[1] + cam_C[2];
                            if(elev_incpt>elev_incpt_los){
                                msldemc_imUFOVmask[s][l] = 0;
                                flg = 0;
                                break;
                            }
                        }
                    }
                }
                
//                 while(flg){
//                     ll = (int32_t) ll_dbl;
//                     if(ll>-1 && ll<msldemc_lines){
//                         s_incpt_ll = (d[1]*ll_dbl+plane_const)/d[0];/* intercept at x=ll */
//                         if(s_incpt_ll>0 && s_incpt_ll<msldemc_samplesm1){
//                             s_ref_1_dbl = floor(s_incpt_ll);
//                             s_ref_2_dbl = ceil(s_incpt_ll);
//                             s_ref_1 = (int32_t) s_ref_1_dbl;
//                             s_ref_2 = (int32_t) s_ref_2_dbl;
//                             elev_incpt = msldemc_img[s_ref_1][ll] * (s_ref_2_dbl-s_incpt_ll) + msldemc_img[s_ref_2][ll] * (s_incpt_ll-s_ref_1_dbl);
//                             elev_incpt_los = d[2]*(ll_dbl-cam_C[0])/d[0] + cam_C[2];
//                             if(elev_incpt>elev_incpt_los){
//                                 msldemc_imUFOVmask[s][l] = 0;
//                                 flg = 0;
//                             } else {
//                                 ll_dbl += incr;
//                             }
//                         } else {
//                             ll_dbl += incr;
//                         }
//                     } else {
//                         ll_dbl += incr;
//                     }
//                     if(ll==l){
//                         flg = 0;
//                     }
//                 }
//                 
//                 /* parse intercept with y = ss */
//                 if(msldemc_imUFOVmask[s][l] > 0){
//                     flg = 1;
//                 }
//                 
//                 if(d[1]>0){
//                     ss_dbl = ceil(cam_C[1]);
//                     incr = 1;
//                 } else if(d[1]<0){
//                     ss_dbl = floor(cam_C[1]);
//                     incr = -1;
//                 }
                
//                 while(flg){
//                     ss = (int32_t) ss_dbl;
//                     if(ss>-1 && ss<msldemc_samples){
//                         l_incpt_ss = (d[0]*ss_dbl-plane_const)/d[1];/* intercept at x=ll */
//                         if(l_incpt_ss>0 && l_incpt_ss<msldemc_linesm1){
//                             l_ref_1_dbl = floor(l_incpt_ss);
//                             l_ref_2_dbl = ceil(l_incpt_ss);
//                             l_ref_1 = (int32_t) l_ref_1_dbl;
//                             l_ref_2 = (int32_t) l_ref_2_dbl;
//                             elev_incpt = msldemc_img[ss][l_ref_1] * (l_ref_2_dbl-l_incpt_ss) + msldemc_img[ss][l_ref_2] * (l_incpt_ss-l_ref_1_dbl);
//                             elev_incpt_los = d[2]*(ss_dbl-cam_C[1])/d[1] + cam_C[2];
//                             if(elev_incpt>elev_incpt_los){
//                                 msldemc_imUFOVmask[s][l] = 0;
//                                 flg = 0;
//                             } else {
//                                 ss_dbl += incr;
//                             }
//                         } else {
//                             ss_dbl += incr;
//                         }
//                     } else {
//                         ss_dbl += incr;
//                     }
//                     
//                     if(ss==s){
//                         flg=0;
//                     }
//                 }
            }


        }
    }
    
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double **msldemc_img;
    mwSize msldemc_samples,msldemc_lines;
    CAHV_MODEL cahv_mdl;
    int8_t **msldemc_imFOVmask;
    int8_t **msldemt_imUFOVmask;
    
    clock_t strt_time, end_time;
    double cpu_time_used;
    mwSize si,li;

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
    msldemc_samples = mxGetN(prhs[0]);
    msldemc_lines   = mxGetM(prhs[0]);
    
    /* INPUT 3 msldemc imFOV */
    msldemc_imFOVmask = set_mxInt8Matrix(prhs[1]);
    
    /* INPUT 4 CAHV model */
    cahv_mdl = mxGet_CAHV_MODEL(prhs[2]);
    
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
    printf("a\n");
    /* Main computation routine */
    get_viewshed_R3(msldemc_img, (int32_t) msldemc_samples, (int32_t) msldemc_lines, cahv_mdl, 
        msldemc_imFOVmask, msldemt_imUFOVmask);
    
    
    mxFree(msldemc_img);
    mxFree(msldemc_imFOVmask);
    mxFree(msldemt_imUFOVmask);

    
}