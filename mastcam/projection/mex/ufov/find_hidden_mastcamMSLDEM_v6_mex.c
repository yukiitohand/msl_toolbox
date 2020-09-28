/* =====================================================================
 * find_hidden_mastcamMSLDEM_v6_mex.c
 * Evaluate any pixels in MSLDEM image whether or not they exist in the 
 * 
 * INPUTS:
 * 0 msldemc_img           Float array [L_demc x S_demc]
 * 1 msldemc_northing      Double array [L_demc]
 * 2 msldemc_easting       Double array [S_demc]
 * 3 msldemc_imx           Doubles [L_demc x S_demc]
 * 4 msldemc_imy           Doubles [L_demc x S_demc]
 * 5 msldemc_imFOVmask     int8_T [L_demc x S_demc]
 * 6 S_im                  int
 * 7 L_im                  int
 *
 * The origin of msldemc_img, msldemc_northing, and msldemc_easting is 
 * cam_C.
 * 
 * 
 * OUTPUTS:
 * 0  msldemc_inImage     [L_demc x S_demc]  Boolean
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
void find_hidden(int32_T msldemc_samples, int32_T msldemc_lines,
        double **msldemc_img,
        double *msldemc_northing, double *msldemc_easting,
        double **msldemc_imx, double **msldemc_imy,
        int8_T **msldemc_imFOVmask, 
        int32_T S_im, int32_T L_im,
        int8_T **msldemc_inImage)
{
    int32_T c,l,cc,ll;
    int32_T cv1,cv2,cv3,cv4,lv1,lv2,lv3,lv4; /* v indicates vertex */
    int32_T L_demcm1,S_demcm1;
    int16_T ti; /* triangle index */
    double tppvx,tppvy,tppvz,tppvgx,tppvgy,tppvgz;
    double ppv1x,ppv1y,ppv2x,ppv2y,ppv3x,ppv3y; /* Plane Position Vectors */
    double ppv1gx,ppv1gy,ppv1gz,ppv2gx,ppv2gy,ppv2gz,ppv3gx,ppv3gy,ppv3gz,ppv4gx,ppv4gy,ppv4gz;
    double pdv1x,pdv1y,pdv2x,pdv2y;
    double pdv1gx,pdv1gy,pdv1gz,pdv2gx,pdv2gy,pdv2gz;
    double detM;
    double M[2][2];
    double Minv[2][2];
    double Minvp[2][3];
    double x_min,y_min,x_max,y_max;
    double pipvx,pipvy;
    double pipvgx,pipvgy,pipvgz;
    double pipvgppv1x,pipvgppv1y,pipvgppv1z;
    double pprm_sd,pprm_td,pprm_1std; /* plane parameter for projected image plane */
    double pprm_s,pprm_t,pprm_1st;
    bool isinFOVd,isinFOV;
    
    double pnx,pny,pnz; /* Plane Normal vectors */
    // double pc; /* Plane Constant */
    double lprm_nume;
    double lprm; /* line parameters */
    double S_imm05,L_imm05;
    
    int32_T **bin_count_im, *bin_count_im_base;
    int32_T **bin_count_im_idx_ctr, *bin_count_im_idx_ctr_base;
    int32_T ***bin_im_c, **bin_im_c_base;
    int32_T ***bin_im_l, **bin_im_l_base;
    int32_T xi,yi;
    int32_T x_min_int,x_max_int,y_min_int,y_max_int;
    int32_T n;
    int32_T S_imm1,L_imm1;
    
    /* Dynamic memory allocations */
    
    bin_count_im = (int32_T**) malloc(sizeof(int32_T*) * (size_t) S_im);
    bin_count_im_base = (int32_T*) malloc(sizeof(int32_T) * (size_t) S_im * (size_t) L_im);
    bin_count_im[0] = &bin_count_im_base[0];
    for(xi=1;xi<S_im;xi++){
        bin_count_im[xi] = bin_count_im[xi-1] + L_im;
    }
    
    bin_count_im_idx_ctr = (int32_T**) malloc(sizeof(int32_T*) * (size_t) S_im);
    bin_count_im_idx_ctr_base = (int32_T*) malloc(sizeof(int32_T) * (size_t) S_im * (size_t) L_im);
    bin_count_im_idx_ctr[0] = &bin_count_im_idx_ctr_base[0];
    for(xi=1;xi<S_im;xi++){
        bin_count_im_idx_ctr[xi] = bin_count_im_idx_ctr[xi-1] + L_im;
    }
    
    bin_im_c = (int32_T***) malloc(sizeof(int32_T**) * (size_t) S_im);
    bin_im_c_base = (int32_T**) malloc(sizeof(int32_T*) * (size_t) S_im * (size_t) L_im);
    bin_im_c[0] = &bin_im_c_base[0];
    for(xi=1;xi<S_im;xi++){
        bin_im_c[xi] = bin_im_c[xi-1] + L_im;
    }
    
    bin_im_l = (int32_T***) malloc(sizeof(int32_T**) * (size_t) S_im);
    bin_im_l_base = (int32_T**) malloc(sizeof(int32_T*) * (size_t) S_im * (size_t) L_im);
    bin_im_l[0] = &bin_im_l_base[0];
    for(xi=1;xi<S_im;xi++){
        bin_im_l[xi] = bin_im_l[xi-1] + L_im;
    }
            
    
    
    // printf("%d\n",msldemc_samples);
    L_demcm1 = msldemc_lines-1;
    S_demcm1 = msldemc_samples-1;
    S_imm1 = S_im - 1;
    L_imm1 = L_im - 1;
    S_imm05 = (double) S_im - 0.5;
    L_imm05 = (double) L_im - 0.5;
    printf("%f,%f\n",S_imm05,L_imm05);
    // printf("%d\n",msldemc_samples);
    
    /* create an bin image counting the number of demc pixels that falls
     * within the */
    
    /* Initialize */
    for(xi=0;xi<S_im;xi++){
        for(yi=0;yi<L_im;yi++){
            bin_count_im[xi][yi] = 0;
            bin_count_im_idx_ctr[xi][yi] = 0;
        }
    }
    // printf("%d\n",msldemc_samples);
    for(ll=0;ll<msldemc_lines;ll++){
        for(cc=0;cc<msldemc_samples;cc++){
            if(msldemc_imFOVmask[cc][ll]){
                xi = (int32_T) floor(msldemc_imx[cc][ll]);
                yi = (int32_T) floor(msldemc_imy[cc][ll]);
                if(xi<0)
                    xi=0;
                else if(xi>S_imm1)
                    xi = S_imm1;
                
                if(yi<0)
                    yi=0;
                else if(yi>L_imm1)
                    yi = L_imm1;
                    
                ++bin_count_im[xi][yi];
            }
        }
    }
    // printf("%d\n",msldemc_samples);
    for(xi=0;xi<S_im;xi++){
        for(yi=0;yi<L_im;yi++){
            if(bin_count_im[xi][yi]>0){
                bin_im_c[xi][yi] = (int32_T*) malloc(bin_count_im[xi][yi]*sizeof(int32_T));
                bin_im_l[xi][yi] = (int32_T*) malloc(bin_count_im[xi][yi]*sizeof(int32_T));
            } else {
                bin_im_c[xi][yi] = NULL;
                bin_im_l[xi][yi] = NULL;
            }
        }
    }
    // printf("%d\n",msldemc_samples);
    for(ll=0;ll<msldemc_lines;ll++){
        for(cc=0;cc<msldemc_samples;cc++){
            if(msldemc_imFOVmask[cc][ll]){
                xi = (int32_T) floor(msldemc_imx[cc][ll]);
                yi = (int32_T) floor(msldemc_imy[cc][ll]);
                if(xi<0)
                    xi=0;
                else if(xi>S_imm1)
                    xi = S_imm1;
                
                if(yi<0)
                    yi=0;
                else if(yi>L_imm1)
                    yi = L_imm1;
                
                bin_im_c[xi][yi][bin_count_im_idx_ctr[xi][yi]] = cc;
                bin_im_l[xi][yi][bin_count_im_idx_ctr[xi][yi]] = ll;
                ++bin_count_im_idx_ctr[xi][yi];
            }
        }
    }
    
    
    free(bin_count_im_idx_ctr);
    free(bin_count_im_idx_ctr_base);
    

    //printf("%d\n",msldemc_samples);
    
    //L_demcm1_dbl = (double) L_demcm1_dbl;
    //S_demcm1_dbl = (double) S_demcm1_dbl;
    // printf("%d\n",msldemc_samples);
    
    for(c=0;c<msldemc_samples;c++){
        for(l=0;l<msldemc_lines;l++){
            if((msldemc_imx[c][l]<-0.5) || (msldemc_imx[c][l]>S_imm05) || (msldemc_imy[c][l]<-0.5) || (msldemc_imy[c][l]>L_imm05)){
                msldemc_inImage[c][l] = 0;
            }
        }
    }
    
                
    for(l=0;l<L_demcm1;l++){
        //printf("l=%d\n",l);
        // decide the first and last indexes to be assessed.
        for(c=0;c<S_demcm1;c++){
            // process if 
            //printf("c=%d,mask_lc = %d,mask_lp1c = %d\n",c,msldemc_imFOVmask[c][l],msldemc_imFOVmask[c][l+1]);
            // printf("c=%d/%d\n",c,S_demcm1);
            if((msldemc_imFOVmask[c][l]==2) || (msldemc_imFOVmask[c][l+1]==2)){
                //printf("c=%d\n",c);
                for(ti=0;ti<2;ti++){
                    if(ti==0){
//                         ppv1x = msldemc_imx[c][l];
//                         ppv1y = msldemc_imy[c][l];
//                         ppv2x = msldemc_imx[c+1][l];
//                         ppv2y = msldemc_imy[c+1][l];
//                         ppv3x = msldemc_imx[c][l+1];
//                         ppv3y = msldemc_imy[c][l+1];
//                         ppv1gx = msldemc_northing[l];
//                         ppv1gy = msldemc_easting[c];
//                         ppv1gz = msldemc_img[c][l];
//                         ppv2gx = msldemc_northing[l];
//                         ppv2gy = msldemc_easting[c+1];
//                         ppv2gz = msldemc_img[c+1][l];
//                         ppv3gx = msldemc_northing[l+1];
//                         ppv3gy = msldemc_easting[c];
//                         ppv3gz = msldemc_img[c][l+1];
//                         ppv4gx = msldemc_northing[l+1];
//                         ppv4gy = msldemc_easting[c+1];
//                         ppv4gz = msldemc_img[c+1][l+1];
                        // ppvx[0] = c; ppvx[1] = c+1; ppvx[2] = c; ppvx[3] = c+1;
                        // ppvy[0] = l; ppvy[1] = l; ppvy[2] = l+1; ppvy[3] = l+1;
                        isinFOVd = ((msldemc_imFOVmask[c][l]==2) && (msldemc_imFOVmask[c+1][l]==2) && (msldemc_imFOVmask[c][l+1]==2));
                        isinFOV = ((msldemc_imFOVmask[c][l]>0) && (msldemc_imFOVmask[c+1][l]>0) && (msldemc_imFOVmask[c][l+1]>0));
                        cv1 = c;   lv1 = l;
                        cv2 = c+1; lv2 = l;
                        cv3 = c;   lv3 = l+1;
                        cv4 = c+1; lv4 = l+1;
                    }
                    else{
                        cv1 = c+1; lv1 = l;
                        cv2 = c+1; lv2 = l+1;
                        cv3 = c;   lv3 = l+1;
                        cv4 = c;   lv4 = l;
//                         ppv1x = msldemc_imx[c+1][l];
//                         ppv1y = msldemc_imy[c+1][l];
//                         ppv2x = msldemc_imx[c+1][l+1];
//                         ppv2y = msldemc_imy[c+1][l+1];
//                         ppv3x = msldemc_imx[c][l+1];
//                         ppv3y = msldemc_imy[c][l+1];
//                         ppv1gx = msldemc_northing[l];
//                         ppv1gy = msldemc_easting[c+1];
//                         ppv1gz = msldemc_img[c+1][l];
//                         ppv2gx = msldemc_northing[l+1];
//                         ppv2gy = msldemc_easting[c+1];
//                         ppv2gz = msldemc_img[c+1][l+1];
//                         ppv3gx = msldemc_northing[l+1];
//                         ppv3gy = msldemc_easting[c];
//                         ppv3gz = msldemc_img[c][l+1];
//                         ppv4gx = msldemc_northing[l];
//                         ppv4gy = msldemc_easting[c];
//                         ppv4gz = msldemc_img[c][l];
                        //ppvx[0] = c+1; ppvx[1] = c+1; ppvx[2] = c; ppvx[3] = c;
                        //ppvy[0] = l; ppvy[1] = l+1; ppvy[2] = l+1; ppvy[3] = l;
                        isinFOVd = ((msldemc_imFOVmask[c+1][l]==2) && (msldemc_imFOVmask[c+1][l+1]==2) && (msldemc_imFOVmask[c][l+1]==2));
                        isinFOV = ((msldemc_imFOVmask[c+1][l]>0) && (msldemc_imFOVmask[c+1][l+1]>0) && (msldemc_imFOVmask[c][l+1]>0));
                    }
//                     isinFOVd = (msldemc_imFOVmaskd[cv1][lv1] && msldemc_imFOVmaskd[cv2][lv2] && msldemc_imFOVmaskd[cv3][lv3]);
//                     isinFOV = (msldemc_imFOVmask[cv1][lv1] && msldemc_imFOVmask[cv2][lv2] && msldemc_imFOVmask[cv3][lv3]);
                    
                    if(isinFOVd){
                        ppv1x = msldemc_imx[cv1][lv1];
                        ppv1y = msldemc_imy[cv1][lv1];
                        ppv2x = msldemc_imx[cv2][lv2];
                        ppv2y = msldemc_imy[cv2][lv2];
                        ppv3x = msldemc_imx[cv3][lv3];
                        ppv3y = msldemc_imy[cv3][lv3];
                        ppv1gx = msldemc_northing[lv1];
                        ppv1gy = msldemc_easting[cv1];
                        ppv1gz = msldemc_img[cv1][lv1];
                        ppv2gx = msldemc_northing[lv2];
                        ppv2gy = msldemc_easting[cv2];
                        ppv2gz = msldemc_img[cv2][lv2];
                        ppv3gx = msldemc_northing[lv3];
                        ppv3gy = msldemc_easting[cv3];
                        ppv3gz = msldemc_img[cv3][lv3];
                        ppv4gx = msldemc_northing[lv4];
                        ppv4gy = msldemc_easting[cv4];
                        ppv4gz = msldemc_img[cv4][lv4];

                        // define some plane parameters
                        pdv1x = ppv2x - ppv1x; pdv1y = ppv2y - ppv1y;
                        pdv2x = ppv3x - ppv1x; pdv2y = ppv3y - ppv1y;
                        detM = pdv1x*pdv2y - pdv1y*pdv2x;
                        Minv[0][0] = pdv2y/detM;
                        Minv[0][1] = -pdv2x/detM;
                        Minv[1][0] = -pdv1y/detM;
                        Minv[1][1] = pdv1x/detM;

                        
                        pdv1gx = ppv2gx - ppv1gx;
                        pdv1gy = ppv2gy - ppv1gy;
                        pdv1gz = ppv2gz - ppv1gz;
                        pdv2gx = ppv3gx - ppv1gx;
                        pdv2gy = ppv3gy - ppv1gy;
                        pdv2gz = ppv3gz - ppv1gz;
                        /* parameters for plane equations
                         * plane normal vector (pn)
                         * plane constant (pc)
                        */
                        pnx = pdv1gy*pdv2gz - pdv1gz*pdv2gy;
                        pny = pdv1gz*pdv2gx - pdv1gx*pdv2gz;
                        pnz = pdv1gx*pdv2gy - pdv1gy*pdv2gx;
                        lprm_nume = pnx*ppv1gx+pny*ppv1gy+pnz*ppv1gz;
                        
                        /* for pre-screening */
                        x_min = fmin(fmin(ppv1x,ppv2x),ppv3x);
                        y_min = fmin(fmin(ppv1y,ppv2y),ppv3y);
                        x_max = fmax(fmax(ppv1x,ppv2x),ppv3x);
                        y_max = fmax(fmax(ppv1y,ppv2y),ppv3y);
                        
                        x_min_int = (int32_T) floor(x_min);
                        y_min_int = (int32_T) floor(y_min);
                        x_max_int = (int32_T) ceil(x_max);
                        y_max_int = (int32_T) ceil(y_max);
                        
                        if(x_min_int<0){
                            x_min_int=0;   
                        }else if(x_min_int>S_imm1){
                            x_min_int=S_imm1;
                        }
                        if(x_max_int<1){
                            x_max_int=1;
                        }else if(x_max_int>S_im){
                            x_max_int=S_im;
                        }
                        
                        if(y_min_int<0){
                            y_min_int=0;
                        }else if(y_min_int>L_imm1){
                            y_min_int=L_imm1;
                        }
                        if(y_max_int<1){
                            y_max_int=1;
                        }else if(y_max_int>L_im){
                            y_max_int=L_im;
                        }
                        
                        for(xi=x_min_int;xi<x_max_int;xi++){
                            //printf("ll=%d/%d\n",ll,msldemc_lines);
                            for(yi=y_min_int;yi<y_max_int;yi++){
                                for (n=0;n<bin_count_im[xi][yi];n++){
                                    cc = bin_im_c[xi][yi][n];
                                    ll = bin_im_l[xi][yi][n];
                                    if(msldemc_imx[cc][ll] > x_min && msldemc_imx[cc][ll] < x_max && msldemc_imy[cc][ll] > y_min && msldemc_imy[cc][ll] < y_max)
                                    {
                                        /* evaluate line param */
                                        // tppvgx = msldemc_northing[ll];
                                        // tppvgy = msldemc_easting[cc];
                                        // tppvgz = msldemc_img[cc][ll];
                                        lprm = lprm_nume/(pnx*msldemc_northing[ll]+pny*msldemc_easting[cc]+pnz*msldemc_img[cc][ll]);
                                        if(lprm<1 && lprm>0){
                                            /* evaluate the test vector is inside the triangle. */
                                            //tppvx = msldemc_imx[cc][ll];
                                            //tppvy = msldemc_imxy[cc][ll];
                                            pipvx = msldemc_imx[cc][ll] - ppv1x; pipvy = msldemc_imy[cc][ll] - ppv1y; 
                                            pprm_sd = Minv[0][0]*pipvx+Minv[0][1]*pipvy;
                                            pprm_td = Minv[1][0]*pipvx+Minv[1][1]*pipvy;
                                            pprm_1std = 1-pprm_sd-pprm_td;
                                            if(pprm_sd>0 && pprm_td>0 && pprm_1std>0){
                                                if((cc==cv1 && ll==lv1) || (cc==cv2 && ll==lv2) || (cc==cv3 && ll==lv3)){
                                                } else {
                                                msldemc_inImage[cc][ll] = 0;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    /* free dynamically allocated memories */
    for(xi=0;xi<S_im;xi++){
        for(yi=0;yi<L_im;yi++){
            if(bin_count_im[xi][yi]>0){
                free(bin_im_c[xi][yi]);
                free(bin_im_l[xi][yi]);
            }
        }
    }
    
    free(bin_im_c);
    free(bin_im_l);
    free(bin_im_c_base);
    free(bin_im_l_base);
    
    free(bin_count_im);
    free(bin_count_im_base);
    
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

float** set_mxSingleMatrix(const mxArray *pmi){
    mwSize M,N;
    mwIndex j;
    float **pm;
    M = mxGetM(pmi); N = mxGetN(pmi);
    pm = (float **) mxMalloc(N*sizeof(float*));
    pm[0] = mxGetSingles(pmi);
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
    double **msldemc_img;
    double *msldemc_northing;
    double *msldemc_easting;
    double **msldemc_imx;
    double **msldemc_imy;
    int8_T **msldemc_imFOVmask;
    mwSize S_im,L_im;
    double **PmCx, **PmCy, **PmCz;
    
    int8_T **msldemc_inImage;
    
    mwIndex si,li;

    
    mwSize S_demc, L_demc;

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
    L_demc = mxGetM(prhs[0]);
    S_demc = mxGetN(prhs[0]);
    
    /* INPUT 1/2 msldem northing easting */
    msldemc_northing = mxGetDoubles(prhs[1]);
    msldemc_easting = mxGetDoubles(prhs[2]);
    
    /* INPUT 3/4 msldem imxy */
    msldemc_imx = set_mxDoubleMatrix(prhs[3]);    
    msldemc_imy = set_mxDoubleMatrix(prhs[4]);
    
    /* INPUT 5 msldem imFOV */
    msldemc_imFOVmask = set_mxInt8Matrix(prhs[5]);
    
    /* INPUT 7/8 image S_im, L_im */
    S_im = (mwSize) mxGetScalar(prhs[6]);
    L_im = (mwSize) mxGetScalar(prhs[7]);
    
    
    /* OUTPUT 0/1/2 north-east-elevation */
    plhs[0] = mxCreateNumericMatrix(L_demc,S_demc,mxINT8_CLASS,mxREAL);
    msldemc_inImage = set_mxInt8Matrix(plhs[0]);
    
    // Initialize matrices
    for(si=0;si<S_demc;si++){
        for(li=0;li<L_demc;li++){
            msldemc_inImage[si][li] = msldemc_imFOVmask[si][li];
        }
    }
    // printf("%d\n",S_im);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    find_hidden((int32_T) S_demc, (int32_T) L_demc, msldemc_img,
       msldemc_northing, msldemc_easting,
       msldemc_imx, msldemc_imy,
       msldemc_imFOVmask, 
       (int32_T) S_im, (int32_T) L_im,
       msldemc_inImage);
    
    /* free memories */
    mxFree(msldemc_img);
    mxFree(msldemc_imx);
    mxFree(msldemc_imy);
    mxFree(msldemc_imFOVmask);
    mxFree(msldemc_inImage);
    
}
