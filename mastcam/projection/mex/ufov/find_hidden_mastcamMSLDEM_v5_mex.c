/* =====================================================================
 * find_hidden_mastcamMSLDEM_v5_mex.c
 * Evaluate any pixels in MSLDEM image whether or not they exist in the 
 * 
 * INPUTS:
 * 0 msldemc_img           Float array [L_demc x S_demc]
 * 1 msldemc_northing      Double array [L_demc]
 * 2 msldemc_easting       Double array [S_demc]
 * 3 msldemc_imx           Doubles [L_demc x S_demc]
 * 4 msldemc_imy           Doubles [L_demc x S_demc]
 * 5 msldemc_imFOVmask     Boolean [L_demc x S_demc]
 * 6 msldemc_imFOVmaskd    Boolean [L_demc x S_demc]
 * 7 S_im                  int
 * 8 L_im                  int
 * 9 cam_C_pxl             Double [2 x 1] camera model center in the DEM pixel coordinate.
 *                         first element is north (line) and the second is east (sample) 
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
        bool **msldemc_imFOVmask, bool **msldemc_imFOVmaskd, 
        int32_T S_im, int32_T L_im, double *cam_C_pxl,
        bool **msldemc_inImage)
{
    int32_T c,l,cc,ll;
    int32_T cv1,cv2,cv3,cv4,lv1,lv2,lv3,lv4; /* v indicates vertex */
    int32_T L_demcm1,S_demcm1;
    double L_demcm1_dbl,S_demcm1_dbl;
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
    
    int32_T cstrtc,cendc,lstrtc,lendc;
    int32_T cstrtc_ii,cendc_ii,lstrtc_ii,lendc_ii;
    double y_xend;
    double vec_c[3],vec_l[3],coef_c[3];
    double ca[3],cb[3],slp_cl[3],slp_lc[3];
    double l_val,l_min,l_max;
    double c_val,c_min,c_max;
    int32_T i_min,i_max;
    double ca_min,ca_max,mba_min,mba_max;
    int32_T ll_min,ll_max;
    double cb_min,cb_max,mab_min,mab_max;
    int32_T cc_min,cc_max,ii;
    
    double lv1_dbl,ll_dbl;
    double cv1_dbl,cc_dbl;
    
    //printf("%d\n",msldemc_samples);
    L_demcm1 = msldemc_lines-1;
    S_demcm1 = msldemc_samples-1;
    L_demcm1_dbl = (double) L_demcm1_dbl;
    S_demcm1_dbl = (double) S_demcm1_dbl;
    //printf("%d\n",msldemc_samples);
    
                
    for(l=10000;l<10100;l++){
        //printf("l=%d\n",l);
        // decide the first and last indexes to be assessed.
        for(c=0;c<S_demcm1;c++){
            // process if 
            //printf("c=%d,mask_lc = %d,mask_lp1c = %d\n",c,msldemc_imFOVmask[c][l],msldemc_imFOVmask[c][l+1]);
            // printf("c=%d/%d\n",c,S_demcm1);
            if(msldemc_imFOVmask[c][l] || msldemc_imFOVmask[c][l+1]){
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
                        isinFOVd = (msldemc_imFOVmaskd[c][l] && msldemc_imFOVmaskd[c+1][l] && msldemc_imFOVmaskd[c][l+1]);
                        isinFOV = (msldemc_imFOVmask[c][l] && msldemc_imFOVmask[c+1][l] && msldemc_imFOVmask[c][l+1]);
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
                        isinFOVd = (msldemc_imFOVmaskd[c+1][l] && msldemc_imFOVmaskd[c+1][l+1] && msldemc_imFOVmaskd[c][l+1]);
                        isinFOV = (msldemc_imFOVmask[c+1][l] && msldemc_imFOVmask[c+1][l+1] && msldemc_imFOVmask[c][l+1]);
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
                        
                        /* line parameters 
                         * coef_a * l + coef_b * c = coef_c
                         * vec_c * l - vec_l * c = coef_c
                         */
                        vec_c[0] = (double) cv1 - cam_C_pxl[1];
                        vec_c[1] = (double) cv2 - cam_C_pxl[1];
                        vec_c[2] = (double) cv3 - cam_C_pxl[1];
                        
                        vec_l[0] = (double) lv1 - cam_C_pxl[0];
                        vec_l[1] = (double) lv2 - cam_C_pxl[0];
                        vec_l[2] = (double) lv3 - cam_C_pxl[0];
                        
                        coef_c[0] = vec_c[0]*cam_C_pxl[0] - vec_l[0]*cam_C_pxl[1];
                        coef_c[1] = vec_c[1]*cam_C_pxl[0] - vec_l[1]*cam_C_pxl[1];
                        coef_c[2] = vec_c[2]*cam_C_pxl[0] - vec_l[2]*cam_C_pxl[1];
                        
                        ca[0] = coef_c[0]/vec_c[0];
                        ca[1] = coef_c[1]/vec_c[1];
                        ca[2] = coef_c[2]/vec_c[2];
                        
                        cb[0] = -coef_c[0]/vec_l[0];
                        cb[1] = -coef_c[1]/vec_l[1];
                        cb[2] = -coef_c[2]/vec_l[2];
                        
                        slp_cl[0] = vec_c[0]/vec_l[0];
                        slp_cl[1] = vec_c[1]/vec_l[1];
                        slp_cl[2] = vec_c[2]/vec_l[2];
                        
                        slp_lc[0] = vec_l[0]/vec_c[0];
                        slp_lc[1] = vec_l[1]/vec_c[1];
                        slp_lc[2] = vec_l[2]/vec_c[2];

                        /* evaluate effective range of line and sample */
                        /* original range [0 msldemc_lines) [0 msldemc_samples) */
                        cstrtc = 0; cendc = msldemc_samples;
                        lstrtc = 0; lendc = msldemc_lines;
                        for(ii=0;ii<3;ii++){
                            y_xend = cb[ii] + slp_cl[ii]*L_demcm1_dbl;
                            if(slp_cl>0){
                                if(cb[ii]>0 && cb[ii]<S_demcm1_dbl){
                                    lstrtc_ii = 0;
                                    cstrtc_ii = (int32_T) ceil(cb[ii]); /* y intercept */
                                }else{
                                    lstrtc_ii = (int32_T) ceil(ca[ii]);  /* x intercept */
                                    cstrtc_ii = 0;
                                }
                                if(y_xend > 0 && y_xend < S_demcm1_dbl){
                                    lendc_ii = msldemc_lines;
                                    cendc_ii = (int32_T) ceil(y_xend); /* y intercept */
                                }else{
                                    lendc_ii = (int32_T) ceil(ca[ii] + slp_lc[ii]*S_demcm1_dbl);  /* x intercept */
                                    cendc_ii = msldemc_samples;
                                }
                            } else {
                                if(cb[ii]>0 && cb[ii]<S_demcm1_dbl){
                                    lstrtc_ii = 0;
                                    cendc_ii = (int32_T) ceil(cb[ii]); /* y intercept */
                                }else{
                                    lstrtc_ii = (int32_T) ceil(ca[ii]+slp_lc[ii]*S_demcm1_dbl);  /* x intercept */
                                    cendc_ii = msldemc_samples;
                                }
                                if(y_xend > 0 && y_xend < S_demcm1_dbl){
                                    lendc_ii = msldemc_lines;
                                    cstrtc_ii = (int32_T) ceil(y_xend); /* y intercept */
                                }else{
                                    lendc_ii = (int32_T) ceil(ca[ii]);  /* x intercept */
                                    cstrtc_ii = 0;
                                }
                            }
                            if(cstrtc_ii < cstrtc)
                                cstrtc = cstrtc_ii;
                            if(cendc_ii > cendc)
                                cendc = cendc_ii;
                            if(lstrtc_ii<lstrtc)
                                lstrtc = lstrtc_ii;
                            if(lendc_ii>lendc)
                                lendc = lendc_ii;
                        }
                        
                        /* further restricting the range */
                        if(vec_l[0]<0 && vec_l[1]<0 && vec_l[2]<0){
                            lendc = l;
                        } else if(vec_l[0]>0 && vec_l[1]>0 && vec_l[2]>0){
                            lstrtc = l;
                        }
                        if(vec_c[0]<0 && vec_c[1]<0 && vec_c[2]<0){
                            cendc = c;
                        } else if(vec_c[0]>0 && vec_c[1]>0 && vec_c[2]>0){
                            cstrtc = c;
                        }
                        
                        /* different processing for column mode and line mode */
                        if ((lendc-lstrtc) > (cendc-cstrtc)){
                            
                            cv1_dbl = (double) cv1;
                            l_val = ca[0]+slp_lc[0]*cv1_dbl;
                            i_min = 0; i_max = 0; l_min = l_val; l_max = l_val;
                            for(ii=1;ii<3;ii++){
                                l_val = ca[ii]+slp_lc[ii]*cv1_dbl;
                                if(l_val<l_min){
                                    i_min = ii;
                                }
                                else if(l_val>l_max){
                                    i_max = ii;
                                }
                            }
                            ca_min  = ca[i_min];
                            mba_min = slp_lc[i_min];
                            ca_max  = ca[i_max];
                            mba_max = slp_lc[i_max];
                            
                            
                            for(cc=cstrtc;cc<cendc;cc++){
                                //printf("ll=%d/%d\n",ll,msldemc_lines);
                                cc_dbl = (double) cc;
                                ll_min = (int32_T) ceil(ca_min + mba_min*cc_dbl);
                                ll_max = (int32_T) ceil(ca_max + mba_max*cc_dbl);
                                if(ll_min < 0)
                                    ll_min = 0;
                                else if(ll_min > msldemc_lines)
                                    ll_min = msldemc_lines;

                                if(ll_max < 0)
                                    ll_max = 0;
                                else if(ll_max > msldemc_lines)
                                    ll_max = msldemc_lines;
                                
                                for(ll=ll_min;ll<ll_max;ll++){
                                        /* If  */
                                    tppvy = msldemc_imy[cc][ll];
                                    if(tppvy > y_min && tppvy < y_max)
                                    {
                                        /* evaluate line param */
                                        // tppvx = msldemc_imx[cc][ll];
                                        // tppvgx = msldemc_northing[ll];
                                        // tppvgy = msldemc_easting[cc];
                                        // tppvgz = msldemc_img[cc][ll];
                                        pipvx = msldemc_imx[cc][ll] - ppv1x; pipvy = tppvy - ppv1y; 
                                        pprm_sd = Minv[0][0]*pipvx+Minv[0][1]*pipvy;
                                        pprm_td = Minv[1][0]*pipvx+Minv[1][1]*pipvy;
                                        pprm_1std = 1-pprm_sd-pprm_td;
                                        if(pprm_sd>0 && pprm_td>0 && pprm_1std>0){
                                            lprm = lprm_nume/(pnx*msldemc_northing[ll]+pny*msldemc_easting[cc]+pnz*msldemc_img[cc][ll]);
                                            if(lprm<1 && lprm>0){
                                                /* evaluate the test vector is inside the triangle. */
                                                //tppvx = msldemc_imx[cc][ll];
                                                //tppvy = msldemc_imxy[cc][ll];
                                                msldemc_inImage[cc][ll] = false;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        else{
                            
                            lv1_dbl = (double) lv1;
                            c_val = cb[0]+slp_cl[0]*lv1_dbl;
                            i_min = 0; i_max = 0; c_min = c_val; c_max = c_val;
                            for(ii=1;ii<3;ii++){
                                c_val = cb[ii]+slp_cl[ii]*lv1_dbl;
                                if(c_val<c_min){
                                    i_min = ii;
                                }
                                else if(c_val>c_max){
                                    i_max = ii;
                                }
                            }

                            cb_min  = cb[i_min];
                            mab_min = slp_cl[i_min];
                            cb_max  = cb[i_max];
                            mab_max = slp_cl[i_max];
                            
                            /* evaluate line loop mode or column loop mode */
                            for(ll=lstrtc;ll<lendc;ll++){
                                //printf("ll=%d/%d\n",ll,msldemc_lines);
                                ll_dbl = (double) ll;
                                cc_min = (int32_T) ceil(cb_min + mab_min*ll_dbl);
                                cc_max = (int32_T) ceil(cb_max + mab_max*ll_dbl);
                                if(cc_min < 0)
                                    cc_min = 0;
                                else if(cc_min > msldemc_samples)
                                    cc_min = msldemc_samples;

                                if(cc_max < 0)
                                    cc_max = 0;
                                else if(cc_max > msldemc_samples)
                                    cc_max = msldemc_samples;
                                
                                //if(c==8000 && l==10001 && ll%5000==0)
                                //printf("ll %d, cc:%d - %d\n",ll,cc_min,cc_max);
                                for(cc=cc_min;cc<cc_max;cc++){
                                    /* If  */
                                    tppvy = msldemc_imy[cc][ll];
                                    if(tppvy > y_min && tppvy < y_max)
                                    {
                                        /* evaluate line param */
                                        // tppvx = msldemc_imx[cc][ll];
                                        // tppvgx = msldemc_northing[ll];
                                        // tppvgy = msldemc_easting[cc];
                                        // tppvgz = msldemc_img[cc][ll];
                                        pipvx = msldemc_imx[cc][ll] - ppv1x; pipvy = tppvy - ppv1y; 
                                        pprm_sd = Minv[0][0]*pipvx+Minv[0][1]*pipvy;
                                        pprm_td = Minv[1][0]*pipvx+Minv[1][1]*pipvy;
                                        pprm_1std = 1-pprm_sd-pprm_td;
                                        if(pprm_sd>0 && pprm_td>0 && pprm_1std>0){
                                            lprm = lprm_nume/(pnx*msldemc_northing[ll]+pny*msldemc_easting[cc]+pnz*msldemc_img[cc][ll]);
                                            if(lprm<1 && lprm>0){
                                                /* evaluate the test vector is inside the triangle. */
                                                //tppvx = msldemc_imx[cc][ll];
                                                //tppvy = msldemc_imxy[cc][ll];
                                                msldemc_inImage[cc][ll] = false;
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
    bool **msldemc_imFOVmask;
    bool **msldemc_imFOVmaskd;
    mwSize S_im,L_im;
    double *cam_C_pxl;
    double **PmCx, **PmCy, **PmCz;
    
    bool **msldemc_inImage;
    
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
    msldemc_imFOVmask = set_mxLogicalMatrix(prhs[5]);

    /* INPUT 6 msldem imFOV */
    msldemc_imFOVmaskd = set_mxLogicalMatrix(prhs[6]);
    
    /* INPUT 7/8 image S_im, L_im */
    S_im = (mwSize) mxGetScalar(prhs[7]);
    L_im = (mwSize) mxGetScalar(prhs[8]);
    
    /* INPUT 9 camera model */
    cam_C_pxl = mxGetDoubles(prhs[9]);
    
    
    /* OUTPUT 0/1/2 north-east-elevation */
    plhs[0] = mxCreateLogicalMatrix(L_demc,S_demc);
    msldemc_inImage = set_mxLogicalMatrix(plhs[0]);
    
    // Initialize matrices
    for(si=0;si<S_demc;si++){
        for(li=0;li<L_demc;li++){
            msldemc_inImage[si][li] = msldemc_imFOVmask[si][li];
        }
    }
    //printf("%d\n",S_im);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    find_hidden((int32_T) S_demc, (int32_T) L_demc, msldemc_img,
       msldemc_northing, msldemc_easting,
       msldemc_imx, msldemc_imy,
       msldemc_imFOVmask,msldemc_imFOVmaskd, 
       (int32_T) S_im, (int32_T) L_im, cam_C_pxl,
       msldemc_inImage);
    
    /* free memories */
    mxFree(msldemc_img);
    mxFree(msldemc_imx);
    mxFree(msldemc_imy);
    mxFree(msldemc_imFOVmask);
    mxFree(msldemc_imFOVmaskd);
    mxFree(msldemc_inImage);
    
}
