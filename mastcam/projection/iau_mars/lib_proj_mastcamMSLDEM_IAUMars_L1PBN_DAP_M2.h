
#ifndef LIB_PROJ_MASTCAMMSLDEM_IAUMARS_H
#define LIB_PROJ_MASTCAMMSLDEM_IAUMARS_H

#include <stdint.h>
#include "io64.h"
#include "math.h"
#include "matrix.h"
#include <string.h>
#include <stdio.h>

#include <stdlib.h>
#include "envi.h"
#include "cahvor.h"
#include "mex_create_array.h"

/* main computation routine */
void mask_obstructed_pts_in_msldemt_using_msldemc_iaumars(
        double **msldemc_img,
        int32_T msldemc_samples, int32_T msldemc_lines,
        double *msldemc_latitude, double *msldemc_longitude, int8_T **msldemc_imFOVmask,
        int32_T msldemt_samples, int32_T msldemt_lines,
        double *msldemt_latitude, double *msldemt_longitude, int8_T **msldemt_inImage,
        int32_T **bin_count_im, int32_T ***bin_im_c, int32_T ***bin_im_l,
        double ***bin_imx, double ***bin_imy, double ***bin_rad,
        int32_T count_napmc, int32_T *c_napmc, int32_T *l_napmc, double *rad_napmc,
        int32_T S_im, int32_T L_im, CAHV_MODEL cahv_mdl)
{
    int32_T c,l,cc,ll;
    int32_T cv1,cv2,cv3,lv1,lv2,lv3; /* v indicates vertex */
    int32_T L_demcm1,S_demcm1;
    int16_T ti; /* triangle index */
    double ppv1x,ppv1y,ppv2x,ppv2y,ppv3x,ppv3y; /* Plane Position Vectors */
    double ppv1gx,ppv1gy,ppv1gz,ppv2gx,ppv2gy,ppv2gz,ppv3gx,ppv3gy,ppv3gz;
    double pdv1x,pdv1y,pdv2x,pdv2y;
    double pdv1gx,pdv1gy,pdv1gz,pdv2gx,pdv2gy,pdv2gz;
    double detM;
    double M[2][2];
    double Minv[2][2];
    double Minvp[2][3];
    double x_min,y_min,x_max,y_max;
    double pipvx,pipvy;
    double pipvgx,pipvgy,pipvgz;
    double pdv1z,pdv2z;
    double pprm_sd,pprm_td,pprm_1std; /* plane parameter for projected image plane */
    double pipvgppv1x,pipvgppv1y,pipvgppv1z;
    double pprm_s,pprm_t,pprm_1st;
    bool isinFOVd,isinFOV;
    
    double pnx,pny,pnz; /* Plane Normal vectors */
    double lprm_nume;
    double lprm; /* line parameters */
    
    int32_T xi,yi;
    int32_T x_min_int,x_max_int,y_min_int,y_max_int;
    int32_T n;
    int32_T S_imm1,L_imm1;
    
    double apmc,pmcx,pmcy,pmcz;
    double ppvgx,ppvgy,ppvgz,ppvx,ppvy;
    double *cam_C,*cam_A,*cam_H,*cam_V;
    
    double *cos_clon, *sin_clon, *cos_tlon, *sin_tlon;
    double *cos_tlat, *sin_tlat;
    double radius_tmp;
    double cos_clatl, sin_clatl, cos_clatlp1, sin_clatlp1, cos_tlatl, sin_tlatl;
    double x_iaumars, y_iaumars, z_iaumars;
    
    cos_clon = (double*) malloc(sizeof(double) * (size_t) msldemc_samples);
    sin_clon = (double*) malloc(sizeof(double) * (size_t) msldemc_samples);
    for(c=0;c<msldemc_samples;c++){
        cos_clon[c] = cos(msldemc_longitude[c]);
        sin_clon[c] = sin(msldemc_longitude[c]);
    }
    
    if(msldemt_latitude==NULL){
        msldemt_latitude = msldemc_latitude;
    }
    if(msldemt_longitude==NULL){
        msldemt_longitude = msldemc_longitude;
        cos_tlon = cos_clon; sin_tlon = sin_clon;
    } else {
        cos_tlon = (double*) malloc(sizeof(double) * (size_t) msldemt_samples);
        sin_tlon = (double*) malloc(sizeof(double) * (size_t) msldemt_samples);
        for(c=0;c<msldemt_samples;c++){
            cos_tlon[c] = cos(msldemt_longitude[c]);
            sin_tlon[c] = sin(msldemt_longitude[c]);
        }
    }
    
    cos_tlat = (double*) malloc(sizeof(double) * (size_t) msldemt_lines);
    sin_tlat = (double*) malloc(sizeof(double) * (size_t) msldemt_lines);
    for(l=0;l<msldemt_lines;l++){
        cos_tlat[l] = cos(msldemt_latitude[l]);
        sin_tlat[l] = sin(msldemt_latitude[l]);
    }
    
    cam_C = cahv_mdl.C; cam_A = cahv_mdl.A; cam_H = cahv_mdl.H; cam_V = cahv_mdl.V;
    
    L_demcm1 = msldemc_lines-1;
    S_demcm1 = msldemc_samples-1;
    S_imm1 = S_im - 1;
    L_imm1 = L_im - 1;
    
    /*********************************************************************/
    
    
    /*********************************************************************/
    /*** Pre-binning of the msldem pixels ********************************/
    /* create an bin image counting the number of demc pixels that falls
     * within the 
     */
    /* Dynamic memory allocations */
    
    
    
    
    /* Main Loop *********************************************************/
    
//     find_hidden_main_loop(msldemc_samples, msldemc_lines, msldemc_imFOVmask,
//         msldemc_northing, msldemc_easting, msldemc_img,
//         cam_A, cam_H, cam_V,
//         S_im, L_im, bin_count_im, bin_im_c, bin_im_l, bin_imx, bin_imy,
//         msldemc_northing, msldemc_easting, msldemc_img, msldemc_inImage);
    
    
    for(l=0;l<L_demcm1;l++){
        // decide the first and last indexes to be assessed.
        cos_clatl   = cos(msldemc_latitude[l]);
        sin_clatl   = sin(msldemc_latitude[l]);
        cos_clatlp1 = cos(msldemc_latitude[l+1]);
        sin_clatlp1 = sin(msldemc_latitude[l+1]);
        for(c=0;c<S_demcm1;c++){
            // process if 
            //printf("c=%d,mask_lc = %d,mask_lp1c = %d\n",c,msldemc_imFOVmask[c][l],msldemc_imFOVmask[c][l+1]);
            // printf("c=%d/%d\n",c,S_demcm1);
            if((msldemc_imFOVmask[c][l]>0) || (msldemc_imFOVmask[c][l+1]>0)){
                //printf("c=%d\n",c);
                for(ti=0;ti<2;ti++){
                    if(ti==0){
                        radius_tmp = msldemc_img[c][l];
                        ppv1gx = radius_tmp * cos_clatl * cos_clon[c];
                        ppv1gy = radius_tmp * cos_clatl * sin_clon[c];
                        ppv1gz = radius_tmp * sin_clatl;
                        radius_tmp = msldemc_img[c+1][l];
                        ppv2gx = radius_tmp * cos_clatl * cos_clon[c+1];
                        ppv2gy = radius_tmp * cos_clatl * sin_clon[c+1];
                        ppv2gz = radius_tmp * sin_clatl;
                        radius_tmp = msldemc_img[c][l+1];
                        ppv3gx = radius_tmp * cos_clatlp1 * cos_clon[c];
                        ppv3gy = radius_tmp * cos_clatlp1 * sin_clon[c];
                        ppv3gz = radius_tmp * sin_clatlp1;
//                         ppv1gx = msldemc_xmc[l];
//                         ppv1gy = msldemc_ymc[c];
//                         ppv1gz = ((double) -elevl[c]) - cam_C[2];
//                         ppv2gx = msldemc_xmc[l];
//                         ppv2gy = msldemc_ymc[c+1];
//                         ppv2gz = ((double) -elevl[c+1]) - cam_C[2];
//                         ppv3gx = msldemc_xmc[l+1];
//                         ppv3gy = msldemc_ymc[c];
//                         ppv3gz = ((double) -elevlp1[c]) - cam_C[2];
                        isinFOVd = ((msldemc_imFOVmask[c][l]>1) && (msldemc_imFOVmask[c+1][l]>1) && (msldemc_imFOVmask[c][l+1]>1));
                        isinFOV = ((msldemc_imFOVmask[c][l]>0) && (msldemc_imFOVmask[c+1][l]>0) && (msldemc_imFOVmask[c][l+1]>0));
                        cv1 = c;   lv1 = l;
                        cv2 = c+1; lv2 = l;
                        cv3 = c;   lv3 = l+1;
                        // cv4 = c+1; lv4 = l+1;
                    }
                    else{
                        radius_tmp = msldemc_img[c+1][l];
                        ppv1gx = radius_tmp * cos_clatl * cos_clon[c+1];
                        ppv1gy = radius_tmp * cos_clatl * sin_clon[c+1];
                        ppv1gz = radius_tmp * sin_clatl;
                        radius_tmp = msldemc_img[c+1][l+1];
                        ppv2gx = radius_tmp * cos_clatlp1 * cos_clon[c+1];
                        ppv2gy = radius_tmp * cos_clatlp1 * sin_clon[c+1];
                        ppv2gz = radius_tmp * sin_clatlp1;
                        radius_tmp = msldemc_img[c][l+1];
                        ppv3gx = radius_tmp * cos_clatlp1 * cos_clon[c];
                        ppv3gy = radius_tmp * cos_clatlp1 * sin_clon[c];
                        ppv3gz = radius_tmp * sin_clatlp1;
//                         ppv1gx = msldemc_xmc[l];
//                         ppv1gy = msldemc_ymc[c+1];
//                         ppv1gz = ((double) -elevl[c+1]) - cam_C[2];
//                         ppv2gx = msldemc_xmc[l+1];
//                         ppv2gy = msldemc_ymc[c+1];
//                         ppv2gz = ((double) -elevlp1[c+1]) - cam_C[2];
//                         ppv3gx = msldemc_xmc[l+1];
//                         ppv3gy = msldemc_ymc[c];
//                         ppv3gz = ((double) -elevlp1[c]) - cam_C[2];
                        cv1 = c+1; lv1 = l;
                        cv2 = c+1; lv2 = l+1;
                        cv3 = c;   lv3 = l+1;
                        // cv4 = c;   lv4 = l;
                        isinFOVd = ((msldemc_imFOVmask[c+1][l]>1) && (msldemc_imFOVmask[c+1][l+1]>1) && (msldemc_imFOVmask[c][l+1]>1));
                        isinFOV = ((msldemc_imFOVmask[c+1][l]>0) && (msldemc_imFOVmask[c+1][l+1]>0) && (msldemc_imFOVmask[c][l+1]>0));
                    }
                    
                    if(isinFOVd){
                        /* Evaluate the projection */
                        // pmcx = ppv1gx; pmcy = ppv1gy; pmcz = ppv1gz;
                        pmcx  = ppv1gx - cam_C[0];
                        pmcy  = ppv1gy - cam_C[1];
                        pmcz  = ppv1gz - cam_C[2];
                        apmc  =  pmcx*cam_A[0] + pmcy*cam_A[1] + pmcz*cam_A[2];
                        ppv1x = (pmcx*cam_H[0] + pmcy*cam_H[1] + pmcz*cam_H[2])/apmc;
                        ppv1y = (pmcx*cam_V[0] + pmcy*cam_V[1] + pmcz*cam_V[2])/apmc;
                        
                        pmcx  = ppv2gx - cam_C[0];
                        pmcy  = ppv2gy - cam_C[1];
                        pmcz  = ppv2gz - cam_C[2];
                        apmc  =  pmcx*cam_A[0] + pmcy*cam_A[1] + pmcz*cam_A[2];
                        ppv2x = (pmcx*cam_H[0] + pmcy*cam_H[1] + pmcz*cam_H[2])/apmc;
                        ppv2y = (pmcx*cam_V[0] + pmcy*cam_V[1] + pmcz*cam_V[2])/apmc;
                        
                        
                        // pmcx = ppv3gx; pmcy = ppv3gy; pmcz = ppv3gz;
                        pmcx  = ppv3gx - cam_C[0];
                        pmcy  = ppv3gy - cam_C[1];
                        pmcz  = ppv3gz - cam_C[2];
                        apmc  =  pmcx*cam_A[0] + pmcy*cam_A[1] + pmcz*cam_A[2];
                        ppv3x = (pmcx*cam_H[0] + pmcy*cam_H[1] + pmcz*cam_H[2])/apmc;
                        ppv3y = (pmcx*cam_V[0] + pmcy*cam_V[1] + pmcz*cam_V[2])/apmc;
                        
                        //printf("c=%d\n",l);

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
                        lprm_nume = pnx*(ppv1gx-cam_C[0])+pny*(ppv1gy-cam_C[1])+pnz*(ppv1gz-cam_C[2]);
                        
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
                            for(yi=y_min_int;yi<y_max_int;yi++){
                                for (n=0;n<bin_count_im[xi][yi];n++){
                                    cc = bin_im_c[xi][yi][n];
                                    ll = bin_im_l[xi][yi][n];
                                    /* evaluate line param */
                                    radius_tmp = bin_rad[xi][yi][n];
                                    x_iaumars  = radius_tmp * cos_tlat[ll] * cos_tlon[cc];
                                    y_iaumars  = radius_tmp * cos_tlat[ll] * sin_tlon[cc];
                                    z_iaumars  = radius_tmp * sin_tlat[ll];
                                    pmcx = x_iaumars - cam_C[0];
                                    pmcy = y_iaumars - cam_C[1];
                                    pmcz = z_iaumars - cam_C[2];
                                    lprm = lprm_nume/(pnx*pmcx+pny*pmcy+pnz*pmcz);
                                    if(lprm<1 && lprm>0){
                                        /* evaluate the test vector is inside the triangle. */
                                        ppvx = bin_imx[xi][yi][n];
                                        ppvy = bin_imy[xi][yi][n];
                                        pipvx = ppvx - ppv1x; pipvy = ppvy - ppv1y; 
                                        pprm_sd = Minv[0][0]*pipvx+Minv[0][1]*pipvy;
                                        pprm_td = Minv[1][0]*pipvx+Minv[1][1]*pipvy;
                                        pprm_1std = 1-pprm_sd-pprm_td;
                                        if(pprm_sd>0 && pprm_td>0 && pprm_1std>0){
                                            if((cc==cv1 && ll==lv1) || (cc==cv2 && ll==lv2) || (cc==cv3 && ll==lv3)){
                                            } else {
                                            msldemt_inImage[cc][ll] = 0;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    } else if(isinFOV){
                        pdv1x = ppv2gx - ppv1gx;
                        pdv1y = ppv2gy - ppv1gy;
                        pdv1z = ppv2gz - ppv1gz;
                        pdv2x = ppv3gx - ppv1gx;
                        pdv2y = ppv3gy - ppv1gy;
                        pdv2z = ppv3gz - ppv1gz;
                        pnx = pdv1y*pdv2z - pdv1z*pdv2y;
                        pny = pdv1z*pdv2x - pdv1x*pdv2z;
                        pnz = pdv1x*pdv2y - pdv1y*pdv2x;
                        
                        /* Get Plane parameters */
                        M[0][0] = pdv1x*pdv1x + pdv1y*pdv1y + pdv1z*pdv1z;
                        M[0][1] = pdv1x*pdv2x + pdv1y*pdv2y + pdv1z*pdv2z;
                        M[1][0] = M[0][1];
                        M[1][1] = pdv2x*pdv2x + pdv2y*pdv2y + pdv2z*pdv2z;
                        detM = M[0][0]*M[1][1] - M[0][1]*M[0][1];
                        Minv[0][0] = M[1][1]/detM;
                        Minv[0][1] = -M[0][1]/detM;
                        Minv[1][0] = -M[1][0]/detM;
                        Minv[1][1] = M[0][0]/detM;
                        Minvp[0][0] = Minv[0][0]*pdv1x+Minv[0][1]*pdv2x;
                        Minvp[0][1] = Minv[0][0]*pdv1y+Minv[0][1]*pdv2y;
                        Minvp[0][2] = Minv[0][0]*pdv1z+Minv[0][1]*pdv2z;
                        Minvp[1][0] = Minv[1][0]*pdv1x+Minv[1][1]*pdv2x;
                        Minvp[1][1] = Minv[1][0]*pdv1y+Minv[1][1]*pdv2y;
                        Minvp[1][2] = Minv[1][0]*pdv1z+Minv[1][1]*pdv2z;
                        
                        /* parameters for plane equations
                         * plane normal vector (pn)
                         * plane constant (pc)
                        */
                        pnx = pdv1gy*pdv2gz - pdv1gz*pdv2gy;
                        pny = pdv1gz*pdv2gx - pdv1gx*pdv2gz;
                        pnz = pdv1gx*pdv2gy - pdv1gy*pdv2gx;
                        lprm_nume = pnx*(ppv1gx-cam_C[0])+pny*(ppv1gy-cam_C[1])+pnz*(ppv1gz-cam_C[2]);
                        
                        for(xi=0;xi<S_im;xi++){
                            for(yi=0;yi<L_im;yi++){
                                for (n=0;n<bin_count_im[xi][yi];n++){
                                    cc = bin_im_c[xi][yi][n];
                                    ll = bin_im_l[xi][yi][n];
                                    /* evaluate line param */
                                    radius_tmp = bin_rad[xi][yi][n];
                                    x_iaumars  = radius_tmp * cos_tlat[ll] * cos_tlon[cc];
                                    y_iaumars  = radius_tmp * cos_tlat[ll] * sin_tlon[cc];
                                    z_iaumars  = radius_tmp * sin_tlat[ll];
                                    pmcx = x_iaumars - cam_C[0];
                                    pmcy = y_iaumars - cam_C[1];
                                    pmcz = z_iaumars - cam_C[2];
                                    lprm = lprm_nume/(pnx*pmcx+pny*pmcy+pnz*pmcz);
                                    if(lprm<1 && lprm>0){
                                        /* evaluate the test vector is inside the triangle. */
                                        pipvgppv1x = lprm*pmcx - ppv1gx;
                                        pipvgppv1y = lprm*pmcy - ppv1gy;
                                        pipvgppv1z = lprm*pmcz - ppv1gy;
                                        pprm_s = Minvp[0][0]*pipvgppv1x+Minvp[0][1]*pipvgppv1y+Minvp[0][2]*pipvgppv1z;
                                        pprm_t = Minvp[1][0]*pipvgppv1x+Minvp[1][1]*pipvgppv1y+Minvp[1][2]*pipvgppv1z;
                                        pprm_1st = 1 - pprm_s - pprm_t;
                                        if(pprm_s>0 && pprm_t>0 && pprm_1st>0){
                                            if((cc==cv1 && ll==lv1) || (cc==cv2 && ll==lv2) || (cc==cv3 && ll==lv3)){
                                            } else {
                                            msldemt_inImage[cc][ll] = 0;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        
                        /* test vectors with napmc<0 */
                        if(count_napmc>0){
                            for(n=0;n<count_napmc;n++){
                                cc = c_napmc[n]; ll = l_napmc[n];
                                radius_tmp = rad_napmc[n];
                                x_iaumars  = radius_tmp * cos_tlat[ll] * cos_tlon[cc];
                                y_iaumars  = radius_tmp * cos_tlat[ll] * sin_tlon[cc];
                                z_iaumars  = radius_tmp * sin_tlat[ll];
                                pmcx = x_iaumars - cam_C[0];
                                pmcy = y_iaumars - cam_C[1];
                                pmcz = z_iaumars - cam_C[2];
                                lprm = lprm_nume/(pnx*pmcx+pny*pmcy+pnz*pmcz);
                                if(lprm<1 && lprm>0){
                                    /* evaluate the test vector is inside the triangle. */
                                    pipvgppv1x = lprm*pmcx - ppv1gx;
                                    pipvgppv1y = lprm*pmcy - ppv1gy;
                                    pipvgppv1z = lprm*pmcz - ppv1gy;
                                    pprm_s = Minvp[0][0]*pipvgppv1x+Minvp[0][1]*pipvgppv1y+Minvp[0][2]*pipvgppv1z;
                                    pprm_t = Minvp[1][0]*pipvgppv1x+Minvp[1][1]*pipvgppv1y+Minvp[1][2]*pipvgppv1z;
                                    pprm_1st = 1 - pprm_s - pprm_t;
                                    if(pprm_s>0 && pprm_t>0 && pprm_1st>0){
                                        if((cc==cv1 && ll==lv1) || (cc==cv2 && ll==lv2) || (cc==cv3 && ll==lv3)){
                                        } else {
                                        msldemt_inImage[cc][ll] = 0;
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
    free(cos_clon);
    free(sin_clon);
    if(cos_tlon)
        free(cos_tlon);
    if(sin_tlon)
        free(sin_tlon);
    free(cos_tlat);
    free(sin_tlat);
    
}



#endif