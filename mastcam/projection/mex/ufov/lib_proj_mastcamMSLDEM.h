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

#ifndef LIB_PROH_MASTCAMMSLDEM_H
#define LIB_PROH_MASTCAMMSLDEM_H

void bin_msldemt_xyz_wAHVint(int32_T S_im, int32_T L_im,
        double *cam_A, double *cam_H, double *cam_V,
        int32_T msldemt_samples, int32_T msldemt_lines,
        double *msldemt_x, double *msldemt_y, double **msldemt_z,
        int8_T **msldemt_inImage,
        int32_T **bin_count_im,
        int32_T ***bin_im_c, int32_T ***bin_im_l, double ***bin_imx, double ***bin_imy,
        int32_T *count_napmc, int32_T **c_napmc, int32_T **l_napmc)
{
    int32_T xi,yi,c,l;
    int32_T S_imm1,L_imm1;
    double pmcx,pmcy,pmcz,apmc,ppvx,ppvy;
    
    S_imm1 = S_im - 1;
    L_imm1 = L_im - 1;
    
    *count_napmc = 0;
    
    /* Initialize the counting matrix */
    for(xi=0;xi<S_im;xi++){
        for(yi=0;yi<L_im;yi++){
            bin_count_im[xi][yi] = 0;
        }
    }
    // printf("%d\n",msldemc_samples);
    for(l=0;l<msldemt_lines;l++){
        for(c=0;c<msldemt_samples;c++){
            if(msldemt_inImage[c][l]>1){
                pmcx = msldemt_x[l];
                pmcy = msldemt_y[c];
                pmcz = msldemt_z[c][l];
                // pmcx = ppvgx-cam_C[0]; pmcy = ppvgy-cam_C[1]; pmcz = ppvgz-cam_C[2];
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
            } else if(msldemt_inImage[c][l]==1) {
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
            } else {
                bin_im_c[xi][yi] = NULL;
                bin_im_l[xi][yi] = NULL;
                bin_imx[xi][yi] = NULL;
                bin_imy[xi][yi] = NULL;
            }
        }
    }
    
    /* */
    if(*count_napmc>0){
        *c_napmc = (int32_T*) malloc(sizeof(int32_T) * (size_t) *count_napmc);
        *l_napmc = (int32_T*) malloc(sizeof(int32_T) * (size_t) *count_napmc);
    } else {
        *c_napmc = NULL;
        *l_napmc = NULL;
    }
    *count_napmc = 0;
    
    for(xi=0;xi<S_im;xi++){
        for(yi=0;yi<L_im;yi++){
            bin_count_im[xi][yi] = 0;
        }
    }
    // printf("%d\n",msldemc_samples);
    for(l=0;l<msldemt_lines;l++){
        for(c=0;c<msldemt_samples;c++){
            if(msldemt_inImage[c][l]>1){
                pmcx = msldemt_x[l];
                pmcy = msldemt_y[c];
                pmcz = msldemt_z[c][l];
                //pmcx = ppvgx-cam_C[0]; pmcy = ppvgy-cam_C[1]; pmcz = ppvgz-cam_C[2];
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
                ++bin_count_im[xi][yi];

            } else if(msldemt_inImage[c][l]==1){
                (*c_napmc)[*count_napmc] = c;
                (*l_napmc)[*count_napmc] = l;
                (*count_napmc)++;
            }
        }
    }
}

/* main computation routine */
void mask_obstructed_pts_in_msldemt_using_msldemc(
        int32_T msldemc_samples, int32_T msldemc_lines,
        double *msldemc_x, double *msldemc_y, double **msldemc_z, int8_T **msldemc_imFOVmask, 
        int32_T msldemt_samples, int32_T msldemt_lines,
        double *msldemt_x, double *msldemt_y, double **msldemt_z, int8_T **msldemt_inImage,
        int32_T S_im, int32_T L_im, CAHV_MODEL cahv_mdl
        )
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
    
    int32_T **bin_count_im, *bin_count_im_base;
    int32_T ***bin_im_c, **bin_im_c_base;
    int32_T ***bin_im_l, **bin_im_l_base;
    int32_T xi,yi;
    int32_T x_min_int,x_max_int,y_min_int,y_max_int;
    int32_T n;
    int32_T S_imm1,L_imm1;
    int32_T *c_napmc, *l_napmc;
    int32_T count_napmc;
    
    double ***bin_imx, **bin_imx_base;
    double ***bin_imy, **bin_imy_base;
    
    double apmc,pmcx,pmcy,pmcz;
    double ppvgx,ppvgy,ppvgz,ppvx,ppvy;
    double *cam_A,*cam_H,*cam_V;
    
    if(msldemt_x==NULL)
        msldemt_x = msldemc_x;
    if(msldemt_y==NULL)
        msldemt_y = msldemc_y;
    if(msldemt_z==NULL)
        msldemt_z = msldemc_z;
    
    cam_A = cahv_mdl.A; cam_H = cahv_mdl.H; cam_V = cahv_mdl.V;
    
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
    createInt32Matrix(&bin_count_im,&bin_count_im_base,(size_t) S_im, (size_t) L_im);
    createInt32PMatrix(&bin_im_c, &bin_im_c_base, (size_t) S_im, (size_t) L_im);
    createInt32PMatrix(&bin_im_l, &bin_im_l_base, (size_t) S_im, (size_t) L_im);
    createDoublePMatrix(&bin_imx, &bin_imx_base, (size_t) S_im, (size_t) L_im);
    createDoublePMatrix(&bin_imy, &bin_imy_base, (size_t) S_im, (size_t) L_im);
    
    bin_msldemt_xyz_wAHVint(S_im, L_im, cam_A, cam_H, cam_V,
            msldemt_samples, msldemt_lines, msldemt_x, msldemt_y, msldemt_z, msldemt_inImage,
            bin_count_im, bin_im_c, bin_im_l, bin_imx, bin_imy,&count_napmc,&c_napmc,&l_napmc);
    
    
    
    /* Main Loop *********************************************************/
    
//     find_hidden_main_loop(msldemc_samples, msldemc_lines, msldemc_imFOVmask,
//         msldemc_northing, msldemc_easting, msldemc_img,
//         cam_A, cam_H, cam_V,
//         S_im, L_im, bin_count_im, bin_im_c, bin_im_l, bin_imx, bin_imy,
//         msldemc_northing, msldemc_easting, msldemc_img, msldemc_inImage);

               
    for(l=0;l<L_demcm1;l++){
        
        // decide the first and last indexes to be assessed.
        for(c=0;c<S_demcm1;c++){
            // process if 
            //printf("c=%d,mask_lc = %d,mask_lp1c = %d\n",c,msldemc_imFOVmask[c][l],msldemc_imFOVmask[c][l+1]);
            // printf("c=%d/%d\n",c,S_demcm1);
            if((msldemc_imFOVmask[c][l]>0) || (msldemc_imFOVmask[c][l+1]>0)){
                //printf("c=%d\n",c);
                for(ti=0;ti<2;ti++){
                    if(ti==0){
                        ppv1gx = msldemc_x[l];
                        ppv1gy = msldemc_y[c];
                        ppv1gz = msldemc_z[c][l];
                        ppv2gx = msldemc_x[l];
                        ppv2gy = msldemc_y[c+1];
                        ppv2gz = msldemc_z[c+1][l];
                        ppv3gx = msldemc_x[l+1];
                        ppv3gy = msldemc_y[c];
                        ppv3gz = msldemc_z[c][l+1];
                        isinFOVd = ((msldemc_imFOVmask[c][l]>1) && (msldemc_imFOVmask[c+1][l]>1) && (msldemc_imFOVmask[c][l+1]>1));
                        isinFOV = ((msldemc_imFOVmask[c][l]>0) && (msldemc_imFOVmask[c+1][l]>0) && (msldemc_imFOVmask[c][l+1]>0));
                        cv1 = c;   lv1 = l;
                        cv2 = c+1; lv2 = l;
                        cv3 = c;   lv3 = l+1;
                        // cv4 = c+1; lv4 = l+1;
                    }
                    else{
                        ppv1gx = msldemc_x[l];
                        ppv1gy = msldemc_y[c+1];
                        ppv1gz = msldemc_z[c+1][l];
                        ppv2gx = msldemc_x[l+1];
                        ppv2gy = msldemc_y[c+1];
                        ppv2gz = msldemc_z[c+1][l+1];
                        ppv3gx = msldemc_x[l+1];
                        ppv3gy = msldemc_y[c];
                        ppv3gz = msldemc_z[c][l+1];
                        cv1 = c+1; lv1 = l;
                        cv2 = c+1; lv2 = l+1;
                        cv3 = c;   lv3 = l+1;
                        // cv4 = c;   lv4 = l;
                        isinFOVd = ((msldemc_imFOVmask[c+1][l]>1) && (msldemc_imFOVmask[c+1][l+1]>1) && (msldemc_imFOVmask[c][l+1]>1));
                        isinFOV = ((msldemc_imFOVmask[c+1][l]>0) && (msldemc_imFOVmask[c+1][l+1]>0) && (msldemc_imFOVmask[c][l+1]>0));
                    }
                    
                    if(isinFOVd){
                        // ppv1gx = msldemc_x[lv1];
                        // ppv1gy = msldemc_y[cv1];
                        // ppv1gz = msldemc_z[cv1][lv1];
                        // ppv2gx = msldemc_x[lv2];
                        // ppv2gy = msldemc_y[cv2];
                        // ppv2gz = msldemc_z[cv2][lv2];
                        // ppv3gx = msldemc_x[lv3];
                        // ppv3gy = msldemc_y[cv3];
                        // ppv3gz = msldemc_z[cv3][lv3];
                        
                        /* Evaluate the projection */
                        // pmcx = ppv1gx; pmcy = ppv1gy; pmcz = ppv1gz;
                        apmc  =  ppv1gx*cam_A[0] + ppv1gy*cam_A[1] + ppv1gz*cam_A[2];
                        ppv1x = (ppv1gx*cam_H[0] + ppv1gy*cam_H[1] + ppv1gz*cam_H[2])/apmc;
                        ppv1y = (ppv1gx*cam_V[0] + ppv1gy*cam_V[1] + ppv1gz*cam_V[2])/apmc;
                        
                        // pmcx = ppv2gx; pmcy = ppv2gy; pmcz = ppv2gz;
                        apmc  =  ppv2gx*cam_A[0] + ppv2gy*cam_A[1] + ppv2gz*cam_A[2];
                        ppv2x = (ppv2gx*cam_H[0] + ppv2gy*cam_H[1] + ppv2gz*cam_H[2])/apmc;
                        ppv2y = (ppv2gx*cam_V[0] + ppv2gy*cam_V[1] + ppv2gz*cam_V[2])/apmc;
                        
                        // pmcx = ppv3gx; pmcy = ppv3gy; pmcz = ppv3gz;
                        apmc  =  ppv3gx*cam_A[0] + ppv3gy*cam_A[1] + ppv3gz*cam_A[2];
                        ppv3x = (ppv3gx*cam_H[0] + ppv3gy*cam_H[1] + ppv3gz*cam_H[2])/apmc;
                        ppv3y = (ppv3gx*cam_V[0] + ppv3gy*cam_V[1] + ppv3gz*cam_V[2])/apmc;
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
                            for(yi=y_min_int;yi<y_max_int;yi++){
                                for (n=0;n<bin_count_im[xi][yi];n++){
                                    cc = bin_im_c[xi][yi][n];
                                    ll = bin_im_l[xi][yi][n];
                                    /* evaluate line param */
                                    pmcx = msldemt_x[ll];
                                    pmcy = msldemt_y[cc];
                                    pmcz = msldemt_z[cc][ll];
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
                        lprm_nume = pnx*ppv1gx+pny*ppv1gy+pnz*ppv1gz;
                        
                        for(xi=0;xi<S_im;xi++){
                            for(yi=0;yi<L_im;yi++){
                                for (n=0;n<bin_count_im[xi][yi];n++){
                                    cc = bin_im_c[xi][yi][n];
                                    ll = bin_im_l[xi][yi][n];
                                    /* evaluate line param */
                                    pmcx = msldemt_x[ll];
                                    pmcy = msldemt_y[cc];
                                    pmcz = msldemt_z[cc][ll];
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
                                pmcx = msldemt_x[ll];
                                pmcy = msldemt_y[cc];
                                pmcz = msldemt_z[cc][ll];
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
    
    /* */
    //for(c=0;c<msldemt_samples;c++){
    //    for(l=0;l<msldemt_lines;l++){
    //        if(msldemt_inImage[c][l]<3)
    //            msldemt_inImage[c][l] = 0;
    //    }
    //}
    
    /* free dynamically allocated memories */
    for(xi=0;xi<S_im;xi++){
        for(yi=0;yi<L_im;yi++){
            if(bin_count_im[xi][yi]>0){
                free(bin_im_c[xi][yi]);
                free(bin_im_l[xi][yi]);
                free(bin_imx[xi][yi]);
                free(bin_imy[xi][yi]);
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
    
    free(bin_count_im);
    free(bin_count_im_base);
    
}

#endif