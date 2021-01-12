
#ifndef LIB_PROH_MASTCAMMSLDEM_L2_H
#define LIB_PROH_MASTCAMMSLDEM_L2_H

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
void mask_obstructed_pts_in_msldemt_using_msldemc_L2(char *msldem_imgpath, EnviHeader msldem_header, 
        int32_T msldemc_imxy_sample_offset, int32_T msldemc_imxy_line_offset,
        int32_T msldemc_samples, int32_T msldemc_lines,
        double *msldemc_xmc, double *msldemc_ymc, int8_T **msldemc_imFOVmask, 
        double *msldemt_xmc, double *msldemt_ymc, int8_T **msldemt_inImage,
        int32_T **bin_count_im, int32_T ***bin_im_c, int32_T ***bin_im_l,
        double ***bin_imx, double ***bin_imy, double ***bin_demzmc,
        int32_T count_napmc, int32_T *c_napmc, int32_T *l_napmc, double *zmc_napmc,
        int32_T S_im, int32_T L_im, CAHV_MODEL cahv_mdl)
{
    int32_T c,l,cc,ll;
    int32_T cv1,cv2,cv3,lv1,lv2,lv3; /* v indicates vertex */
    int32_T L_demcm1,S_demcm1;
    int16_T ti; /* triangle index */
    long skip_pri;
    long skip_l, skip_r;
    float *elevl,*elevlp1;
    long ncpy;
    const int sz=sizeof(float);
    FILE *fid;
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
    
    if(msldemt_xmc==NULL)
        msldemt_xmc = msldemc_xmc;
    if(msldemt_ymc==NULL)
        msldemt_ymc = msldemc_ymc;
    
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
    
    fid = fopen(msldem_imgpath,"rb");
    
    /* skip lines */
    skip_pri = (long) msldem_header.samples * (long) msldemc_imxy_line_offset * (long) sz;
    // printf("%d*%d*%d=%ld\n",msldem_header.samples,msldemc_imxy_line_offset,s,skip_pri);
    fseek(fid,skip_pri,SEEK_CUR);
    
    elevl = (float*) malloc(sz*msldemc_samples);
    elevlp1 = (float*) malloc(sz*msldemc_samples);
    skip_l = (long) sz * (long) msldemc_imxy_sample_offset;
    skip_r = ((long) msldem_header.samples - (long) msldemc_samples)* (long) sz - skip_l;
    ncpy = (long) msldemc_samples* (long)sz;
    fseek(fid,skip_l,SEEK_CUR);
    fread(elevlp1,sz,msldemc_samples,fid);
    fseek(fid,skip_r,SEEK_CUR);
    
               
    for(l=0;l<L_demcm1;l++){
        memcpy(elevl,elevlp1,ncpy);
        fseek(fid,skip_l,SEEK_CUR);
        fread(elevlp1,sz,msldemc_samples,fid);
        fseek(fid,skip_r,SEEK_CUR);
        // decide the first and last indexes to be assessed.
        for(c=0;c<S_demcm1;c++){
            // process if 
            //printf("c=%d,mask_lc = %d,mask_lp1c = %d\n",c,msldemc_imFOVmask[c][l],msldemc_imFOVmask[c][l+1]);
            // printf("c=%d/%d\n",c,S_demcm1);
            if((msldemc_imFOVmask[c][l]>0) || (msldemc_imFOVmask[c][l+1]>0)){
                //printf("c=%d\n",c);
                for(ti=0;ti<2;ti++){
                    if(ti==0){
                        ppv1gx = msldemc_xmc[l];
                        ppv1gy = msldemc_ymc[c];
                        ppv1gz = ((double) -elevl[c]) - cam_C[2];
                        ppv2gx = msldemc_xmc[l];
                        ppv2gy = msldemc_ymc[c+1];
                        ppv2gz = ((double) -elevl[c+1]) - cam_C[2];
                        ppv3gx = msldemc_xmc[l+1];
                        ppv3gy = msldemc_ymc[c];
                        ppv3gz = ((double) -elevlp1[c]) - cam_C[2];
//                         ppv1gx = msldemc_x[l];
//                         ppv1gy = msldemc_y[c];
//                         ppv1gz = msldemc_z[c][l];
//                         ppv2gx = msldemc_x[l];
//                         ppv2gy = msldemc_y[c+1];
//                         ppv2gz = msldemc_z[c+1][l];
//                         ppv3gx = msldemc_x[l+1];
//                         ppv3gy = msldemc_y[c];
//                         ppv3gz = msldemc_z[c][l+1];
                        isinFOVd = ((msldemc_imFOVmask[c][l]>1) && (msldemc_imFOVmask[c+1][l]>1) && (msldemc_imFOVmask[c][l+1]>1));
                        isinFOV = ((msldemc_imFOVmask[c][l]>0) && (msldemc_imFOVmask[c+1][l]>0) && (msldemc_imFOVmask[c][l+1]>0));
                        cv1 = c;   lv1 = l;
                        cv2 = c+1; lv2 = l;
                        cv3 = c;   lv3 = l+1;
                        // cv4 = c+1; lv4 = l+1;
                    }
                    else{
                        ppv1gx = msldemc_xmc[l];
                        ppv1gy = msldemc_ymc[c+1];
                        ppv1gz = ((double) -elevl[c+1]) - cam_C[2];
                        ppv2gx = msldemc_xmc[l+1];
                        ppv2gy = msldemc_ymc[c+1];
                        ppv2gz = ((double) -elevlp1[c+1]) - cam_C[2];
                        ppv3gx = msldemc_xmc[l+1];
                        ppv3gy = msldemc_ymc[c];
                        ppv3gz = ((double) -elevlp1[c]) - cam_C[2];
//                         ppv1gx = msldemc_x[l];
//                         ppv1gy = msldemc_y[c+1];
//                         ppv1gz = msldemc_z[c+1][l];
//                         ppv2gx = msldemc_x[l+1];
//                         ppv2gy = msldemc_y[c+1];
//                         ppv2gz = msldemc_z[c+1][l+1];
//                         ppv3gx = msldemc_x[l+1];
//                         ppv3gy = msldemc_y[c];
//                         ppv3gz = msldemc_z[c][l+1];
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
                                    pmcx = msldemt_xmc[ll];
                                    pmcy = msldemt_ymc[cc];
                                    // pmcz = msldemt_z[cc][ll];
                                    pmcz = bin_demzmc[xi][yi][n];
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
                                    pmcx = msldemt_xmc[ll];
                                    pmcy = msldemt_ymc[cc];
                                    // pmcz = msldemt_z[cc][ll];
                                    pmcz = bin_demzmc[xi][yi][n];
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
                                pmcx = msldemt_xmc[ll];
                                pmcy = msldemt_ymc[cc];
                                // pmcz = msldemt_z[cc][ll];
                                pmcz = zmc_napmc[n];
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
    free(elevl);
    free(elevlp1);
    fclose(fid);
    
}



#endif