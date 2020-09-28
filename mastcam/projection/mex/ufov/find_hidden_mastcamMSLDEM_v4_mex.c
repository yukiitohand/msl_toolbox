/* =====================================================================
 * find_hidden_mastcamMSLDEM_v4_mex.c
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
        int32_T S_im, int32_T L_im,
        bool **msldemc_inImage)
{
    int32_T c,l,cc,ll;
    int32_T cs1,cs2,cs3,cs4,ls1,ls2,ls3,ls4;
    int32_T L_demcm1,S_demcm1;
    int16_T sxy;
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
    
    printf("%d\n",msldemc_samples);
    L_demcm1 = msldemc_lines-1;
    S_demcm1 = msldemc_samples-1;
    printf("%d\n",msldemc_samples);
    
                
    for(l=0;l<1;l++){
        printf("l=%d\n",l);
        // decide the first and last indexes to be assessed.
        for(c=0;c<S_demcm1;c++){
            // process if 
            //printf("c=%d,mask_lc = %d,mask_lp1c = %d\n",c,msldemc_imFOVmask[c][l],msldemc_imFOVmask[c][l+1]);
            // printf("c=%d/%d\n",c,S_demcm1);
            if(msldemc_imFOVmask[c][l] || msldemc_imFOVmask[c][l+1]){
                //printf("c=%d\n",c);
                for(sxy=0;sxy<2;sxy++){
                    if(sxy==0){
                        cs1 = c; ls1 = l;
                        cs2 = c+1; ls2 = l;
                        cs3 = c; ls3 = l+1;
                        cs4 = c+1; ls4 = l+1;
                    }
                    else{
                        cs1 = c+1; ls1 = l;
                        cs2 = c+1; ls2 = l+1;
                        cs3 = c; ls3 = l+1;
                        cs4 = c; ls4 = l;
                    }
                    // isinFOVd = (msldemc_imFOVmaskd[cs1][ls1] && msldemc_imFOVmaskd[cs2][ls2] && msldemc_imFOVmaskd[cs3][ls3]);
                    // isinFOV = (msldemc_imFOVmask[cs1][ls1] && msldemc_imFOVmask[cs2][ls2] && msldemc_imFOVmask[cs3][ls3]);

                    //printf("sxy=%d\n",sxy);
                    //printf("isinFOVd=%d\n",isinFOVd);
                    //printf("%d\n",c);
                    if(msldemc_imFOVmaskd[cs1][ls1] && msldemc_imFOVmaskd[cs2][ls2] && msldemc_imFOVmaskd[cs3][ls3])
                    {
                        ppv1x = msldemc_imx[cs1][ls1];
                        ppv1y = msldemc_imy[cs1][ls1];
                        ppv2x = msldemc_imx[cs2][ls2];
                        ppv2y = msldemc_imy[cs2][ls2];
                        ppv3x = msldemc_imx[cs3][ls3];
                        ppv3y = msldemc_imy[cs3][ls3];
                        ppv1gx = msldemc_northing[ls1];
                        ppv1gy = msldemc_easting[cs1];
                        ppv1gz = msldemc_img[cs1][ls1];
                        ppv2gx = msldemc_northing[ls2];
                        ppv2gy = msldemc_easting[cs2];
                        ppv2gz = msldemc_img[cs2][ls2];
                        ppv3gx = msldemc_northing[ls3];
                        ppv3gy = msldemc_easting[cs3];
                        ppv3gz = msldemc_img[cs3][ls3];
                        ppv4gx = msldemc_northing[ls4];
                        ppv4gy = msldemc_easting[cs4];
                        ppv4gz = msldemc_img[cs4][ls4];

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
                        
                        // pre-screening
                        x_min = fmin(fmin(ppv1x,ppv2x),ppv3x);
                        y_min = fmin(fmin(ppv1y,ppv2y),ppv3y);
                        x_max = fmax(fmax(ppv1x,ppv2x),ppv3x);
                        y_max = fmax(fmax(ppv1y,ppv2y),ppv3y);

                       for(ll=0;ll<msldemc_lines;ll++){
                            //printf("ll=%d/%d\n",ll,msldemc_lines);
                            for(cc=0;cc<msldemc_samples;cc++){
                                if(msldemc_inImage[cc][ll]){
                                    /* If  */
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
                                                msldemc_inImage[cc][ll] = false;
                                            }
                                        }
                                    }
                               }
                           }
                        }
                    }
                    else if(msldemc_imFOVmask[cs1][ls1] && msldemc_imFOVmask[cs2][ls2] && msldemc_imFOVmask[cs3][ls3]) // if(isinFOV) This mode is rarely invoked and not debugged well.
                    {
                        /* parameters for plane equations
                         * plane normal vector (pn)
                         * plane constant (pc)
                        */
                        printf("t");
                        ppv1gx = msldemc_northing[ls1];
                        ppv1gy = msldemc_easting[cs1];
                        ppv1gz = -msldemc_img[cs1][ls1];
                        ppv2gx = msldemc_northing[ls2];
                        ppv2gy = msldemc_easting[cs2];
                        ppv2gz = -msldemc_img[cs2][ls2];
                        ppv3gx = msldemc_northing[ls3];
                        ppv3gy = msldemc_easting[cs3];
                        ppv3gz = -msldemc_img[cs3][ls3];
                        ppv4gx = msldemc_northing[ls4];
                        ppv4gy = msldemc_easting[cs4];
                        ppv4gz = -msldemc_img[cs4][ls4];

                        pdv1gx = ppv2gx - ppv1gx;
                        pdv1gy = ppv2gy - ppv1gy;
                        pdv1gz = ppv2gz - ppv1gz;
                        pdv2gx = ppv3gx - ppv1gx;
                        pdv2gy = ppv3gy - ppv1gy;
                        pdv2gz = ppv3gz - ppv1gz;
                        pnx = pdv1gy*pdv2gz - pdv1gz*pdv2gy;
                        pny = pdv1gz*pdv2gx - pdv1gx*pdv2gz;
                        pnz = pdv1gx*pdv2gy - pdv1gy*pdv2gx;
                        //pc = pnx*ppv1gx + pny*ppv1gy + pnz*ppv1gz;

                        lprm = (pnx*ppv1gx+pny*ppv1gy+pnz*ppv1gz)
                                       /(pnx*tppvgx+pny*tppvgy+pnz*tppvgz);

                        if(lprm<1 && lprm>0){
                            /* Get Plane parameters */
                            M[0][0] = pow(pdv1gx,2) + pow(pdv1gy,2) + pow(pdv1gz,2);
                            M[0][1] = pdv1gx*pdv2gx + pdv1gy*pdv2gy + pdv1gz*pdv2gz;
                            M[1][0] = M[0][1];
                            M[1][1] = pow(pdv2gx,2) + pow(pdv2gy,2) + pow(pdv2gz,2);
                            detM = M[0][0]*M[1][1] - pow(M[0][1],2);
                            Minv[0][0] = M[1][1]/detM;
                            Minv[0][1] = -M[0][1]/detM;
                            Minv[1][0] = -M[1][0]/detM;
                            Minv[1][1] = M[0][0]/detM;
                            Minvp[0][0] = Minv[0][0]*pdv1gx+Minv[0][1]*pdv2gx;
                            Minvp[0][1] = Minv[0][0]*pdv1gy+Minv[0][1]*pdv2gy;
                            Minvp[0][2] = Minv[0][0]*pdv1gz+Minv[0][1]*pdv2gz;
                            Minvp[1][0] = Minv[1][0]*pdv1gx+Minv[1][1]*pdv2gx;
                            Minvp[1][1] = Minv[1][0]*pdv1gy+Minv[1][1]*pdv2gy;
                            Minvp[1][2] = Minv[1][0]*pdv1gz+Minv[1][1]*pdv2gz;

                            // plane intersection pointing vector
                            pipvgx = lprm*tppvgx;
                            pipvgy = lprm*tppvgy;
                            pipvgz = lprm*tppvgz;

                            pipvgppv1x = pipvgx - ppv1gx;
                            pipvgppv1y = pipvgy - ppv1gy;
                            pipvgppv1z = pipvgz - ppv1gz;

                            // Get plane coefficiets
                            pprm_s = Minvp[0][0]*pipvgppv1x+Minvp[0][1]*pipvgppv1y+Minvp[0][2]*pipvgppv1z;
                            pprm_t = Minvp[1][0]*pipvgppv1x+Minvp[1][1]*pipvgppv1y+Minvp[1][2]*pipvgppv1z;
                            pprm_1st = 1 - pprm_s - pprm_t;

                            if(pprm_s>0 && pprm_t>0 && pprm_1st>0){
                                msldemc_inImage[cc][ll] = false;
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
    double *cam_C;
    double *cam_A;
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
    
    
    /* OUTPUT 0/1/2 north-east-elevation */
    plhs[0] = mxCreateLogicalMatrix(L_demc,S_demc);
    msldemc_inImage = set_mxLogicalMatrix(plhs[0]);
    
    // Initialize matrices
    for(si=0;si<S_demc;si++){
        for(li=0;li<L_demc;li++){
            msldemc_inImage[si][li] = msldemc_imFOVmask[si][li];
        }
    }
    printf("%d\n",S_im);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    find_hidden((int32_T) S_demc, (int32_T) L_demc, msldemc_img,
       msldemc_northing, msldemc_easting,
       msldemc_imx, msldemc_imy,
       msldemc_imFOVmask,msldemc_imFOVmaskd, 
       (int32_T) S_im, (int32_T) L_im,
       msldemc_inImage);
    
    /* free memories */
    mxFree(msldemc_img);
    mxFree(msldemc_imx);
    mxFree(msldemc_imy);
    mxFree(msldemc_imFOVmask);
    mxFree(msldemc_imFOVmaskd);
    mxFree(msldemc_inImage);
    
}
