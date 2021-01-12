/* =====================================================================
 * cahvor_get_imxy_MSLDEM_mex.c
 * Read MSLDEM image data
 * Perform projection of mastcam pixels onto MSLDEM data
 * 
 * INPUTS:
 * 0 msldem_imgpath        char*
 * 1 msldem_header         struct
 * 2 msldemc_header        struct
 * 3 msldemc_northing      Double array [L_demc]
 * 4 msldemc_easting       Double array [S_demc]
 * 5 msldemc_imFOVmaskd    int8 [L_demc x S_demc]
 * 6 cammdl               CAHVOR model class
 * 
 * 
 * OUTPUTS:
 * 0 msldemc_imx           Doubles [L_demc x S_demc]
 * 1 msldemc_imy           Doubles [L_demc x S_demc]
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
#include "cahvor.h"
#include "mex_create_array.h"

/* main computation routine */
void get_imxy_MSLDEM(char *msldem_imgpath, EnviHeader msldem_hdr, 
        int32_T msldemc_imxy_sample_offset, int32_T msldemc_imxy_line_offset,
        int32_T msldemc_samples, int32_T msldemc_lines,
        double *msldemc_northing, double *msldemc_easting,
        int8_T **msldemc_imFOVmaskd, 
        CAHVOR_MODEL cahvor_mdl,
        double **msldemc_imx, double **msldemc_imy)
{
    int32_T c,l;
    int16_T sxy;
    long skip_pri;
    long skip_l, skip_r;
    float *elevl;
    size_t ncpy;
    size_t sz=sizeof(float);
    FILE *fid;
    double dem_cl;
    double pmcx,pmcy,pmcz;
    double pmc_o_mag;
    double pmcnox,pmcnoy,pmcnoz;
    double lamx,lamy,lamz;
    double lam_mag2;
    double mup1;
    double pdmcx,pdmcy,pdmcz;
    double apdmc,hpdmc,vpdmc;
    double x_im,y_im;
    double *PmCys;
    double *cam_C, *cam_A, *cam_H, *cam_V, *cam_O, *cam_R;
    
    cam_C = cahvor_mdl.C; cam_A = cahvor_mdl.A;
    cam_H = cahvor_mdl.H; cam_V = cahvor_mdl.V;
    cam_O = cahvor_mdl.O; cam_R = cahvor_mdl.R;
    
    PmCys = (double*) malloc((size_t) msldemc_samples * sizeof(double));
    for(c=0;c<msldemc_samples;c++){
        pmcy  = msldemc_easting[c] - cam_C[1];
        PmCys[c] = pmcy;
    }

    fid = fopen(msldem_imgpath,"rb");
    
    /* skip lines */
    skip_pri = (long) msldem_hdr.samples * (long) msldemc_imxy_line_offset * (long) sz;
    // printf("%d*%d*%d=%ld\n",msldem_header.samples,msldemc_imxy_line_offset,s,skip_pri);
    fseek(fid,skip_pri,SEEK_CUR);
    
    /* read the data */
    ncpy = sz * (size_t) msldemc_samples;
    elevl = (float*) malloc(ncpy);
    skip_l = (long) sz * (long) msldemc_imxy_sample_offset;
    skip_r = ((long) msldem_hdr.samples - (long) msldemc_samples)* (long) sz - skip_l;
    
    
    // printf("%d,%d,%d\n",skip_l,msldemc_samples*s,skip_r);
    for(l=0;l<msldemc_lines;l++){
        fseek(fid,skip_l,SEEK_CUR);
        fread(elevl,sz,msldemc_samples,fid);
        fseek(fid,skip_r,SEEK_CUR);
        pmcx  = msldemc_northing[l] - cam_C[0];
        for(c=0;c<msldemc_samples;c++){
            if(msldemc_imFOVmaskd[c][l]>1){
                dem_cl = (double) elevl[c];
                pmcz  = -dem_cl-cam_C[2];
                pmcy   = PmCys[c];
                
                pmc_o_mag = pmcx*cam_O[0]+pmcy*cam_O[1]+pmcz*cam_O[2];
                pmcnox = pmcx / pmc_o_mag;
                pmcnoy = pmcy / pmc_o_mag;
                pmcnoz = pmcz / pmc_o_mag;
                /* then get the vector "lambda" */
                lamx = pmcnox - cam_O[0];
                lamy = pmcnoy - cam_O[1];
                lamz = pmcnoz - cam_O[2];
                /* next get the parameter "mu" */
                lam_mag2 = lamx*lamx + lamy*lamy + lamz*lamz;
                mup1 = 1+cam_R[0] + cam_R[1]*lam_mag2 + cam_R[2]*lam_mag2*lam_mag2;
                pdmcx = cam_O[0] + mup1*lamx;
                pdmcy = cam_O[1] + mup1*lamy;
                pdmcz = cam_O[2] + mup1*lamz;

                /* perform projection */
                apdmc = cam_A[0] * pdmcx + cam_A[1] * pdmcy + cam_A[2] * pdmcz;
                hpdmc = cam_H[0] * pdmcx + cam_H[1] * pdmcy + cam_H[2] * pdmcz;
                vpdmc = cam_V[0] * pdmcx + cam_V[1] * pdmcy + cam_V[2] * pdmcz;
                x_im = hpdmc/apdmc; y_im = vpdmc/apdmc;
                
                msldemc_imx[c][l] = x_im;
                msldemc_imy[c][l] = y_im;
            }

        }

    }
    free(elevl);
    free(PmCys);
    fclose(fid);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *msldem_imgpath;
    EnviHeader msldem_hdr;
    CAHVOR_MODEL cahvor_mdl;
    mwSize msldemc_imxy_sample_offset,msldemc_imxy_line_offset;
    double *msldemc_northing;
    double *msldemc_easting;
    int8_T **msldemc_imFOVmaskd;
    mwSize S_im,L_im;
    
    double **msldemc_imx;
    double **msldemc_imy;
    

    
    mwSize si,li;
    mwSize msldemc_samples, msldemc_lines;

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
    msldem_hdr = mxGetEnviHeader(prhs[1]);
    
    /* INPUT 2 msldemc_sheader*/
    msldemc_imxy_sample_offset = (mwSize) mxGetScalar(mxGetField(prhs[2],0,"sample_offset"));
    msldemc_imxy_line_offset = (mwSize) mxGetScalar(mxGetField(prhs[2],0,"line_offset"));
    msldemc_samples = (mwSize) mxGetScalar(mxGetField(prhs[2],0,"samples"));
    msldemc_lines = (mwSize) mxGetScalar(mxGetField(prhs[2],0,"lines"));
    
    /* INPUT 3/4 msldem northing easting */
    msldemc_northing = mxGetDoubles(prhs[3]);
    msldemc_easting = mxGetDoubles(prhs[4]);
    

    /* INPUT 5 msldem imFOV */
    msldemc_imFOVmaskd = set_mxInt8Matrix(prhs[5]);
    
    /* INPUT 6 camera model */
    cahvor_mdl = mxGet_CAHVOR_MODEL(prhs[6]);
    
    /* OUTPUT 0/1 msldem imxy */
    plhs[0] = mxCreateDoubleMatrix(msldemc_lines,msldemc_samples,mxREAL);
    msldemc_imx = set_mxDoubleMatrix(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(msldemc_lines,msldemc_samples,mxREAL);
    msldemc_imy = set_mxDoubleMatrix(plhs[1]);  
    

    
    // Initialize matrices
    for(si=0;si<msldemc_samples;si++){
        for(li=0;li<msldemc_lines;li++){
            msldemc_imx[si][li] = NAN;
            msldemc_imy[si][li] = NAN;            
        }
    }
    
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    get_imxy_MSLDEM(msldem_imgpath, msldem_hdr, 
        (int32_T) msldemc_imxy_sample_offset, (int32_T) msldemc_imxy_line_offset,
        (int32_T) msldemc_samples, (int32_T) msldemc_lines,
        msldemc_northing, msldemc_easting,
        msldemc_imFOVmaskd, 
        cahvor_mdl,
        msldemc_imx, msldemc_imy);
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldemc_imx);
    mxFree(msldemc_imy);
    mxFree(msldemc_imFOVmaskd);
    
}
