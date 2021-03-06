/* =====================================================================
 * get_imFOVmask_MSLDEM_v2_mex.c
 * Evaluate if pixels in the MSL DEM image are potentially in 
 * 
 * INPUTS:
 * 0 msldem_imgpath        char*
 * 1 msldem_header         struct
 * 2 msldem_northing      Double array [L_dem]
 * 3 msldem_easting       Double array [S_dem]
 * 4 S_im                 int
 * 5 L_im                 int
 * 6 cammdl               CAHVOR model class
 * 7 coef_mrgn            coefficient for the margin
 * 
 * 
 * OUTPUTS:
 * 0 msldemc_imFOVmaskd    int8 [L_dem x S_dem]
 * // 1 msldemc_imFOVmaskd    Boolean [L_dem x S_dem]
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
void get_imFOVmask_MSLDEM(char *msldem_imgpath, EnviHeader msldem_hdr,
        double *msldem_northing, double *msldem_easting, 
        int32_T S_im, int32_T L_im, 
        CAHV_MODEL cahv_mdl,
        int8_T **msldem_imFOVmaskd, double coef_mrgn)
{
    int32_T c,l;
    
    float *elevlm1, *elevl ,*elevlp1;
    size_t nbytes_l;
    int32_T S_dem, L_dem;
    int32_T L_demm1;
    int32_T S_demp1;
    size_t sz = sizeof(float);
    //size_t sz1 = sizeof(size_t);
    FILE *fid;
    double dem_cl;
    double pmcx,pmcy,pmcz;
    double apmcx,apmcy,apmcz;
    double hpmcx,hpmcy,hpmcz;
    double vpmcx,vpmcy,vpmcz;
    double *APmCys,*HPmCys,*VPmCys;
    double apmc,hpmc,vpmc;
    double x_im,y_im;
    double resolhxy,resolvxy;
    double resolz;
    float elev_amax,elev_amin;
    float resol_amax, resol_amin;
    double mrgnh,mrgnv;
    double S_imm05,L_imm05;
    float data_ignore_value_float;
    double Hd2_abs,Vd2_abs;
    double *cam_C, *cam_A, *cam_H, *cam_V, *cam_Hd, *cam_Vd;
    double hc,vc,hs,vs;
    
    cam_C = cahv_mdl.C; cam_A = cahv_mdl.A; cam_H = cahv_mdl.H; cam_V = cahv_mdl.V;
    hs = cahv_mdl.hs; vs = cahv_mdl.vs; hc = cahv_mdl.hc; vc = cahv_mdl.vc;
    cam_Hd = cahv_mdl.Hdash; cam_Vd = cahv_mdl.Vdash;
    
    S_dem = (int32_T) msldem_hdr.samples;
    L_dem = (int32_T) msldem_hdr.lines;
    L_demm1 = L_dem - 1;
    S_demp1 = S_dem + 1;
    
    // printf("%d \n",L_dem);
    
    Hd2_abs = fabs(cam_Hd[2]);
    Vd2_abs = fabs(cam_Vd[2]);
    
    // printf("%d \n",L_dem);
    
    data_ignore_value_float = (float) msldem_hdr.data_ignore_value + 1.0;
    // printf("ignore: %f\n",data_ignore_value_float);
    
    // printf("%d \n",L_dem);
    
    S_imm05 = (double) S_im - 0.5;
    L_imm05 = (double) L_im - 0.5;
    
    // printf("size of size_t: %d\n",sz1);
    
    APmCys = (double*) malloc((size_t) S_dem * sizeof(double));
    HPmCys = (double*) malloc((size_t) S_dem * sizeof(double));
    VPmCys = (double*) malloc((size_t) S_dem * sizeof(double));
    for(c=0;c<S_dem;c++){
        pmcy  = msldem_easting[c] - cam_C[1];
        APmCys[c] = cam_A[1] * pmcy;
        HPmCys[c] = cam_H[1] * pmcy;
        VPmCys[c] = cam_V[1] * pmcy;
    }
    
    resolhxy = fabs(cam_Hd[0]) + fabs(cam_Hd[1]);
    resolvxy = fabs(cam_Vd[0]) + fabs(cam_Vd[1]);
    
    /* read the data */
    nbytes_l = sz * ((size_t) S_dem+2);
    elevlm1  = (float*) malloc(nbytes_l);
    elevl    = (float*) malloc(nbytes_l);
    elevlp1  = (float*) malloc(nbytes_l);
    
    // printf("%d \n",L_dem);
    
    for(c=0;c<S_dem+2;c++){
        elevlm1[c] = NAN;
        elevl[c]   = NAN;
        elevlp1[c] = NAN;
    }
    
    // printf("%d \n",L_dem);
    
    fid = fopen(msldem_imgpath,"rb");
    
    fread(&elevlp1[1],sz,S_dem,fid);
    for(c=1;c<S_demp1;c++){
        if(elevlp1[c]<data_ignore_value_float)
            elevlp1[c] = NAN;
    }
    
    // printf("%d \n",L_dem);
    
    
    // printf("%d,%d,%d\n",skip_l,msldemc_samples*s,skip_r);
    for(l=0;l<L_dem;l++){
        // printf("l=%d \n",l);
        memcpy(elevlm1,elevl,nbytes_l);
        memcpy(elevl,elevlp1,nbytes_l);
        if(l<L_demm1){
            fread(&elevlp1[1],sz,S_dem,fid);
            for(c=1;c<S_demp1;c++){
                if(elevlp1[c]<data_ignore_value_float)
                    elevlp1[c] = NAN;
            }
        } else {
            for(c=1;c<S_demp1;c++){
                elevlp1[c] = NAN;
            }
        }
        
        pmcx  = msldem_northing[l] - cam_C[0];
        apmcx = cam_A[0] * pmcx;
        hpmcx = cam_H[0] * pmcx;
        vpmcx = cam_V[0] * pmcx;
        for(c=0;c<S_dem;c++){
            dem_cl = (double) elevl[c+1];
            pmcz  = -dem_cl-cam_C[2];
            apmcz = cam_A[2] * pmcz;
            apmcy = APmCys[c];
            apmc = apmcx + apmcy + apmcz;
            //apmc = apmcx + APmCys[c] + cam_A[2] * pmcz;
            if(apmc>0){
                hpmcy = HPmCys[c];
                vpmcy = VPmCys[c];
                hpmcz = cam_H[2] * pmcz;
                vpmcz = cam_V[2] * pmcz;
                hpmc = hpmcx + hpmcy + hpmcz;
                vpmc = vpmcx + vpmcy + vpmcz;
                x_im = hpmc / apmc;
                y_im = vpmc / apmc;
                //x_im = (hpmcx + HPmCys[c] + cam_H[2] * pmcz) / apmc;
                //y_im = (vpmcx + VPmCys[c] + cam_V[2] * pmcz) / apmc;
                
                /* Evaluate resolution */
                elev_amax = elevl[c+1];
                elev_amin = elevl[c+1];
                
                if (elev_amax<elevlm1[c])
                    elev_amax = elevlm1[c];
                else if(elev_amin>elevlm1[c])
                    elev_amin = elevlm1[c];
                
                if (elev_amax<elevlm1[c+1])
                    elev_amax = elevlm1[c+1];
                else if(elev_amin>elevlm1[c+1])
                    elev_amin = elevlm1[c+1];
                
                if (elev_amax<elevlm1[c+2])
                    elev_amax = elevlm1[c+2];
                else if(elev_amin>elevlm1[c+2])
                    elev_amin = elevlm1[c+2];
                
                if (elev_amax<elevl[c])
                    elev_amax = elevl[c];
                else if(elev_amin>elevl[c])
                    elev_amin = elevl[c];
                
                if (elev_amax<elevl[c+2])
                    elev_amax = elevl[c+2];
                else if(elev_amin>elevl[c+2])
                    elev_amin = elevl[c+2];
                
                if (elev_amax<elevlp1[c])
                    elev_amax = elevlp1[c];
                else if(elev_amin>elevlp1[c])
                    elev_amin = elevlp1[c];
                
                if (elev_amax<elevlp1[c+1])
                    elev_amax = elevlp1[c+1];
                else if(elev_amin>elevlp1[c+1])
                    elev_amin = elevlp1[c+1];
                
                if (elev_amax<elevlp1[c+2])
                    elev_amax = elevlp1[c+2];
                else if(elev_amin>elevlp1[c+2])
                    elev_amin = elevlp1[c+2];
                
                resol_amax = fabsf(elev_amax-elevl[c+1]);
                resol_amin = fabsf(elev_amin-elevl[c+1]);
                
                if(resol_amax>resol_amin)
                    resolz = (double) resol_amax;
                else
                    resolz = (double) resol_amin;
                 
                mrgnh = coef_mrgn*hs/apmc * (resolhxy + Hd2_abs * resolz);
                mrgnv = coef_mrgn*vs/apmc * (resolvxy + Vd2_abs * resolz);
                
                if (x_im>-0.5-mrgnh && x_im<S_imm05+mrgnh && 
                        y_im>-0.5-mrgnv && y_im<L_imm05+mrgnv){
                    msldem_imFOVmaskd[c][l] = 4;
                }                
            } else {
                /* Evaluate */
                
                
            }

        }

    }
    free(elevlm1);
    free(elevl);
    free(elevlp1);
    free(APmCys);
    free(HPmCys);
    free(VPmCys);
    fclose(fid);
    
    /* safeguarding not implemented yet */
    
    
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *msldem_imgpath;
    EnviHeader msldem_hdr;
    CAHV_MODEL cahv_mdl;
    double *msldem_northing;
    double *msldem_easting;
    // bool **msldem_imFOVmask;
    int8_T **msldem_imFOVmaskd;
    mwSize S_im,L_im;
    
    double coef_mrgn;
    
    int32_T si,li;
    int32_T msldem_samples, msldem_lines;

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
    
    /* INPUT 2/3 msldem northing easting */
    msldem_northing = mxGetDoubles(prhs[2]);
    msldem_easting = mxGetDoubles(prhs[3]);

    /* INPUT 4/5 image S_im, L_im */
    S_im = (mwSize) mxGetScalar(prhs[4]);
    L_im = (mwSize) mxGetScalar(prhs[5]);
    //printf("sim = %d\n",S_im);
    
    /* INPUT 6 camera model */
    cahv_mdl = mxGet_CAHV_MODEL(prhs[6]);
    
    // printf("cam_Hd %f %f %f\n",cam_Hd[0],cam_Hd[1],cam_Hd[2]);
    coef_mrgn = mxGetScalar(prhs[7]);
    
    /* OUTPUT 0 msldem imFOV */
    // plhs[0] = mxCreateLogicalMatrix((mwSize) msldem_hdr.lines,(mwSize) msldem_hdr.samples);
    // msldem_imFOVmask = set_mxLogicalMatrix(plhs[0]);

    /* OUTPUT 1 msldem imFOV */
    plhs[0] = mxCreateNumericMatrix((mwSize) msldem_hdr.lines,(mwSize) msldem_hdr.samples,mxINT8_CLASS,mxREAL);
    msldem_imFOVmaskd = set_mxInt8Matrix(plhs[0]);
    
    // Initialize matrices
    msldem_samples = (int32_T) msldem_hdr.samples;
    msldem_lines = (int32_T) msldem_hdr.lines;
    for(si=0;si<msldem_samples;si++){
        for(li=0;li<msldem_lines;li++){
            // msldem_imFOVmask[si][li] = false;
            msldem_imFOVmaskd[si][li] = 0;
        }
    }
    // printf("sim = %d\n",S_im);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    get_imFOVmask_MSLDEM(msldem_imgpath, msldem_hdr,
        msldem_northing, msldem_easting, 
        (int32_T) S_im, (int32_T) L_im, 
        cahv_mdl,
        msldem_imFOVmaskd,coef_mrgn);
    
    /* free memories */
    mxFree(msldem_imgpath);
    // mxFree(msldem_imFOVmask);
    mxFree(msldem_imFOVmaskd);
    
}
