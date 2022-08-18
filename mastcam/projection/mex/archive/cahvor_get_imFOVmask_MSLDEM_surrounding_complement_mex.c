/* =====================================================================
 * cahvor_get_imFOVmask_MSLDEM_surrounding_complement_mex.c
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
 *    5: inside the FOV (-0.5<x<S_im-0.5 && -0.5<y<L_im-0.5)
 *    4: tightest margin of the FOV (any of 8 neighbors is marked as 5)
 *    3: potential FOV at the short distance where pixel resolution is 
 *       comparable to the camera full FOV.
 *    2: supplemental marginal region (used later for the pointing vector 
 *       intersection.)
 *    1: not evaluated here (FOV with apmc<0)
 *    0: outside FOV
 * 
 * 2 is marked for 
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
        int32_T S_im, int32_T L_im, double *srange, double *lrange, 
        CAHVOR_MODEL cahvor_mdl,
        int8_T **msldem_imFOVmaskd, double coef_mrgn)
{
    int32_T c,l;
    
    double *cam_C, *cam_A, *cam_H, *cam_V, *cam_O, *cam_R, *cam_Hd, *cam_Vd;
    double hc,vc,hs,vs;
    
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
    double *PmCys,*APmCys,*HPmCys,*VPmCys;
    double apmc,hpmc,vpmc;
    
    double pmc_o_mag;
    double pmcnox,pmcnoy,pmcnoz;
    double lamx,lamy,lamz;
    double lam_mag2;
    double mup1;
    double pdmcx,pdmcy,pdmcz;
    double apdmc,hpdmc,vpdmc;
    double S_imm05,L_imm05;
    double x_im,y_im;
    
    double s_wdth, l_wdth;
    int32_T l_min,l_max,c_min,c_max;
    int32_T cc,ll;
    
    double x_im_ap,y_im_ap;
    double resolhxy,resolvxy;
    double resolz;
    float elev_amax,elev_amin;
    float resol_amax, resol_amin;
    double mrgnh,mrgnv;
    float data_ignore_value_float;
    double Hd2_abs,Vd2_abs;
    
    cam_C = cahvor_mdl.C; cam_A = cahvor_mdl.A; cam_H = cahvor_mdl.H; cam_V = cahvor_mdl.V;
    cam_O = cahvor_mdl.O; cam_R = cahvor_mdl.R;
    hs = cahvor_mdl.hs; vs = cahvor_mdl.vs; hc = cahvor_mdl.hc; vc = cahvor_mdl.vc;
    cam_Hd = cahvor_mdl.Hdash; cam_Vd = cahvor_mdl.Vdash;
    
    S_dem = (int32_T) msldem_hdr.samples;
    L_dem = (int32_T) msldem_hdr.lines;
    L_demm1 = L_dem - 1;
    S_demp1 = S_dem + 1;
    
    S_imm05 = (double)S_im - 0.5;
    L_imm05 = (double)L_im - 0.5;
    s_wdth = srange[1] - srange[0];
    l_wdth = lrange[1] - lrange[0];
    
    // printf("%d \n",L_dem);
    
    Hd2_abs = fabs(cam_Hd[2]);
    Vd2_abs = fabs(cam_Vd[2]);
    
    // printf("%d \n",L_dem);
    
    data_ignore_value_float = (float) msldem_hdr.data_ignore_value + 1.0;
    // printf("ignore: %f\n",data_ignore_value_float);
    
    // printf("%d \n",L_dem);
    
    
    // printf("size of size_t: %d\n",sz1);
    PmCys = (double*) malloc((size_t) S_dem * sizeof(double));
    APmCys = (double*) malloc((size_t) S_dem * sizeof(double));
    HPmCys = (double*) malloc((size_t) S_dem * sizeof(double));
    VPmCys = (double*) malloc((size_t) S_dem * sizeof(double));
    for(c=0;c<S_dem;c++){
        pmcy  = msldem_easting[c] - cam_C[1];
        PmCys[c] = pmcy;
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
            if(isnan(dem_cl)){
                msldem_imFOVmaskd[c][l] = -1;
            } else {
                pmcz   = -dem_cl-cam_C[2];
                apmcz = cam_A[2] * pmcz;
                apmcy = APmCys[c];
                apmc = apmcx + apmcy + apmcz;
                if(apmc>0){
                    // hpmcy = HPmCys[c];
                    // vpmcy = VPmCys[c];
                    // hpmcz = cam_H[2] * pmcz;
                    // vpmcz = cam_V[2] * pmcz;
                    // hpmc = hpmcx + hpmcy + hpmcz;
                    // vpmc = vpmcx + vpmcy + vpmcz;
                    // x_im_ap = hpmc / apmc;
                    // y_im_ap = vpmc / apmc;
                    x_im_ap = (hpmcx + HPmCys[c] + cam_H[2] * pmcz) / apmc;
                    y_im_ap = (vpmcx + VPmCys[c] + cam_V[2] * pmcz) / apmc;
                    
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

                    //if(c==6140 && l==53979){
                    //    
                    //    printf("dem_c-1l-1 = %f\n",elevlm1[c]);
                    //    printf("dem_cl-1 = %f\n",elevlm1[c+1]);
                    //    printf("dem_c+1l-1 = %f\n",elevlm1[c+2]);
                    //    printf("dem_c-1l = %f\n",elevl[c]);
                    //    printf("dem_cl = %f\n",elevl[c+1]);
                    //    printf("dem_c+1l = %f\n",elevl[c+2]);
                    //    printf("dem_c-1l+1 = %f\n",elevlp1[c]);
                    //    printf("dem_cl+1 = %f\n",elevlp1[c+1]);
                    //    printf("dem_c+1l+1 = %f\n",elevlp1[c+2]);
                    //    printf("resolz = %f\n",resolz);
                    //}

                    mrgnh = coef_mrgn*hs/apmc * (resolhxy + Hd2_abs * resolz);
                    mrgnv = coef_mrgn*vs/apmc * (resolvxy + Vd2_abs * resolz);
                    
                    if (x_im_ap>srange[0]-mrgnh && x_im_ap<srange[1]+mrgnh && 
                            y_im_ap>lrange[0]-mrgnv && y_im_ap<lrange[1]+mrgnv){
                        if(mrgnh>s_wdth || mrgnv>l_wdth){
                            msldem_imFOVmaskd[c][l] = 3;
                        } else {
                            msldem_imFOVmaskd[c][l] = 2;
                        }
                    }
                    
                    if(msldem_imFOVmaskd[c][l]>1){
                        pmcy   = PmCys[c];
                        /*************************************************************/
                        /********** Projection with the CAHVOR model *****************/
                        /*************************************************************/
                        /* Calculation of the apparent P vector */
                        /* first normalization with respect to O */
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
                        if(lam_mag2<1){ /* not consider wide angle */
                            mup1 = 1+cam_R[0] + cam_R[1]*lam_mag2 + cam_R[2]*lam_mag2*lam_mag2;
                            pdmcx = cam_O[0] + mup1*lamx;
                            pdmcy = cam_O[1] + mup1*lamy;
                            pdmcz = cam_O[2] + mup1*lamz;

                            /* perform projection */
                            apdmc = cam_A[0] * pdmcx + cam_A[1] * pdmcy + cam_A[2] * pdmcz;
                            hpdmc = cam_H[0] * pdmcx + cam_H[1] * pdmcy + cam_H[2] * pdmcz;
                            vpdmc = cam_V[0] * pdmcx + cam_V[1] * pdmcy + cam_V[2] * pdmcz;
                            x_im = hpdmc/apdmc; y_im = vpdmc/apdmc;
                            /*************************************************************/

                            if(x_im>-0.5 && x_im<S_imm05 && y_im>-0.5 && y_im<L_imm05){
                                /* if the pixel is within the range, there needs no more 
                                 * evaluation */
                                msldem_imFOVmaskd[c][l] = 5;
                            }
                        }
                    }
                } else {
                    /* Evaluation for apmc<0 */
                }
            }
        }

    }
    
    
    /* Last step: complementation of the surroundings */
    /* This step is costly, since it is performed on the whole image 
     * can be speeded up in a more elaborated way. */
    for(l=0;l<L_dem;l++){
        l_min = (l-1>0) ? (l-1) : 0;
        l_max = (l+2<L_dem) ? (l+2) : L_dem;
        for(c=0;c<S_dem;c++){
            if(msldem_imFOVmaskd[c][l]==0 || msldem_imFOVmaskd[c][l]==2 || msldem_imFOVmaskd[c][l]==3){
                c_min = (c-1>0) ? (c-1) : 0;
                c_max = (c+2<S_dem) ? (c+2) : S_dem;
                for(cc=c_min;cc<c_max;cc++){
                    for(ll=l_min;ll<l_max;ll++){
                        if(msldem_imFOVmaskd[cc][ll]==5)
                            msldem_imFOVmaskd[c][l] = 4;
                    }
                }
            }
            
        }
    }
    
    /* replace -1 with 0 */
    /* This step is moderetely costly, since it is performed on the whole image 
     * can be speeded up in a more elaborated way. */
    for(l=0;l<L_dem;l++){
        for(c=0;c<S_dem;c++){
            if(msldem_imFOVmaskd[c][l]==-1){
                msldem_imFOVmaskd[c][l] = 0;
            }
        }
    }
    
    
    free(elevlm1);
    free(elevl);
    free(elevlp1);
    free(PmCys);
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
    CAHVOR_MODEL cahvor_mdl;
    double *msldem_northing;
    double *msldem_easting;
    // bool **msldem_imFOVmask;
    int8_T **msldem_imFOVmaskd;
    int32_T S_im,L_im;
    double *srange,*lrange;
//     mxArray *cam_C_mxar, *cam_C_mxard;
//     mxArray *cam_A_mxar, *cam_A_mxard;
//     mxArray *cam_H_mxar, *cam_H_mxard;
//     mxArray *cam_V_mxar, *cam_V_mxard;
//     mxArray *cam_Hd_mxar, *cam_Hd_mxard;
//     mxArray *cam_Vd_mxar, *cam_Vd_mxard;
//     mxArray *cam_hc_mxar, *cam_hc_mxard;
//     mxArray *cam_vc_mxar, *cam_vc_mxard;
//     mxArray *cam_hs_mxar, *cam_hs_mxard;
//     mxArray *cam_vs_mxar, *cam_vs_mxard;
//     mxArray *cam_O_mxar, *cam_O_mxard;
//     mxArray *cam_R_mxar, *cam_R_mxard;
//     double *cam_C, *cam_A, *cam_H, *cam_V, *cam_O, *cam_R, *cam_Hd, *cam_Vd;
//     double cam_hc,cam_vc,cam_hs,cam_vs;
    
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
    S_im = (int32_T) mxGetScalar(prhs[4]);
    L_im = (int32_T) mxGetScalar(prhs[5]);
    srange = mxGetDoubles(prhs[6]);
    lrange = mxGetDoubles(prhs[7]);
    //printf("sim = %d\n",S_im);
    
    /* INPUT 6 camera model */
    cahvor_mdl = mxGet_CAHVOR_MODEL(prhs[8]);
    
    // printf("cam_Hd %f %f %f\n",cam_Hd[0],cam_Hd[1],cam_Hd[2]);
    coef_mrgn = mxGetScalar(prhs[9]);
    

    /* OUTPUT 1 msldem imFOV */
    plhs[0] = mxCreateNumericMatrix((mwSize) msldem_hdr.lines,(mwSize) msldem_hdr.samples,mxINT8_CLASS,mxREAL);
    msldem_imFOVmaskd = set_mxInt8Matrix(plhs[0]);
    
    // Initialize matrices
    msldem_samples = (int32_T) msldem_hdr.samples;
    msldem_lines = (int32_T) msldem_hdr.lines;
    //for(si=0;si<msldem_samples;si++){
    //    for(li=0;li<msldem_lines;li++){
    //        msldem_imFOVmaskd[si][li] = 0;
    //    }
    //}
    // printf("sim = %d\n",S_im);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    get_imFOVmask_MSLDEM(msldem_imgpath, msldem_hdr,
        msldem_northing, msldem_easting, 
        S_im, L_im, srange, lrange,cahvor_mdl,
        msldem_imFOVmaskd,coef_mrgn);
    
    /* free memories */
    mxFree(msldem_imgpath);
    // mxFree(msldem_imFOVmask);
    mxFree(msldem_imFOVmaskd);
    
}
