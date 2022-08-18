/* =====================================================================
 * get_msldemtUFOVmask_d_wmsldemc_L2_mex.c
 * Evaluate any pixels in MSLDEM image whether or not they exist in the 
 * 
 * INPUTS:
 * 0 msldemc_img           Double array [L_demc x S_demc]
 * 1 msldemc_northing      Double array [L_demc]
 * 2 msldemc_easting       Double array [S_demc]
 * 3 msldemc_imFOVmask     int8_T [L_demc x S_demc]
 * 4 S_im                  int
 * 5 L_im                  int
 * 6 cahv_mdl              CAHV_MODEL
 * 7 msldemt_img           Double array [L_demt x S_demt]
 * 8 msldemt_northing      Double array [L_demt]
 * 9 msldemt_easting       Double array [S_demt]
 * 10 msldemt_imFOVmask    int8_T [L_demt x S_demt]
 *
 * The origin of msldemc_img, msldemc_northing, and msldemc_easting is 
 * cam_C.
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
#include "lib_proj_mastcamMSLDEM_L2.h"

void bin_msldemt_xyz_wAHVint_L2_d(int32_T S_im, int32_T L_im, CAHV_MODEL cahv_mdl,
        char *msldem_imgpath, EnviHeader msldem_hdr, 
        int32_T msldemc_imxy_sample_offset, int32_T msldemc_imxy_line_offset,
        int32_T msldemc_samples, int32_T msldemc_lines,
        double *msldemc_xmc, double *msldemc_ymc, int8_T **msldemc_imFOVmask,
        double *msldemt_xmc, double *msldemt_ymc, int8_T **msldemt_imFOVmask,
        int32_T **bin_count_im,int32_T ***bin_im_c, int32_T ***bin_im_l,
        double ***bin_imx, double ***bin_imy,double ***bin_demzmc,
        int32_T *count_napmc, int32_T **c_napmc, int32_T **l_napmc, double **zmc_napmc)
{
    int32_T xi,yi,c,l;
    int32_T S_imm1,L_imm1;
    double pmcx,pmcy,pmcz,apmc,ppvx,ppvy;
    double apmcx,hpmcx,vpmcx;
    double *cam_C, *cam_A, *cam_H, *cam_V;
    
    long skip_pri;
    long skip_l, skip_r;
    float *elevl,*elevlp1;
    size_t ncpy;
    size_t sz=sizeof(float);
    FILE *fid;
    float data_ignore_value_float;
    double dem_cl;
    int32_T msldemc_samplesm1, msldemc_linesm1;
    double *APmCys,*HPmCys,*VPmCys;
    
    msldemc_samplesm1 = msldemc_samples - 1;
    msldemc_linesm1 = msldemc_lines - 1;
    
    S_imm1 = S_im - 1;
    L_imm1 = L_im - 1;
    
    cam_C = cahv_mdl.C; cam_A = cahv_mdl.A; cam_H = cahv_mdl.H; cam_V = cahv_mdl.V;
    
    *count_napmc = 0;
    
    // APmCys = (double*) malloc((size_t) msldemc_samplesm1 * sizeof(double));
    // for(c=0;c<msldemc_samplesm1;c++){
    //     msldemt_easting[c] = 0.5*(msldemc_easting[c]+msldemc_easting[c+1]);
    //     pmcy  = msldemt_easting[c] - cam_C[1];
    //     APmCys[c] = cam_A[1] * pmcy;
    // }
    
    /* Initialize the counting matrix */
    for(xi=0;xi<S_im;xi++){
        for(yi=0;yi<L_im;yi++){
            bin_count_im[xi][yi] = 0;
        }
    }
    
    APmCys = (double*) malloc((size_t) msldemc_samplesm1 * sizeof(double));
    HPmCys = (double*) malloc((size_t) msldemc_samplesm1 * sizeof(double));
    VPmCys = (double*) malloc((size_t) msldemc_samplesm1 * sizeof(double));
    for(c=0;c<msldemc_samplesm1;c++){
        pmcy = 0.5*(msldemc_ymc[c]+msldemc_ymc[c+1]);
        msldemt_ymc[c] = pmcy;
        APmCys[c] = cam_A[1] * pmcy;
        HPmCys[c] = cam_H[1] * pmcy;
        VPmCys[c] = cam_V[1] * pmcy;
    }
    
    fid = fopen(msldem_imgpath,"rb");
    /* skip lines */
    skip_pri = (long) msldem_hdr.samples * (long) msldemc_imxy_line_offset * (long) sz;
    // printf("%d*%d*%d=%ld\n",msldem_header.samples,msldemc_imxy_line_offset,s,skip_pri);
    fseek(fid,skip_pri,SEEK_CUR);
    /* read the data */
    ncpy = sz * (size_t) msldemc_samples;
    elevl = (float*) malloc(ncpy);
    elevlp1 = (float*) malloc(ncpy);
    skip_l = (long) sz * (long) msldemc_imxy_sample_offset;
    skip_r = ((long) msldem_hdr.samples - (long) msldemc_samples)* (long) sz - skip_l;
    data_ignore_value_float = (float) msldem_hdr.data_ignore_value + 1.0;
    
    fseek(fid,skip_l,SEEK_CUR);
    fread(elevlp1,sz,msldemc_samples,fid);
    fseek(fid,skip_r,SEEK_CUR);
    for(c=0;c<msldemc_samples;c++){
        if(elevlp1[c]<data_ignore_value_float)
            elevlp1[c] = NAN;
    }
    
    // printf("%d\n",msldemc_samples);
    for(l=0;l<msldemc_linesm1;l++){
        /* read a line */ 
        memcpy(elevl,elevlp1,ncpy);
        
        fseek(fid,skip_l,SEEK_CUR);
        fread(elevlp1,sz,msldemc_samples,fid);
        fseek(fid,skip_r,SEEK_CUR);
        for(c=0;c<msldemc_samples;c++){
            if(elevlp1[c]<data_ignore_value_float)
                elevlp1[c] = NAN;
        }
        
        pmcx = 0.5*(msldemc_xmc[l] + msldemc_xmc[l+1]);
        msldemt_xmc[l] = pmcx;
        apmcx = cam_A[0]*pmcx;
        hpmcx = cam_H[0]*pmcx;
        vpmcx = cam_V[0]*pmcx;
        for(c=0;c<msldemc_samplesm1;c++){
            if((msldemc_imFOVmask[c][l]>0) || (msldemc_imFOVmask[c+1][l]>0) || (msldemc_imFOVmask[c][l+1]>0) || (msldemc_imFOVmask[c+1][l+1]>0)){
                pmcz = -0.5*( (double) elevl[c+1] + (double) elevlp1[c] ) - cam_C[2];
                // pmcx = ppvgx-cam_C[0]; pmcy = ppvgy-cam_C[1]; pmcz = ppvgz-cam_C[2];
                apmc =  apmcx + APmCys[c] + pmcz*cam_A[2];
                ppvx = (hpmcx + HPmCys[c] + pmcz*cam_H[2])/apmc;
                ppvy = (vpmcx + VPmCys[c] + pmcz*cam_V[2])/apmc;
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
                
                if(apmc>0){
                    ++bin_count_im[xi][yi];
                    msldemt_imFOVmask[c][l] = 2;
                } else if(apmc<0) {
                    (*count_napmc)++;
                    msldemt_imFOVmask[c][l] = 1;
                }
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
                bin_demzmc[xi][yi] = (double*) malloc(bin_count_im[xi][yi]*sizeof(double));
            } else {
                bin_im_c[xi][yi] = NULL;
                bin_im_l[xi][yi] = NULL;
                bin_imx[xi][yi] = NULL;
                bin_imy[xi][yi] = NULL;
                bin_demzmc[xi][yi] = NULL;
            }
        }
    }
    
    /* */
    if(*count_napmc>0){
        *c_napmc = (int32_T*) malloc(sizeof(int32_T) * (size_t) *count_napmc);
        *l_napmc = (int32_T*) malloc(sizeof(int32_T) * (size_t) *count_napmc);
        *zmc_napmc = (double*) malloc(sizeof(double) * (size_t) *count_napmc);
    } else {
        *c_napmc = NULL;
        *l_napmc = NULL;
        *zmc_napmc = NULL;
    }
    *count_napmc = 0;
    
    for(xi=0;xi<S_im;xi++){
        for(yi=0;yi<L_im;yi++){
            bin_count_im[xi][yi] = 0;
        }
    }
    
    fseek(fid,skip_pri,SEEK_SET);
    fseek(fid,skip_l,SEEK_CUR);
    fread(elevlp1,sz,msldemc_samples,fid);
    fseek(fid,skip_r,SEEK_CUR);
    for(c=0;c<msldemc_samples;c++){
        if(elevlp1[c]<data_ignore_value_float)
            elevlp1[c] = NAN;
    }
    
    for(l=0;l<msldemc_linesm1;l++){
        memcpy(elevl,elevlp1,ncpy);
        
        fseek(fid,skip_l,SEEK_CUR);
        fread(elevlp1,sz,msldemc_samples,fid);
        fseek(fid,skip_r,SEEK_CUR);
        for(c=0;c<msldemc_samples;c++){
            if(elevlp1[c]<data_ignore_value_float)
                elevlp1[c] = NAN;
        }
        pmcx = msldemt_xmc[l];
        apmcx = cam_A[0]*pmcx;
        hpmcx = cam_H[0]*pmcx;
        vpmcx = cam_V[0]*pmcx;
        for(c=0;c<msldemc_samplesm1;c++){
            if(msldemt_imFOVmask[c][l]>1){
                pmcz = -0.5*( (double) elevl[c+1] + (double) elevlp1[c] ) - cam_C[2];
                // pmcx = ppvgx-cam_C[0]; pmcy = ppvgy-cam_C[1]; pmcz = ppvgz-cam_C[2];
                apmc =  apmcx + APmCys[c] + pmcz*cam_A[2];
                ppvx = (hpmcx + HPmCys[c] + pmcz*cam_H[2])/apmc;
                ppvy = (vpmcx + VPmCys[c] + pmcz*cam_V[2])/apmc;
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
                bin_demzmc[xi][yi][bin_count_im[xi][yi]]  = pmcz;
                ++bin_count_im[xi][yi];

            } else if(msldemt_imFOVmask[c][l]==1){
                pmcz = -0.5*( (double) elevl[c+1] + (double) elevlp1[c] ) - cam_C[2];
                (*c_napmc)[*count_napmc] = c;
                (*l_napmc)[*count_napmc] = l;
                (*zmc_napmc)[*count_napmc] = pmcz;
                (*count_napmc)++;
            }
        }
    }
    
    fclose(fid);
    free(elevl);
    free(elevlp1);
    free(APmCys);
    free(HPmCys);
    free(VPmCys);
}
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *msldem_imgpath;
    EnviHeader msldem_header;
    mwSize msldemc_imxy_sample_offset,msldemc_imxy_line_offset;
    mwSize msldemc_samples,msldemc_lines;
    double *msldemc_xmc;
    double *msldemc_ymc;
    CAHV_MODEL cahv_mdl;
    int8_T **msldemc_imFOVmask;
    mwSize S_im,L_im;
    
    double *msldemt_xmc;
    double *msldemt_ymc;
    int8_T **msldemt_imUFOVmask;
    
    mwIndex si,li;
    
    int32_T **bin_count_im, *bin_count_im_base;
    int32_T ***bin_im_c, **bin_im_c_base;
    int32_T ***bin_im_l, **bin_im_l_base;
    double ***bin_imx, **bin_imx_base;
    double ***bin_imy, **bin_imy_base;
    double ***bin_demzmc, **bin_demzmc_base;
    int32_T *c_napmc, *l_napmc;
    int32_T count_napmc;
    double *zmc_napmc;

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
    msldem_header = mxGetEnviHeader(prhs[1]);
    
    /* INPUT 2 msldemc_sheader*/
    msldemc_imxy_sample_offset = (mwSize) mxGetScalar(mxGetField(prhs[2],0,"sample_offset"));
    msldemc_imxy_line_offset = (mwSize) mxGetScalar(mxGetField(prhs[2],0,"line_offset"));
    msldemc_samples = (mwSize) mxGetScalar(mxGetField(prhs[2],0,"samples"));
    msldemc_lines = (mwSize) mxGetScalar(mxGetField(prhs[2],0,"lines"));
    //msldemc_img = set_mxDoubleMatrix(prhs[0]);
    //L_demc = mxGetM(prhs[0]);
    //S_demc = mxGetN(prhs[0]);
    
    /* INPUT 1/2 msldem northing easting */
    msldemc_xmc = mxGetDoubles(prhs[3]);
    msldemc_ymc = mxGetDoubles(prhs[4]);
    
    /* INPUT 3 msldemc imFOV */
    msldemc_imFOVmask = set_mxInt8Matrix(prhs[5]);
    
    /* INPUT 4/5 image S_im, L_im */
    S_im = (mwSize) mxGetScalar(prhs[6]);
    L_im = (mwSize) mxGetScalar(prhs[7]);
    
    /* INPUT 6 CAHV model */
    cahv_mdl = mxGet_CAHV_MODEL(prhs[8]);
    
    
    /* OUTPUT 0 msldemc imFOV */
    plhs[0] = mxCreateNumericMatrix(msldemc_lines-1,msldemc_samples-1,mxINT8_CLASS,mxREAL);
    msldemt_imUFOVmask = set_mxInt8Matrix(plhs[0]);
    
    // Initialize matrices
    for(si=0;si<msldemc_samples-1;si++){
        for(li=0;li<msldemc_lines-1;li++){
            msldemt_imUFOVmask[si][li] = 0;
        }
    }
    
    /* -----------------------------------------------------------------
     * binning the image
     * ----------------------------------------------------------------- */
    createInt32Matrix(&bin_count_im,&bin_count_im_base,(size_t) S_im, (size_t) L_im);
    createInt32PMatrix(&bin_im_c, &bin_im_c_base, (size_t) S_im, (size_t) L_im);
    createInt32PMatrix(&bin_im_l, &bin_im_l_base, (size_t) S_im, (size_t) L_im);
    createDoublePMatrix(&bin_imx, &bin_imx_base, (size_t) S_im, (size_t) L_im);
    createDoublePMatrix(&bin_imy, &bin_imy_base, (size_t) S_im, (size_t) L_im);
    createDoublePMatrix(&bin_demzmc, &bin_demzmc_base, (size_t) S_im, (size_t) L_im);
    
    msldemt_xmc = (double*) malloc(sizeof(double)*(size_t) (msldemc_lines-1));
    msldemt_ymc = (double*) malloc(sizeof(double)*(size_t) (msldemc_samples-1));
    
    bin_msldemt_xyz_wAHVint_L2_d((int32_T) S_im, (int32_T) L_im, cahv_mdl,
            msldem_imgpath, msldem_header, 
            msldemc_imxy_sample_offset, msldemc_imxy_line_offset,
            msldemc_samples, msldemc_lines, msldemc_xmc, msldemc_ymc, msldemc_imFOVmask,
            msldemt_xmc, msldemt_ymc, msldemt_imUFOVmask,
            bin_count_im, bin_im_c, bin_im_l, bin_imx, bin_imy, bin_demzmc,
            &count_napmc,&c_napmc,&l_napmc,&zmc_napmc);
    
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */    
    mask_obstructed_pts_in_msldemt_using_msldemc_L2(msldem_imgpath, msldem_header, 
        (int32_T) msldemc_imxy_sample_offset, (int32_T) msldemc_imxy_line_offset,
        (int32_T) msldemc_samples, (int32_T) msldemc_lines,
        msldemc_xmc, msldemc_ymc, msldemc_imFOVmask,
        msldemt_xmc, msldemt_ymc, msldemt_imUFOVmask,
        bin_count_im, bin_im_c, bin_im_l,
        bin_imx, bin_imy, bin_demzmc,
        count_napmc, c_napmc, l_napmc, zmc_napmc,
        S_im, L_im, cahv_mdl);
        
    
    /* Freeing memories */
    for(si=0;si<S_im;si++){
        for(li=0;li<L_im;li++){
            if(bin_count_im[si][li]>0){
                free(bin_im_c[si][li]);
                free(bin_im_l[si][li]);
                free(bin_imx[si][li]);
                free(bin_imy[si][li]);
                free(bin_demzmc[si][li]);
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
    free(bin_demzmc);
    free(bin_demzmc_base);
    
    free(bin_count_im);
    free(bin_count_im_base);
    
    if(count_napmc>0){
        free(c_napmc);
        free(l_napmc);
        free(zmc_napmc);
    }
    
    free(msldemt_xmc);
    free(msldemt_ymc);
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldemc_imFOVmask);
    mxFree(msldemt_imUFOVmask);
    
}