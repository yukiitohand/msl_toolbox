/* =====================================================================
 * iaumars_get_msldemtUFOVmask_d_L2PBK_LL0_M3_mex.c
 * L2  : msldemc will be read from a file not an input.
 * PBK : Prior Binning into bins with the auxiliary size defined by two parameters K_L and K_S
 * LL0 : Linked List with least basic information (c,l,radius), this is most memory efficient
 * M3  : 3x3 matrix inversion object image coordinate.
 * 
 * INPUTS:
 * 0 msldem_imgpath        char* path to the image
 * 1 msldem_hdr            EnviHeader
 * 2 mslrad_offset         Double Scalar
 * 3 msldemc_imFOVhdr      Struct
 * 4 msldemc_latitude      Double array [L_demc]
 * 5 msldemc_longitude     Double array [S_demc]
 * 6 msldemc_imFOVmask     int8_t [L_demc x S_demc]
 * 7 S_im                  int
 * 8 L_im                  int
 * 9 cahv_mdl              CAHV_MODEL
 * 10 K_L                  double reciprocal of the length of the bin in the image line direction.
 * 11 K_S                  double reciprocal of the length of the bin in the image sample direction.
 * 12 dyu                  int8_t flag for dynamic linked list update.
 *
 * 
 * 
 * 
 * OUTPUTS:
 * 0  msldemt_inImage     [(L_demc-1) x (S_demc-1)]  Boolean
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
#include "lib_proj_mastcamMSLDEM_IAUMars_L2PBK_LL0_M3.h"

#include <time.h>



void bin_msldemt_d_iaumars_L2PBK_LL0(double S_im, double L_im, CAHV_MODEL cahv_mdl,
        char *msldem_imgpath, EnviHeader msldem_hdr, double mslrad_offset,
        int32_t msldemc_imxy_sample_offset, int32_t msldemc_imxy_line_offset,
        int32_t msldemc_samples, int32_t msldemc_lines,
        double *msldemc_latitude, double *msldemc_longitude, 
        int8_T **msldemc_imFOVmask, 
        double *msldemt_latitude, double *msldemt_longitude, 
        int8_T **msldemt_imUFOVmask,
        double K_L, double K_S,
        struct MSLDEMmask_LinkedList ***ll_papmc,
        struct MSLDEMmask_LinkedList **ll_napmc)
{
    int32_t xi,yi,c,l;
    int32_t binL,binS;
    int32_t binLm1,binSm1;
    double pmcx,pmcy,pmcz,apmc,ppvx,ppvy;
    double *cam_C, *cam_A, *cam_H, *cam_V;
    
    int32_t msldemc_samplesm1, msldemc_linesm1;
    
    long skip_pri;
    long skip_l, skip_r;
    float *elevl,*elevlp1;
    size_t ncpy;
    size_t sz=sizeof(float);
    FILE *fid;
    float data_ignore_value_float;
    // double dem_cl;
    
    double *cos_tlon, *sin_tlon;
    double tlonc,tlatl;
    double radius_tmp;
    double cos_tlatl, sin_tlatl;
    double x_iaumars, y_iaumars, z_iaumars;
    struct MSLDEMmask_LinkedList *ll_papmc_next;
    struct MSLDEMmask_LinkedList *ll_napmc_next;
    
    msldemc_samplesm1 = msldemc_samples - 1;
    msldemc_linesm1   = msldemc_lines - 1;
    
    // ll_napmc_next = (*ll_napmc);
    
    // S_imm1 = S_im - 1;
    // L_imm1 = L_im - 1;
    
    binL = (int32_t) (K_L * L_im); binS = (int32_t) (K_S * S_im);
    binLm1 = binL-1; binSm1 = binS-1;
    
    cam_C = cahv_mdl.C; cam_A = cahv_mdl.A; cam_H = cahv_mdl.H; cam_V = cahv_mdl.V;
    
    // printf("msldemc_samples = %d\n",msldemc_samples);
    // printf("msldemc_lines   = %d\n",msldemc_lines);
    
    cos_tlon = (double*) malloc(sizeof(double) * (size_t) msldemc_samplesm1);
    sin_tlon = (double*) malloc(sizeof(double) * (size_t) msldemc_samplesm1);
    for(c=0;c<msldemc_samplesm1;c++){
        // printf("c=%d,%d\n",c,msldemc_samplesm1);
        tlonc = 0.5*(msldemc_longitude[c]+msldemc_longitude[c+1]);
        cos_tlon[c] = cos(tlonc);
        sin_tlon[c] = sin(tlonc);
        msldemt_longitude[c] = tlonc;
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
    
    /* Read the first line  */
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
        // printf("l=%d\n",l);
        tlatl = 0.5*(msldemc_latitude[l]+msldemc_latitude[l+1]);
        cos_tlatl = cos(tlatl);
        sin_tlatl = sin(tlatl);
        msldemt_latitude[l] = tlatl;
        for(c=0;c<msldemc_samplesm1;c++){
            if((msldemc_imFOVmask[c][l]>0) || (msldemc_imFOVmask[c+1][l]>0) || 
                    (msldemc_imFOVmask[c][l+1]>0) || (msldemc_imFOVmask[c+1][l+1]>0)){
                radius_tmp = 0.5*( (double) elevl[c+1] + (double) elevlp1[c] ) + mslrad_offset;
                x_iaumars  = radius_tmp * cos_tlatl * cos_tlon[c];
                y_iaumars  = radius_tmp * cos_tlatl * sin_tlon[c];
                z_iaumars  = radius_tmp * sin_tlatl;
                pmcx = x_iaumars - cam_C[0];
                pmcy = y_iaumars - cam_C[1];
                pmcz = z_iaumars - cam_C[2];
                
                apmc =  pmcx*cam_A[0] + pmcy*cam_A[1] + pmcz*cam_A[2];
                if(apmc>0){
                    ppvx = (pmcx*cam_H[0] + pmcy*cam_H[1] + pmcz*cam_H[2])/apmc;
                    ppvy = (pmcx*cam_V[0] + pmcy*cam_V[1] + pmcz*cam_V[2])/apmc;
                    xi = (int32_T) floor(K_S*(ppvx+0.5));
                    yi = (int32_T) floor(K_L*(ppvy+0.5));
                    //xi = xi<0?0:xi;
                    //xi = xi>S_imm1?:xi;
                    if(xi<0)
                        xi=0;
                    else if(xi>binSm1)
                        xi = binSm1;

                    if(yi<0)
                        yi=0;
                    else if(yi>binLm1)
                        yi = binLm1;

                    ll_papmc_next = malloc(sizeof(struct MSLDEMmask_LinkedList));
                    ll_papmc_next->c = c;
                    ll_papmc_next->l = l;
                    ll_papmc_next->radius = radius_tmp;
                    ll_papmc_next->next = ll_papmc[xi][yi];
                    ll_papmc[xi][yi] = ll_papmc_next;
                
                    msldemt_imUFOVmask[c][l] = 2;
                    
                } else if(apmc<0) {
                    
                    ll_napmc_next = malloc(sizeof(struct MSLDEMmask_LinkedList));
                    ll_napmc_next->c = c;
                    ll_napmc_next->l = l;
                    ll_napmc_next->radius = radius_tmp;
                    ll_napmc_next->next = (*ll_napmc);
                    (*ll_napmc) = ll_napmc_next;
                    
                    msldemt_imUFOVmask[c][l] = 1;
                }
            } 
        }
    }
    
    fclose(fid);
    free(elevl);
    free(elevlp1);
    free(cos_tlon);
    free(sin_tlon);
}
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *msldem_imgpath;
    EnviHeader msldem_header;
    mwSize msldemc_sample_offset,msldemc_line_offset;
    mwSize msldemc_samples,msldemc_lines;
    double *msldemc_latitude;
    double *msldemc_longitude;
    double mslrad_offset;
    CAHV_MODEL cahv_mdl;
    int8_T **msldemc_imFOVmask;
    double S_im,L_im;
    int8_t dyu;
    
    mwIndex si,li;
    
    mwSize msldemt_samples, msldemt_lines;
    double *msldemt_latitude;
    double *msldemt_longitude;
    int8_t **msldemt_imUFOVmask;
    
    struct MSLDEMmask_LinkedList ***ll_papmc_bin, **ll_papmc_bin_base;
    struct MSLDEMmask_LinkedList *ll_napmc;
    struct MSLDEMmask_LinkedList *ll_tmp;
    
    
    double K_L,K_S;
    mwSize binL,binS;
    
    clock_t strt_time, end_time;
    double cpu_time_used;

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
    
    /* INPUT 1 msldem_header and 2 radius offset */
    msldem_header = mxGetEnviHeader(prhs[1]);
    
    mslrad_offset = mxGetScalar(prhs[2]);
    
    /* INPUT 3 msldemc_sheader*/
    msldemc_sample_offset = (mwSize) mxGetScalar(mxGetField(prhs[3],0,"sample_offset"));
    msldemc_line_offset = (mwSize) mxGetScalar(mxGetField(prhs[3],0,"line_offset"));
    msldemc_samples = (mwSize) mxGetScalar(mxGetField(prhs[3],0,"samples"));
    msldemc_lines = (mwSize) mxGetScalar(mxGetField(prhs[3],0,"lines"));
    //msldemc_img = set_mxDoubleMatrix(prhs[0]);
    //L_demc = mxGetM(prhs[0]);
    //S_demc = mxGetN(prhs[0]);
    
    /* INPUT 4/5 msldem northing easting */
    msldemc_latitude  = mxGetDoubles(prhs[4]);
    msldemc_longitude = mxGetDoubles(prhs[5]);
    
    
    /* INPUT 6 msldemc imFOV */
    msldemc_imFOVmask = set_mxInt8Matrix(prhs[6]);
    
    /* INPUT 7/8 image S_im, L_im */
    S_im = mxGetScalar(prhs[7]);
    L_im = mxGetScalar(prhs[8]);
    
    /* INPUT 9 CAHV model */
    cahv_mdl = mxGet_CAHV_MODEL(prhs[9]);
    
    /* INPUT 10/11 bin size parameters*/
    K_L = mxGetScalar(prhs[10]);
    K_S = mxGetScalar(prhs[11]);
    
    /* INPUT 12 Dynamic Linked List Update FLAG */
    dyu = (int8_t) mxGetScalar(prhs[12]);
    
    
    /* OUTPUT 0 msldemc imUFOV */
    msldemt_samples = msldemc_samples-1;
    msldemt_lines = msldemc_lines-1;
    plhs[0] = mxCreateNumericMatrix(msldemt_lines,msldemt_samples,
            mxINT8_CLASS,mxREAL);
    msldemt_imUFOVmask = set_mxInt8Matrix(plhs[0]);
    
    // plhs[1] = mxCreateNumericMatrix(L_im,S_im,mxINT32_CLASS,mxREAL);
    // bin_count_im = set_mxInt32Matrix(plhs[1]);
    
    // Initialize matrices
    for(si=0;si<msldemt_samples;si++){
        for(li=0;li<msldemt_lines;li++){
            msldemt_imUFOVmask[si][li] = 0;
        }
    }
    
    /* -----------------------------------------------------------------
     * binning the image
     * ----------------------------------------------------------------- */
    ll_napmc = NULL;
    ll_tmp   = NULL;
    strt_time = clock();
    binL = (mwSize) (K_L * L_im); binS = (mwSize) (K_S * S_im);
    createMSLDEMmask_LLMatrix(&ll_papmc_bin, &ll_papmc_bin_base, (size_t) binS, (size_t) binL);
    for(si=0;si<binS;si++){
        for(li=0;li<binL;li++){
            ll_papmc_bin[si][li] = NULL;
        }
    }
    
    msldemt_latitude  = (double*) malloc(sizeof(double)*(size_t) msldemt_lines);
    msldemt_longitude = (double*) malloc(sizeof(double)*(size_t) msldemt_samples);
    
    
    bin_msldemt_d_iaumars_L2PBK_LL0(S_im, L_im, cahv_mdl,
            msldem_imgpath, msldem_header, mslrad_offset,
            (int32_t) msldemc_sample_offset, (int32_t) msldemc_line_offset,
            (int32_t) msldemc_samples, (int32_t) msldemc_lines, 
            msldemc_latitude, msldemc_longitude, msldemc_imFOVmask,
            msldemt_latitude, msldemt_longitude, msldemt_imUFOVmask,
            K_L,K_S,
            ll_papmc_bin, &ll_napmc);
    end_time = clock();
    cpu_time_used = ((double) (end_time - strt_time)) / CLOCKS_PER_SEC;
    printf("Preprocessing took %f [sec].\n",cpu_time_used);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    strt_time = clock();
    if(dyu==1){
        mask_obstructed_pts_in_msldemt_using_msldemc_iaumars_L2PBK_LL0DYU_M3(
                msldem_imgpath, msldem_header, mslrad_offset,
            (int32_T) msldemc_sample_offset, (int32_T) msldemc_line_offset,
            (int32_T) msldemc_samples, (int32_T) msldemc_lines,
            msldemc_latitude, msldemc_longitude, msldemc_imFOVmask,
            (int32_T) msldemt_samples, (int32_T) msldemt_lines,
            msldemt_latitude, msldemt_longitude, msldemt_imUFOVmask,
            K_L,K_S,
            ll_papmc_bin, ll_napmc,
            S_im, L_im, cahv_mdl);
    } else if(dyu==0){
        mask_obstructed_pts_in_msldemt_using_msldemc_iaumars_L2PBK_LL0_M3(
                msldem_imgpath, msldem_header, mslrad_offset,
            (int32_T) msldemc_sample_offset, (int32_T) msldemc_line_offset,
            (int32_T) msldemc_samples, (int32_T) msldemc_lines,
            msldemc_latitude, msldemc_longitude, msldemc_imFOVmask,
            (int32_T) msldemt_samples, (int32_T) msldemt_lines,
            msldemt_latitude, msldemt_longitude, msldemt_imUFOVmask,
            K_L,K_S,
            ll_papmc_bin, ll_napmc,
            S_im, L_im, cahv_mdl);
    } else {
        printf("Undefined DYU_FLAG=%d\n",dyu);
    }
        
    end_time = clock();
    cpu_time_used = ((double) (end_time - strt_time)) / CLOCKS_PER_SEC;
    printf("Main operation took %f [sec].\n",cpu_time_used);
    
    /* Freeing memories */
    for(si=0;si<binS;si++){
        for(li=0;li<binL;li++){
            while(ll_papmc_bin[si][li] != NULL){
                ll_tmp   = ll_papmc_bin[si][li];
                ll_papmc_bin[si][li] = ll_papmc_bin[si][li]->next;
                free(ll_tmp);
            }
        }
    }
    free(ll_papmc_bin);
    free(ll_papmc_bin_base);
    
    while(ll_napmc != NULL){
        ll_tmp   = ll_napmc;
        ll_napmc = ll_napmc->next;
        free(ll_tmp);
    }
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldemc_imFOVmask);
    mxFree(msldemt_imUFOVmask);
    free(msldemt_latitude);
    free(msldemt_longitude);
    
}