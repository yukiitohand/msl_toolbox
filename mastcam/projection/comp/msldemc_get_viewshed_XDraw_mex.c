/* =====================================================================
 * viewshed_get_msldemtUFOVmask_ctr_XDraw_mex.c
 * Get viewshed within the FOV using XDraw algorithm
 * 
 * INPUTS:
 * 0 msldem_imgpath        char* path to the image
 * 1 msldem_hdr            EnviHeader
 * 2 msldemc_imFOVhdr      Struct
 * 3 msldemc_imFOVmask     int8_t [L_demc x S_demc]
 * 4 cahv_mdl              CAHV_MODEL
 *
 * 
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
#include "msldem_util.h"
#include <time.h>


void get_viewshed_XDraw(char *msldem_imgpath, EnviHeader msldem_hdr, 
        MSLDEMC_HEADER msldemc_hdr,CAHV_MODEL cahv_mdl, 
        int8_t **msldemc_imFOVmask, int8_t **msldemc_imUFOVmask)
{
    int32_t si,li,c;
    int32_t msldemc_samples, msldemc_lines;
    
    long skip_pri;
    long skip_l, skip_r;
    float *elevl;
    size_t ncpy;
    size_t sz=sizeof(float);
    FILE *fid;
    float data_ignore_value_float;
    
    double **local_rim, *local_rim_base;
    double *cam_C;
    double d[3];
    double si_dbl,li_dbl;
    double slp,y_lm1;
    double si_ref_left,si_ref_right;
    int32_t ref_grid_si_1,ref_grid_si_2;
    double elevcl_dbl,local_rim_sl;
    
    cam_C = cahv_mdl.C;
    
    msldemc_samples = msldemc_hdr.samples;
    msldemc_lines   = msldemc_hdr.lines;
    
    createDoubleMatrix(&local_rim, &local_rim_base, (size_t) msldemc_samples, (size_t) msldemc_lines);
    for(si=0;si<msldemc_samples;si++){
        for(li=0;li<msldemc_lines;li++){
            local_rim[si][li] = -32767.0;
        }
    }
    
    fid = fopen(msldem_imgpath,"rb");
    /* skip lines */
    skip_pri = (long) msldem_hdr.samples * (long) msldemc_hdr.line_offset * (long) sz;
    // printf("%d*%d*%d=%ld\n",msldem_header.samples,msldemc_imxy_line_offset,s,skip_pri);
    fseek(fid,skip_pri,SEEK_CUR);
    /* read the data */
    ncpy = sz * (size_t) msldemc_samples;
    elevl = (float*) malloc(ncpy);
    skip_l = (long) sz * (long) msldemc_hdr.sample_offset;
    skip_r = ((long) msldem_hdr.samples - (long) msldemc_samples)* (long) sz - skip_l;
    data_ignore_value_float = (float) msldem_hdr.data_ignore_value + 1.0;
    
    
    
    for(li=0;li<msldemc_lines;li++){
        fseek(fid,skip_l,SEEK_CUR);
        fread(elevl,sz,msldemc_samples,fid);
        fseek(fid,skip_r,SEEK_CUR);
        for(c=0;c<msldemc_samples;c++){
            if(elevl[c]<data_ignore_value_float)
                elevl[c] = NAN;
        }
        for(si=0;si<msldemc_samples;si++){
            if(msldemc_imFOVmask[si][li]>0){
                elevcl_dbl = (double) elevl[si];
                /* Get the angle between C  */
                li_dbl = (double) li; si_dbl = (double) si;
                d[0] = li_dbl - cam_C[0];
                d[1] = si_dbl - cam_C[1];
                d[2] = elevcl_dbl - cam_C[2];
                
                slp = d[1]/d[0];
                y_lm1 = si_dbl - slp;
                
                si_ref_left  = floor(y_lm1);
                si_ref_right = ceil(y_lm1);
                
                ref_grid_si_1 = (int32_t) si_ref_left;
                ref_grid_si_2 = (int32_t) si_ref_right;
                
                /* Compute the coordiante of the intersection of the 2d line
                 * of sight and grid. */
                
                /* z coordinate of the line of sight at l=li-1 */
                if(li==0 || si==0 || si==msldemc_samples-1){
                    local_rim[si][li] = elevcl_dbl;
                
                } else {
                    // local_rim_sl = fmax(local_rim[ref_grid_si_1][li-1],local_rim[ref_grid_si_2][li-1]);
                    local_rim_sl = local_rim[ref_grid_si_1][li-1] * (si_ref_right-y_lm1) +  local_rim[ref_grid_si_2][li-1] * (y_lm1-si_ref_left);
                    local_rim_sl = (local_rim_sl-cam_C[2]) * d[0] / (d[0]-1) + cam_C[2];

                    if(elevcl_dbl>local_rim_sl){
                        local_rim[si][li] = elevcl_dbl;
                    } else {
                        local_rim[si][li] =local_rim_sl;
                        msldemc_imUFOVmask[si][li] = 0;
                    }
                }
            }


        }
    }
    fclose(fid);
    
    free(local_rim);
    free(local_rim_base);
    free(elevl);
    
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *msldem_imgpath;
    EnviHeader msldem_hdr;
    MSLDEMC_HEADER msldemc_hdr;
    mwSize msldemc_samples,msldemc_lines;
    CAHV_MODEL cahv_mdl;
    int8_t **msldemc_imFOVmask;
    int8_t **msldemt_imUFOVmask;
    
    clock_t strt_time, end_time;
    double cpu_time_used;
    mwSize si,li;

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
    
    /* INPUT 1 msldem_header and 2 msldemc_hdr */
    msldem_hdr = mxGetEnviHeader(prhs[1]);
    mxGet_MSLDEMC_HEADER(prhs[2],&msldemc_hdr);
    msldemc_samples = (mwSize) msldemc_hdr.samples;
    msldemc_lines   = (mwSize) msldemc_hdr.lines;
    
    /* INPUT 3 msldemc imFOV */
    msldemc_imFOVmask = set_mxInt8Matrix(prhs[3]);
    
    /* INPUT 4 CAHV model */
    cahv_mdl = mxGet_CAHV_MODEL(prhs[4]);
    
    /* OUTPUT 0 msldemc imFOV */
    plhs[0] = mxCreateNumericMatrix(msldemc_lines,msldemc_samples,mxINT8_CLASS,mxREAL);
    msldemt_imUFOVmask = set_mxInt8Matrix(plhs[0]);
    
    // Initialize matrices
    for(si=0;si<msldemc_samples;si++){
        for(li=0;li<msldemc_lines;li++){
            if(msldemc_imFOVmask[si][li]>1){
                msldemt_imUFOVmask[si][li] = 2;
            } else {
                msldemt_imUFOVmask[si][li] = msldemc_imFOVmask[si][li];
            }
        }
    }
    printf("a\n");
    /* Main computation routine */
    get_viewshed_XDraw(msldem_imgpath, msldem_hdr, msldemc_hdr, cahv_mdl, 
        msldemc_imFOVmask, msldemt_imUFOVmask);
    
    
    mxFree(msldem_imgpath);
    mxFree(msldemc_imFOVmask);
    mxFree(msldemt_imUFOVmask);

    
}