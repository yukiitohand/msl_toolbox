/* =====================================================================
 * msldem_restore_areoid_fill_imageblock_mex.c
 * Fill out the image block.
 * 
 * INPUTS:
 * 0 output_imgpath        char*
 * 1 msldem_header         struct
 * 2 sample_offset         scalar
 * 3 line_offset           scalar
 * 4 samples               scalar
 * 5 lines                 scalar
 * 6 subimg                [lines x samples]
 * 
 * OUTPUTS:
 * none
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

/* main computation routine */
void msldem_fill_subimage(char *output_imgpath, EnviHeader msldem_hdr, 
        int32_T sample_offset,int32_T line_offset, 
        int32_T samples, int32_T lines, double **subimg)
{
    int32_T c,l;
    float data_ignore_value_float;
    FILE *fid_out;
    float *out_buf;
    int32_T S_dem, L_dem;
    
    size_t s=sizeof(float);
    long skip_pri;
    long skip_l, skip_r;
    size_t ncpy;
    
    S_dem = (int32_T) msldem_hdr.samples;
    L_dem = (int32_T) msldem_hdr.lines;
    data_ignore_value_float = (float) msldem_hdr.data_ignore_value;
    // printf("ignore: %f\n",data_ignore_value_float);
    
    printf("%d \n",L_dem);

    
    
    fid_out = fopen(output_imgpath,"r+b");
    
    skip_pri = (long) S_dem * (long) line_offset * (long) s;
    fseek(fid_out,skip_pri,SEEK_CUR);
    
    ncpy = (size_t) samples * s;
    skip_l = (long) s * (long) sample_offset;
    skip_r = ((long) S_dem - (long) samples) * (long) s - skip_l;
    
    out_buf = (float*) malloc(ncpy);
    for(c=0;c<samples;c++)
        out_buf[c] = data_ignore_value_float;
    
    for(l=0;l<lines;l++){
        for(c=0;c<samples;c++){
            if isnan(subimg[c][l]){
                out_buf[c] = data_ignore_value_float;
            }
            else{
                out_buf[c] = (float) subimg[c][l];
            }
        }
        fseek(fid_out,skip_l,SEEK_CUR);
        fwrite(out_buf,s,samples,fid_out);
        fseek(fid_out,skip_r,SEEK_CUR);
    }
    
    free(out_buf);
    fclose(fid_out);
    
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *output_imgpath;
    EnviHeader msldem_hdr;
    mwSize sample_offset,line_offset;
    mwSize samples, lines;
    double **subimg;
    

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
    
    /* INPUT 2 output_imgpath */
    output_imgpath = mxArrayToString(prhs[0]);
    
    /* INPUT 1 msldem_header */
    msldem_hdr = mxGetEnviHeader(prhs[1]);
    
    
    
    /* INPUT 3/4/5/6 offset and sub image sizes */
    sample_offset = (mwSize) mxGetScalar(prhs[2]);
    line_offset   = (mwSize) mxGetScalar(prhs[3]);
    samples = (mwSize) mxGetScalar(prhs[4]);
    lines   = (mwSize) mxGetScalar(prhs[5]);
    
    subimg = set_mxDoubleMatrix(prhs[6]);
    
    /* No output */
    
    // printf("sim = %d\n",S_im);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    msldem_fill_subimage(output_imgpath,msldem_hdr,
            (int32_T) sample_offset, (int32_T) line_offset,
            (int32_T) samples, (int32_T) lines,subimg);
    
    /* free memories */
    mxFree(output_imgpath);
    mxFree(subimg);
    
}
