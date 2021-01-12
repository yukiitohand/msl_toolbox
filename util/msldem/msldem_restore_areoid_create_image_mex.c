/* =====================================================================
 * msldem_restore_areoid_create_image_mex.c
 * Create a image data filled with data_ignore_value 
 * 
 * INPUTS:
 * 0 output_imgpath        char*
 * 1 msldem_header         struct
 * 
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
void msldem_create_image(char *output_imgpath, EnviHeader msldem_hdr)
{
    int32_T c,l;
    float data_ignore_value_float;
    FILE *fid_out;
    float *out_buf;
    int32_T S_dem, L_dem;
    
    S_dem = (int32_T) msldem_hdr.samples;
    L_dem = (int32_T) msldem_hdr.lines;
    data_ignore_value_float = (float) msldem_hdr.data_ignore_value;
    // printf("ignore: %f\n",data_ignore_value_float);
    
    printf("%d \n",L_dem);

    out_buf = (float*) malloc(sizeof(float) * (size_t) S_dem);
    for(c=0;c<S_dem;c++){
        out_buf[c] = data_ignore_value_float;
    }
    
    fid_out = fopen(output_imgpath,"wb");
    
    printf("%d \n",L_dem);
    
    
    // printf("%d,%d,%d\n",skip_l,msldemc_samples*s,skip_r);
    for(l=0;l<L_dem;l++){
        fwrite(out_buf,sizeof(float),S_dem,fid_out);
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
    
    /* INPUT 0 output_imgpath */
    output_imgpath = mxArrayToString(prhs[0]);
    
    /* INPUT 1 msldem_header */
    msldem_hdr = mxGetEnviHeader(prhs[1]);
    
    
    /* No output */
    
    // printf("sim = %d\n",S_im);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    msldem_create_image(output_imgpath,msldem_hdr);
    
    /* free memories */
    mxFree(output_imgpath);
    
}
