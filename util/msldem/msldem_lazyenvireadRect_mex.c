/* =====================================================================
 * msldem_lazyenvireadRect_mex.c
 * Read the specified rectangle region of the MSLDEM image data. The 
 * rectangle part of the image:
 * MSLDEM[msldemc_sample_offset:msldemc_sample_offset+msldemc_samples,
 *        msldemc_line_offset:msldemc_line_offset+msldemc_line]
 * is read.
 * 
 * INPUTS:
 * 0 msldem_imgpath        char*
 * 1 msldem_header         struct
 * 2 msldemc_sample_offset integer
 * 3 msldemc_line_offset   integer
 * 4 msldemc_samples       integer
 * 5 msldemc_lines         integer
 * 
 * 
 * OUTPUTS:
 * 0  dem_img [S_demc x L_demc]   Float
 * #Note that the image needs to be flipped after this.
 *
 *
 * This is a MEX file for MATLAB.
 * ===================================================================== */
#include "io64.h"
#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "envi.h"


/* main computation routine */
void msldem_lazyenvireadRect(char *msldem_imgpath, EnviHeader msldem_header, 
        int32_T msldemc_imxy_sample_offset, int32_T msldemc_imxy_line_offset,
        float **img, int32_T msldemc_samples, int32_T msldemc_lines)
{
    int32_T i;
    long skip_pri;
    long skip_l, skip_r;
    float *buf;
    size_t s=sizeof(float);
    FILE *fid;
    size_t ncpy;
    
    
    fid = fopen(msldem_imgpath,"rb");
    
    /* skip lines */
    skip_pri = (long) msldem_header.samples * (long) msldemc_imxy_line_offset * (long) s;
    // printf("%d*%d*%d=%ld\n",msldem_header.samples,msldemc_imxy_line_offset,s,skip_pri);
    fseek(fid,skip_pri,SEEK_CUR);
    
    /* read the data */
    ncpy = (size_t) msldemc_samples* s;
    buf = (float*) malloc(ncpy);
    skip_l = (long) s * (long) msldemc_imxy_sample_offset;
    skip_r = ((long) msldem_header.samples - (long) msldemc_samples)* (long) s - skip_l;
    
    // printf("%d,%d,%d\n",skip_l,msldemc_samples*s,skip_r);
    for(i=0;i < msldemc_lines;i++){
        //printf("%d/%d\n",i,msldemc_lines);
        fseek(fid,skip_l,SEEK_CUR);
        fread(buf,s,msldemc_samples,fid);
        memcpy(img[i],buf,ncpy);
        //for(j=0;j<msldemc_samples;j++){
        //    img[j][i] = buf[j];
        //}
        fseek(fid,skip_r,SEEK_CUR);
        //_fseeki64(fp,skips,SEEK_CUR);
    }
    free(buf);
    fclose(fid);
}



/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *msldem_imgpath;
    EnviHeader msldem_header;
    mwSize msldemc_imxy_sample_offset,msldemc_imxy_line_offset;
    
    float **dem_img;
    mwSize msldemc_samples, msldemc_lines;
    mwSize j;

    /* -----------------------------------------------------------------
     * CHECK PROPER NUMBER OF INPUTS AND OUTPUTS
     * ----------------------------------------------------------------- */
    if(nrhs!=6) {
        mexErrMsgIdAndTxt("msldem_lazyenvireadRect_mex:nrhs","Six inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("msldem_lazyenvireadRect_mex:nlhs","Five outputs required.");
    }
    /* make sure the first input argument is scalar */
    if( !mxIsChar(prhs[0]) ) {
        mexErrMsgIdAndTxt("msldem_lazyenvireadRect_mex:notChar","Input 0 needs to be a string.");
    }
    if( !mxIsStruct(prhs[1]) ) {
        mexErrMsgIdAndTxt("msldem_lazyenvireadRect_mex:notStruct","Input 1 needs to be a struct.");
    }
    if( !mxIsScalar(prhs[2]) ) {
        mexErrMsgIdAndTxt("msldem_lazyenvireadRect_mex:notScalar","Input 2 needs to be a scalar.");
    }
    if( !mxIsScalar(prhs[3]) ) {
        mexErrMsgIdAndTxt("msldem_lazyenvireadRect_mex:notScalar","Input 3 needs to be a scalar.");
    }
    if( !mxIsScalar(prhs[4]) ) {
        mexErrMsgIdAndTxt("msldem_lazyenvireadRect_mex:notScalar","Input 4 needs to be a scalar.");
    }
    if( !mxIsScalar(prhs[5]) ) {
        mexErrMsgIdAndTxt("msldem_lazyenvireadRect_mex:notScalar","Input 5 needs to be a scalar.");
    }
    
    /* -----------------------------------------------------------------
     * I/O SETUPs
     * ----------------------------------------------------------------- */
    
    /* INPUT 0 msldem_imgpath */
    msldem_imgpath = mxArrayToString(prhs[0]);
    //printf("%s\n",msldem_imgpath);
    
    /* INPUT 1 msldem_header */
    msldem_header.samples = (int) mxGetScalar(mxGetField(prhs[1],0,"samples"));
    msldem_header.lines = (int) mxGetScalar(mxGetField(prhs[1],0,"lines"));
    msldem_header.bands = (int) mxGetScalar(mxGetField(prhs[1],0,"bands"));
    msldem_header.data_type = (int) mxGetScalar(mxGetField(prhs[1],0,"data_type"));
    msldem_header.byte_order = (int) mxGetScalar(mxGetField(prhs[1],0,"byte_order"));
    msldem_header.header_offset = (int) mxGetScalar(mxGetField(prhs[1],0,"header_offset"));

    // printf("%d\n",msldem_header.samples);
    /* INPUT 2/3 msldemc_sample_offset */
    msldemc_imxy_sample_offset = (mwSize) mxGetScalar(prhs[2]);
    msldemc_imxy_line_offset = (mwSize) mxGetScalar(prhs[3]);
    
    /* INPUT 4/5 msldem rectangle size */
    msldemc_samples = (mwSize) mxGetScalar(prhs[4]);
    msldemc_lines = (mwSize) mxGetScalar(prhs[5]);
    
    plhs[0] = mxCreateNumericMatrix(msldemc_samples,msldemc_lines,mxSINGLE_CLASS,mxREAL);
    
    dem_img = (float **) mxMalloc(msldemc_lines * sizeof(*dem_img) );
    dem_img[0] = mxGetSingles(plhs[0]);
    for( j=1;j<msldemc_lines;j++){
        dem_img[j] = dem_img[j-1] + msldemc_samples;
    }

    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    msldem_lazyenvireadRect(msldem_imgpath,msldem_header,(int32_T) msldemc_imxy_sample_offset, 
            (int32_T) msldemc_imxy_line_offset, dem_img, (int32_T) msldemc_samples, (int32_T) msldemc_lines);
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(dem_img);
    
    
}
