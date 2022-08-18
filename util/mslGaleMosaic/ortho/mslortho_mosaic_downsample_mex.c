/* =====================================================================
 * mslortho_mosaic_downsample_mex.c
 * Read the specified rectangle region of the MOLA MEGDR image data. The 
 * rectangle part of the image:
 * MOLAMEGDR[smpl_offset:smpl_offset+samples,
 *           line_offset:line_offset+lines]
 * is read.
 * 
 * INPUTS:
 * 0 imgpath       char*
 * 1 header        struct
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
#include "math.h"
#include "matrix.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "envi.h"
#include "mex_create_array.h"


/* main computation routine */
void mslortho_mosaic_downsample(char *imgpath, EnviHeader hdr, 
        char *out_imgpath, int32_T ave_wndw_sz)
{
    int32_T samples, lines, samples_ave;
    int32_T i,j,k;
    uint8_T *buf_in;
    uint8_T *buf_out;
    double *buf_ave, *num_ave;
    size_t s_in=sizeof(uint8_T);
    size_t s_ave = sizeof(double);
    FILE *fid_in, *fid_out;
    size_t ncpy_in, ncpy_ave;
    
    
    fid_out = fopen(out_imgpath,"wb");
    fid_in  = fopen(imgpath,"rb");
    
    samples = hdr.samples;
    lines   = hdr.lines;
    
    /* read the data */
    ncpy_in = s_in * (size_t) samples;
    buf_in = (uint8_T*) malloc(ncpy_in);
    
    if(samples % ave_wndw_sz == 0)
        samples_ave = samples/ave_wndw_sz;
    else
        samples_ave = samples/ave_wndw_sz + 1;
    
    ncpy_ave = s_ave * (size_t) samples_ave;
    buf_ave = (double*) malloc(ncpy_ave);
    num_ave = (double*) malloc(ncpy_ave);
    buf_out = (uint8_T*) malloc(s_in * (size_t) samples_ave);
    
    for(k=0;k<samples_ave;k++){
        buf_ave[k] = 0;
        num_ave[k] = 0;
    }
    
    // printf("%d,%d,%d\n",skip_l,msldemc_samples*s,skip_r);
    for(i=0;i<lines;i++){
        if(i % 1000 == 0){
           printf("i=%d/%d\n",i,lines);
        }
        fread(buf_in,s_in,samples,fid_in);
        for(j=0;j<samples;j++){
            k = j / ave_wndw_sz;
            buf_ave[k] += (double) buf_in[j];
            num_ave[k]++;
        }
        if( ( (i+1) % ave_wndw_sz ) == 0){
            // printf("i=%d\n",i);
            for(k=0;k<samples_ave;k++){
                buf_ave[k] /= (double) num_ave[k];
                buf_out[k] = (uint8_T) floor(buf_ave[k]+0.5);
                buf_ave[k] = 0;
                num_ave[k] = 0;
            }
            // printf("buf_out[%d] = %d\n",k,buf_out[k]);
            fwrite(buf_out, s_in, samples_ave, fid_out);
        }
    }
    free(buf_in);
    fclose(fid_in);
    free(buf_ave);
    free(num_ave);
    free(buf_out);
    fclose(fid_out);
}



/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *imgpath;
    EnviHeader hdr;
    char *out_imgpath;
    int32_T ave_wndw_sz;

    /* -----------------------------------------------------------------
     * CHECK PROPER NUMBER OF INPUTS AND OUTPUTS
     * ----------------------------------------------------------------- */
//     if(nrhs!=6) {
//         mexErrMsgIdAndTxt("mola_megdr_lazyenvireadRectInt16_mex:nrhs","Six inputs required.");
//     }
//     if(nlhs!=1) {
//         mexErrMsgIdAndTxt("mola_megdr_lazyenvireadRectInt16_mex:nlhs","Five outputs required.");
//     }
//     /* make sure the first input argument is scalar */
//     if( !mxIsChar(prhs[0]) ) {
//         mexErrMsgIdAndTxt("mola_megdr_lazyenvireadRectInt16_mex:notChar","Input 0 needs to be a string.");
//     }
//     if( !mxIsStruct(prhs[1]) ) {
//         mexErrMsgIdAndTxt("mola_megdr_lazyenvireadRectInt16_mex:notStruct","Input 1 needs to be a struct.");
//     }
//     if( !mxIsScalar(prhs[2]) ) {
//         mexErrMsgIdAndTxt("mola_megdr_lazyenvireadRectInt16_mex:notScalar","Input 2 needs to be a scalar.");
//     }
//     if( !mxIsScalar(prhs[3]) ) {
//         mexErrMsgIdAndTxt("mola_megdr_lazyenvireadRectInt16_mex:notScalar","Input 3 needs to be a scalar.");
//     }
//     if( !mxIsScalar(prhs[4]) ) {
//         mexErrMsgIdAndTxt("mola_megdr_lazyenvireadRectInt16_mex:notScalar","Input 4 needs to be a scalar.");
//     }
//     if( !mxIsScalar(prhs[5]) ) {
//         mexErrMsgIdAndTxt("mola_megdr_lazyenvireadRectInt16_mex:notScalar","Input 5 needs to be a scalar.");
//     }
    
    /* -----------------------------------------------------------------
     * I/O SETUPs
     * ----------------------------------------------------------------- */
    
    /* INPUT 0 msldem_imgpath */
    imgpath = mxArrayToString(prhs[0]);
    //printf("%s\n",msldem_imgpath);
    
    /* INPUT 1 msldem_header */
    hdr = mxGetEnviHeader(prhs[1]);
    
    /* INPUT 2 output image path */
    out_imgpath = mxArrayToString(prhs[2]);
    
    /* INPUT 3 image_ave_wndw_size */
    ave_wndw_sz = (int32_T) mxGetScalar(prhs[3]);
    
    
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    mslortho_mosaic_downsample(imgpath,hdr,out_imgpath,ave_wndw_sz);
    
    /* free memories */
    mxFree(imgpath);
    mxFree(out_imgpath);
    
    
}
