/* =====================================================================
 * msldem_cleanup_numeqnb_step2_mex.c
 * cleanup "msldem_numeqnb_clean", which is the output of 
 * "msldem_cleanup_numeqnb_step1_mex". This function complements low 
 * resolution pixels that are not detected in the first step. Such pixels 
 * lie around the borders of the low resolution pixels where 
 * "msldem_numeqnb" may have small values. By evaluating the value of the 
 * pixels that are marked as high resolution in "msldem_numeqnb_clean"
 * is equal to the value of the surrounding low resolution pixels. If that 
 * is equal, then the pixels are marked as low resolution pixels.
 * Surrounding pixels can be defined as below:
 *        _
 *       |_|
 *    _ _|_|_ _
 *   |_|_|_|_|_|
 *       |_|
 *       |_|
 *   <--------->
 *   (wndw_size): size of the window (five in this example)
 * 
 * 
 * INPUTS:
 * 0 msldem_imgpath        char*
 * 1 msldem_header         struct
 * 2 msldem_numeqnb        int8 [L_dem x S_dem]
 * 3 msldem_numeqnb_clean  int8 [L_dem x S_dem]
 * 4 wndw_size             scalar
 * 
 * 
 * OUTPUTS:
 * 3 msldem_mod            float [L_dem x S_dem]
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
void msldem_numeqnb_clean_2ndstep(char *msldem_imgpath, EnviHeader msldem_hdr,
        int8_T **msldem_numeqneighbors, int8_T **msldem_numeqnb_clean,
        int8_T **msldem_numeqnb_cl2, int32_T wndw_size)
{
    int32_T c,l;
    int32_T cc,ll;
    int32_T ltu,l_elevtmp; // Line of Temporary array for update
    int32_T c_min,c_max,l_min,l_max;
    float **elevtmp;
    float *elevtmp_base;
    int32_T wleft, wright;
    int32_T S_dem, L_dem;
    int32_T L_dem_mmrgn;
    // int32_T S_dem_pmrgn;
    size_t sz = sizeof(float);
    FILE *fid;
    float elevcl;
    float data_ignore_value_float, data_ignore_value_float_p1;
    
    // wndw_size = 20;
    wleft = (int32_T) floor(((double)(wndw_size - 1.)) / 2.);
    wright = wndw_size - wleft - 1;
    
    S_dem = (int32_T) msldem_hdr.samples;
    L_dem = (int32_T) msldem_hdr.lines;
    L_dem_mmrgn = L_dem - wright;
    // S_dem_pmrgn = S_dem + wleft + wright;

    // printf("%d \n",L_dem);
    

    // printf("%d \n",L_dem);
    data_ignore_value_float = (float) msldem_hdr.data_ignore_value;
    data_ignore_value_float_p1 = data_ignore_value_float + 1.0;
    // printf("ignore: %f\n",data_ignore_value_float);
    
    // printf("%d \n",L_dem);

    
    // printf("size of size_t: %d\n",sz1);
    
    // create 20 x (S_dem + 20 - 1) temporary array
    elevtmp = (float**) malloc(sizeof(float*) * (size_t) wndw_size);
    elevtmp_base = (float*) malloc(sizeof(float) * (size_t) (wndw_size * S_dem));
    elevtmp[0] = &elevtmp_base[0];
    for(l=1;l<wndw_size;l++){
        elevtmp[l] = elevtmp[l-1] + S_dem;
    }
    
    // initialize the temporary array
    for(l=0;l<wndw_size;l++){
        for(c=0;c<S_dem;c++)
            elevtmp[l][c] = data_ignore_value_float;
    }
    
    fid = fopen(msldem_imgpath,"rb");
    
    ltu = 0;
    for(l=0;l<wright;l++){
        fread(elevtmp[ltu],sz,S_dem,fid);
        ltu++;
        ltu = ltu % wndw_size;
    }
    
    // printf("%d \n",L_dem);
    
    
    // printf("%d,%d,%d\n",skip_l,msldemc_samples*s,skip_r);
    for(l=0;l<L_dem;l++){
        // printf("l=%d \n",l);
        l_elevtmp = l % wndw_size;
        
        if(l<L_dem_mmrgn){
            fread(elevtmp[ltu],sz,S_dem,fid);
            //printf("%d \n",l);
        } else {
            for(c=0;c<S_dem;c++)
                elevtmp[ltu][c] = data_ignore_value_float;
        }
        ltu++;
        ltu = ltu % wndw_size;
        //printf("%d \n",ltu);
        
        l_min = (l-wleft>0)?(l-wleft):0;
        l_max = (l+wright+1<L_dem)?(l+wright+1):L_dem;
        
        for(c=0;c<S_dem;c++){
            // cp = c+wleft;
            elevcl = elevtmp[l_elevtmp][c];
            //printf("c=%d \n",c);
            /* if pixel (c,l) is lageled as high resolution pixels */
            if(msldem_numeqnb_clean[c][l]==0){
                c_min = (c-wleft>0)?(c-wleft):0;
                c_max = (c+wright+1<S_dem)?(c+wright+1):S_dem;
                for(cc=c_min;cc<c_max;cc++){
                    if( (cc!=c) && msldem_numeqnb_clean[cc][l]>0 
                            && (fabsf(elevcl-elevtmp[l_elevtmp][cc]) < 1e-9))
                        msldem_numeqnb_cl2[c][l] = msldem_numeqneighbors[c][l];
                }
                for(ll=l_min;ll<l_max;ll++){
                    if(ll!=l && msldem_numeqnb_clean[c][ll]>0
                            && (fabsf(elevcl-elevtmp[ll%wndw_size][c]) < 1e-9) )
                        msldem_numeqnb_cl2[c][l] = msldem_numeqneighbors[c][l];
                }
            }
        }
    }
    
    free(elevtmp);
    free(elevtmp_base);
    fclose(fid);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *msldem_imgpath;
    EnviHeader msldem_hdr;
    // bool **msldem_imFOVmask;
    int8_T **msldem_numeqneighbors;
    int8_T **msldem_numeqnb_clean;
    int8_T **msldem_numeqnb_cl2;
    
    int32_T si,li;
    int32_T msldem_samples, msldem_lines;
    
    int32_T wndw_size;

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
    
    /* INPUT 2 msldem_numeqneighbors */
    msldem_numeqneighbors = set_mxInt8Matrix(prhs[2]);
    
    /* INPUT 3 msldem_numeqneighbors */
    msldem_numeqnb_clean = set_mxInt8Matrix(prhs[3]);
    
    /* INPUT 4 wndw_size */
    wndw_size = (int32_T) mxGetScalar(prhs[4]);

    /* OUTPUT 0 msldem imFOV */
    plhs[0] = mxCreateNumericMatrix((mwSize) msldem_hdr.lines,(mwSize) msldem_hdr.samples,mxINT8_CLASS,mxREAL);
    msldem_numeqnb_cl2 = set_mxInt8Matrix(plhs[0]);
    
    // Initialize matrices
    msldem_samples = (int32_T) msldem_hdr.samples;
    msldem_lines = (int32_T) msldem_hdr.lines;
    for(si=0;si<msldem_samples;si++){
        for(li=0;li<msldem_lines;li++){
            // msldem_imFOVmask[si][li] = false;
            if(msldem_numeqneighbors[si][li] > 0){
                msldem_numeqnb_cl2[si][li] = msldem_numeqnb_clean[si][li];
            }else if(msldem_numeqneighbors[si][li] == 0){
                msldem_numeqnb_cl2[si][li] = 0;
            }else if(msldem_numeqneighbors[si][li] == -1){
                msldem_numeqnb_cl2[si][li] = -1;
            }
        }
    }
    // printf("sim = %d\n",S_im);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    msldem_numeqnb_clean_2ndstep(msldem_imgpath,msldem_hdr,
            msldem_numeqneighbors,msldem_numeqnb_clean,msldem_numeqnb_cl2,wndw_size);
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldem_numeqneighbors);
    mxFree(msldem_numeqnb_clean);
    mxFree(msldem_numeqnb_cl2);
    
}
