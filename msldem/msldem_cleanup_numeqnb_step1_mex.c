/* =====================================================================
 * msldem_cleanup_numeqnb_step1_mex.c
 * Clean up "msldem_numeqneighbors", which is the output of 
 * "msldem_get_numeqnb". The idea behind cleaning is that real 
 * low resolution pixels can be traced to a pixel that has the largest 
 * possible value in "msldem_numeqneighbors". The largest possible value
 * would be 8 with wndw_size=5, for example.
 * Exact tracing from each candidate low resolution pixel is approximately
 * implemented in this function.
 * We assume that the pixel is automatically considered as a low resolution
 * pixel if the number of neighbors that have the same value is greater 
 * than a certain value ("eval_num"). A typical of such value would be 3 or
 * 4. Then it is evaluated that a pixel is directly connected to the pixel 
 * that is traced to a previously determined low resolution pixel.
 * We perform this iteratively. We start evaluating the pixel with the 
 * value right below the threshold, then test if there is any neighbor that
 * has value larger than the threshold. If there exists any, the pixel is 
 * determined to be a lower resolution pixl and the value retaine. If not,
 * the pixel value is determined to be a higher reoslution one and replaced
 * with zero. Next we decrement the value to be assessed by one and test if 
 * there exists any neighbor that has a larger value. If does, the pixel is 
 * considered to be a lower resolution one. If not, it is considered to be 
 * a high resolution one. Repeat until the value gets to "eval_num_min".
 * The neighboring pixels are the four pixels below.
 *      _
 *    _|_|_ 
 *   |_|_|_|
 *     |_|
 * 
 * INPUTS:
 * 0 msldem_header         struct, envi header
 * 1 msldem_numeqneighbors int8 [L_dem x S_dem] 
 *    - output of "get_num_equal_neighbors_mex_v2"
 * 2 eval_num              scalar (integer)
 * 3 eval_num_min          scalar (integer)
 * 
 * OUTPUTS:
 * 0 msldem_numeqnb_clean  int8 [L_dem x S_dem]
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
void cleanup_numeqneighbors(EnviHeader msldem_hdr,
        int8_T **msldem_numeqneighbors,int8_T **msldem_numeqnb_clean,
        int8_T eval_num, int8_T eval_num_min)
{
    int32_T c,l;
    int32_T cc,ll;
    int32_T c_min,c_max,l_min,l_max;
    int32_T wndw_size, wleft, wright;
    int32_T S_dem, L_dem;
    
    wndw_size = 3;
    wleft = (int32_T) floor(((double)(wndw_size - 1.)) / 2.);
    wright = wndw_size - wleft - 1;
    
    S_dem = (int32_T) msldem_hdr.samples;
    L_dem = (int32_T) msldem_hdr.lines;
    
    // printf("%d \n",L_dem);
    
    // eval_num = 3; // the number for evaluated.
    // printf("%d,%d,%d\n",skip_l,msldemc_samples*s,skip_r);
    while(eval_num>(eval_num_min-1)){
        for(l=0;l<L_dem;l++){
            // printf("l=%d \n",l); 
            l_min = (l-1>0)?(l-1):0;
            l_max = (l+2<L_dem)?(l+2):L_dem;
            for(c=0;c<S_dem;c++){
                //printf("c=%d \n",c);
                if(msldem_numeqneighbors[c][l]==eval_num){
                    c_min = (c-1>0)?(c-1):0;
                    c_max = (c+2<S_dem)?(c+2):S_dem;
                    for(cc=c_min;cc<c_max;cc++){
                        if(msldem_numeqnb_clean[cc][l]>eval_num)
                            msldem_numeqnb_clean[c][l] = msldem_numeqneighbors[c][l];
                    }
                    for(ll=l_min;ll<l_max;ll++){
                        if(msldem_numeqnb_clean[c][ll]>eval_num)
                            msldem_numeqnb_clean[c][l] = msldem_numeqneighbors[c][l];
                    }
                }
            }
        }
        --eval_num;
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    EnviHeader msldem_hdr;
    // bool **msldem_imFOVmask;
    int8_T **msldem_numeqneighbors;
    int8_T **msldem_numeqnb_clean;
    
    int32_T si,li;
    int32_T msldem_samples, msldem_lines;
    int8_T eval_num, eval_num_min;

    /* -----------------------------------------------------------------
     * CHECK PROPER NUMBER OF INPUTS AND OUTPUTS
     * ----------------------------------------------------------------- */
    /* the number of inputs and outputs checked */
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("cleanup_numeqneighbors:nrhs","Three inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("cleanup_numeqneighbors:nlhs","One output required.");
    }
    /* Input Check */
    if( !mxIsStruct(prhs[0]) ) {
        mexErrMsgIdAndTxt("cleanup_numeqneighbors:notStruct","Input 0 needs to be a header struct.");
    }
    if( !mxIsClass(prhs[1],"int8") ) {
        mexErrMsgIdAndTxt("cleanup_numeqneighbors:notArray","Input 1 needs to be an Int8 array.");
    }
    if( !mxIsScalar(prhs[2]) ) {
        mexErrMsgIdAndTxt("cleanup_numeqneighbors:notScalar","Input 2 needs to be a scalar.");
    }
    if( !mxIsScalar(prhs[3]) ) {
        mexErrMsgIdAndTxt("cleanup_numeqneighbors:notScalar","Input 3 needs to be a scalar.");
    }
    /* -----------------------------------------------------------------
     * I/O SETUPs
     * ----------------------------------------------------------------- */
    /* INPUT 0 msldem_header */
    msldem_hdr = mxGetEnviHeader(prhs[0]);
    
    /* INPUT 1 msldem_numeqneighbors */
    msldem_numeqneighbors = set_mxInt8Matrix(prhs[1]);
    
    /* INPUT 2 eval_num */
    eval_num = (int32_T) mxGetScalar(prhs[2]);
    eval_num_min = (int32_T) mxGetScalar(prhs[3]);
    
    /* OUTPUT 0 cleaned msldem_numeqneighbors */
    plhs[0] = mxCreateNumericMatrix((mwSize) msldem_hdr.lines,(mwSize) msldem_hdr.samples,mxINT8_CLASS,mxREAL);
    msldem_numeqnb_clean = set_mxInt8Matrix(plhs[0]);
    
    // Initialize matrices
    msldem_samples = (int32_T) msldem_hdr.samples;
    msldem_lines = (int32_T) msldem_hdr.lines;
    for(si=0;si<msldem_samples;si++){
        for(li=0;li<msldem_lines;li++){
            // msldem_imFOVmask[si][li] = false;
            if(msldem_numeqneighbors[si][li] > eval_num){
                msldem_numeqnb_clean[si][li] = msldem_numeqneighbors[si][li];
            }else if(msldem_numeqneighbors[si][li] == 0){
                msldem_numeqnb_clean[si][li] = 0;
            }else if(msldem_numeqneighbors[si][li] == -1){
                msldem_numeqnb_clean[si][li] = -1;
            }
        }
    }
    // printf("sim = %d\n",S_im);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    cleanup_numeqneighbors(msldem_hdr,msldem_numeqneighbors,
            msldem_numeqnb_clean,eval_num,eval_num_min);
    
    /* free memories */
    mxFree(msldem_numeqneighbors);
    mxFree(msldem_numeqnb_clean);
    
}
