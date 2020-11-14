/* =====================================================================
 * msldem_cleanup_numeqnb_step2p_mex.c
 * Clean up "msldem_numeqneighbors", which is the output of 
 * "msldem_get_numeqnb". The idea behind cleaning is that real 
 * low resolution pixels can be traced to a pixel that has the largest 
 * possible value in "msldem_numeqneighbors". The largest possible value
 * would be 8 with wndw_size=5, for example.
 * Exact tracing from each candidate low resolution pixel is approximately
 * implemented in this function.
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
        int8_T **msldem_numeqnb_cl2,int8_T **msldem_numeqnb_cl2p,
        int8_T eval_num)
{
    int32_T c,l;
    int32_T cc,ll;
    int32_T c_min,c_max,l_min,l_max;
    int32_T wndw_size, wleft, wright;
    int32_T S_dem, L_dem;
    int8_T flg;
    
    // wndw_size = 3;
    // wleft = (int32_T) floor(((double)(wndw_size - 1.)) / 2.);
    // wright = wndw_size - wleft - 1;
    
    S_dem = (int32_T) msldem_hdr.samples;
    L_dem = (int32_T) msldem_hdr.lines;
    
    // printf("%d \n",L_dem);
    
    // eval_num = 3; // the number for evaluated.
    // printf("%d,%d,%d\n",skip_l,msldemc_samples*s,skip_r);
    for(l=0;l<L_dem;l++){
        // printf("l=%d \n",l); 
        l_min = (l-11>0)?(l-11):0;
        l_max = (l+12<L_dem)?(l+12):L_dem;
        for(c=0;c<S_dem;c++){
            //printf("c=%d \n",c);
            if(msldem_numeqnb_cl2[c][l]>0){
                c_min = (c-11>0)?(c-11):0;
                c_max = (c+12<S_dem)?(c+12):S_dem;
                cc=c_min;
                flg = 1;
                while(flg && cc<c_max){
                    if(msldem_numeqnb_cl2[cc][l]>eval_num){
                        msldem_numeqnb_cl2p[c][l] = msldem_numeqnb_cl2[c][l];
                        flg=0;
                    }
                    cc++;
                }
                ll = l_min;
                while(flg && ll<l_max){
                    if(msldem_numeqnb_cl2[c][ll]>eval_num){
                        msldem_numeqnb_cl2p[c][l] = msldem_numeqnb_cl2[c][l];
                        flg=0;
                    }
                    ll++;
                }
                if(flg){
                    for(ll=l_min;ll<l_max;ll++){
                        for(cc=c_min;cc<c_max;cc++){
                            if(msldem_numeqnb_cl2[cc][ll]>eval_num){
                                msldem_numeqnb_cl2p[c][l] = msldem_numeqnb_cl2[c][l];
                            }
                        }
                    }
                }
            }
        }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    EnviHeader msldem_hdr;
    // bool **msldem_imFOVmask;
    int8_T **msldem_numeqneighbors;
    int8_T **msldem_numeqnb_cl2;
    int8_T **msldem_numeqnb_cl2p;
    
    int32_T si,li;
    int32_T msldem_samples, msldem_lines;
    int8_T eval_num;

    /* -----------------------------------------------------------------
     * CHECK PROPER NUMBER OF INPUTS AND OUTPUTS
     * ----------------------------------------------------------------- */
    /* the number of inputs and outputs checked */
    if(nrhs!=3) {
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
    /* -----------------------------------------------------------------
     * I/O SETUPs
     * ----------------------------------------------------------------- */
    /* INPUT 0 msldem_header */
    msldem_hdr = mxGetEnviHeader(prhs[0]);
    
    
    /* INPUT 2 msldem_numeqneighbors */
    msldem_numeqnb_cl2 = set_mxInt8Matrix(prhs[1]);
    
    /* INPUT 3 eval_num */
    eval_num = (int8_T) mxGetScalar(prhs[2]);
    //eval_num_min = (int32_T) mxGetScalar(prhs[3]);
    
    /* OUTPUT 0 cleaned msldem_numeqneighbors */
    plhs[0] = mxCreateNumericMatrix((mwSize) msldem_hdr.lines,(mwSize) msldem_hdr.samples,mxINT8_CLASS,mxREAL);
    msldem_numeqnb_cl2p = set_mxInt8Matrix(plhs[0]);
    
    // Initialize matrices
    msldem_samples = (int32_T) msldem_hdr.samples;
    msldem_lines = (int32_T) msldem_hdr.lines;
    for(si=0;si<msldem_samples;si++){
        for(li=0;li<msldem_lines;li++){
            if(msldem_numeqnb_cl2[si][li] == 0){
                msldem_numeqnb_cl2p[si][li] = 0;
            }else if(msldem_numeqnb_cl2[si][li] == -1){
                msldem_numeqnb_cl2p[si][li] = -1;
            }else{
                msldem_numeqnb_cl2p[si][li] = 0;
            }
        }
    }
    // printf("sim = %d\n",S_im);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    cleanup_numeqneighbors(msldem_hdr,msldem_numeqnb_cl2,
            msldem_numeqnb_cl2p,eval_num);
    
    /* free memories */
    mxFree(msldem_numeqnb_cl2);
    mxFree(msldem_numeqnb_cl2p);
    
}
