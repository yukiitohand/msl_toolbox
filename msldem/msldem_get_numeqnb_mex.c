/* =====================================================================
 * msldem_get_numeqnb_mex.c
 * Count the number of surrounding pixels that has the same pixel value.
 * The surrouding pixels are ones that share either the column or the row
 * inside the window whose size is specified by an input.
 *        _
 *       |_|
 *    _ _|_|_ _
 *   |_|_|_|_|_|
 *       |_|
 *       |_|
 *   <--------->
 *   (wndw_size): size of the window (five in this example)
 * 
 * The pixels that have invalid values (data_ignore_value) are filled with 
 * -1.
 *      
 * INPUTS:
 * 0 msldem_imgpath        char*   : path to the image file
 * 1 msldem_header         struct  : envi header struct
 * 2 wndw_size             scalar  : size of the window
 * 
 * 
 * OUTPUTS:
 * 0 msldem_numeqneighbors    int8 [L_dem x S_dem]
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
void get_num_equal_neighbors(char *msldem_imgpath, EnviHeader msldem_hdr,
        int8_T **msldem_numeqneighbors,int32_T wndw_size)
{
    int32_T c,l,cp;
    int32_T cc,ll;
    int32_T ltu,l_elevtmp; // Line of Temporary array for update
    float **elevtmp;
    float *elevtmp_base;
    int32_T wleft, wright;
    int32_T S_dem, L_dem;
    int32_T L_dem_mmrgn;
    int32_T S_dem_pmrgn;
    size_t sz = sizeof(float);
    FILE *fid;
    float elevcl;
    float data_ignore_value_float, data_ignore_value_float_p1;
    
    
    // ave_wndw_size = 5;
    wleft = (int32_T) floor(((double)(wndw_size - 1.)) / 2.);
    wright = wndw_size - wleft - 1;
    
    S_dem = (int32_T) msldem_hdr.samples;
    L_dem = (int32_T) msldem_hdr.lines;
    L_dem_mmrgn = L_dem - wright;
    S_dem_pmrgn = S_dem + wleft + wright;

    
    // printf("%d \n",L_dem);
    

    // printf("%d \n",L_dem);
    data_ignore_value_float = (float) msldem_hdr.data_ignore_value;
    data_ignore_value_float_p1 = data_ignore_value_float + 1.0;
    // printf("ignore: %f\n",data_ignore_value_float);
    
    // printf("%d \n",L_dem);

    
    // printf("size of size_t: %d\n",sz1);

    
    // create 20 x (S_dem + 20 - 1) temporary array
    elevtmp = (float**) malloc(sizeof(float*) * (size_t) wndw_size);
    elevtmp_base = (float*) malloc(sizeof(float) * (size_t) (wndw_size * S_dem_pmrgn));
    elevtmp[0] = &elevtmp_base[0];
    for(l=1;l<wndw_size;l++){
        elevtmp[l] = elevtmp[l-1] + S_dem_pmrgn;
    }
    
    // initialize the temporary array
    for(l=0;l<wndw_size;l++){
        for(c=0;c<S_dem_pmrgn;c++)
            elevtmp[l][c] = data_ignore_value_float;
    }
    
    // printf("%d \n",L_dem);
    
    fid = fopen(msldem_imgpath,"rb");
    
    ltu = 0;
    for(l=0;l<wright;l++){
        fread(&elevtmp[ltu][wleft],sz,S_dem,fid);
        ltu++;
        ltu = ltu % wndw_size;
    }
    
    // printf("%d \n",L_dem);
    
    
    // printf("%d,%d,%d\n",skip_l,msldemc_samples*s,skip_r);
    for(l=0;l<L_dem;l++){
        // printf("l=%d \n",l);
        l_elevtmp = l % wndw_size;
        
       if(l<L_dem_mmrgn){
            fread(&elevtmp[ltu][wleft],sz,S_dem,fid);
            //printf("%d \n",l);
        } else {
            for(c=0;c<S_dem_pmrgn;c++)
                elevtmp[ltu][c] = data_ignore_value_float;
        }
        ltu++;
        ltu = ltu % wndw_size;
        
        for(c=0;c<S_dem;c++){
            cp = c+wleft;
            elevcl = elevtmp[l_elevtmp][cp];
            if(elevcl<data_ignore_value_float_p1){
                msldem_numeqneighbors[c][l] = -1;
            } else{
                /* perform only if the elevation value is valid */
                /* Evaluate the surrounding */
                for(cc=c;cc<c+wndw_size;cc++){
                    if( (cc!=cp) && (fabsf(elevcl-elevtmp[l_elevtmp][cc]) < 1e-9))
                        msldem_numeqneighbors[c][l]++;
                }
                for(ll=0;ll<wndw_size;ll++){
                    if ((ll!=l_elevtmp) && (fabsf(elevcl-elevtmp[ll][cp]) < 1e-9))
                        msldem_numeqneighbors[c][l]++;
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
    
    int32_T si,li;
    int32_T msldem_samples, msldem_lines;
    int32_T wndw_size;

    /* -----------------------------------------------------------------
     * CHECK PROPER NUMBER OF INPUTS AND OUTPUTS
     * ----------------------------------------------------------------- */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("get_num_equal_neighbors_mex_v2:nrhs","Three inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("get_num_equal_neighbors_mex_v2:nlhs","One output required.");
    }
    
    /* Input Check */
    if( !mxIsChar(prhs[0]) ) {
        mexErrMsgIdAndTxt("cleanup_numeqneighbors:notChar","Input 9 needs to be a string or characters.");
    }
    if( !mxIsStruct(prhs[1]) ) {
        mexErrMsgIdAndTxt("cleanup_numeqneighbors:notStruct","Input 1 needs to be a header struct.");
    }
    if( !mxIsScalar(prhs[2]) ) {
        mexErrMsgIdAndTxt("cleanup_numeqneighbors:notScalar","Input 2 needs to be a scalar.");
    }
    /* -----------------------------------------------------------------
     * I/O SETUPs
     * ----------------------------------------------------------------- */
    
    /* INPUT 0 msldem_imgpath */
    msldem_imgpath = mxArrayToString(prhs[0]);
    
    /* INPUT 1 msldem_header */
    msldem_hdr = mxGetEnviHeader(prhs[1]);
    
    /* INPUT 2  */
    wndw_size = (int32_T) mxGetScalar(prhs[2]);

    /* OUTPUT 0 msldem imFOV */
    plhs[0] = mxCreateNumericMatrix((mwSize) msldem_hdr.lines,(mwSize) msldem_hdr.samples,mxINT8_CLASS,mxREAL);
    msldem_numeqneighbors = set_mxInt8Matrix(plhs[0]);
    
    // Initialize matrices
    msldem_samples = (int32_T) msldem_hdr.samples;
    msldem_lines = (int32_T) msldem_hdr.lines;
    for(si=0;si<msldem_samples;si++){
        for(li=0;li<msldem_lines;li++){
            msldem_numeqneighbors[si][li] = 0;
        }
    }
    // printf("sim = %d\n",S_im);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    get_num_equal_neighbors(msldem_imgpath,msldem_hdr,msldem_numeqneighbors,
            wndw_size);
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldem_numeqneighbors);
    
}
