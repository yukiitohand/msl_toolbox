/* =====================================================================
 * msldem_average_20x20_lowrespixels_mex.c
 * Apply 20x20 average filter on the low resolution pixels of msldem data.
 * 
 * INPUTS:
 * 0 msldem_imgpath        char*
 * 1 msldem_header         struct
 * 2 msldem_numeqneighbors int8 [L_dem x S_dem]
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
void msldem_average_20x20_lowrespixels(char *msldem_imgpath, EnviHeader msldem_hdr,
        int8_T **msldem_numeqneighbors,float **msldem_mod)
{
    int32_T c,l;
    int32_T cc,ll;
    int32_T ltu,l_elevtmp; // Line of Temporary array for update
    float **elevtmp;
    float *elevtmp_base;
    int32_T ave_wndw_size, wleft, wright;
    int32_T S_dem, L_dem;
    int32_T L_dem_mmrgn;
    int32_T S_dem_pmrgn;
    size_t sz = sizeof(float);
    //size_t sz1 = sizeof(size_t);
    FILE *fid;
    float data_ignore_value_float, data_ignore_value_float_p1;
    float val;
    int32_T N;
    
    ave_wndw_size = 20;
    wleft = (int32_T) floor(((double)(ave_wndw_size - 1.)) / 2.);
    wright = ave_wndw_size - wleft - 1;
    
    S_dem = (int32_T) msldem_hdr.samples;
    L_dem = (int32_T) msldem_hdr.lines;
    L_dem_mmrgn = L_dem - wright;
    S_dem_pmrgn = S_dem + wleft + wright;

    // printf("%d \n",L_dem);
    

    // printf("%d \n",L_dem);
    data_ignore_value_float = (float) msldem_hdr.data_ignore_value;
    data_ignore_value_float_p1 = data_ignore_value_float + 1.0;
    // printf("ignore: %f\n",data_ignore_value_float);
    
    printf("%d \n",L_dem);

    
    // printf("size of size_t: %d\n",sz1);
    
    // create 20 x (S_dem + 20 - 1) temporary array
    elevtmp = (float**) malloc(sizeof(float*) * (size_t) ave_wndw_size);
    elevtmp_base = (float*) malloc(sizeof(float) * (size_t) (ave_wndw_size * S_dem_pmrgn));
    elevtmp[0] = &elevtmp_base[0];
    for(l=1;l<ave_wndw_size;l++){
        elevtmp[l] = elevtmp[l-1] + S_dem_pmrgn;
    }
    
    // initialize the temporary array
    for(l=0;l<ave_wndw_size;l++){
        for(c=0;c<S_dem_pmrgn;c++)
            elevtmp[l][c] = data_ignore_value_float;
    }

    
    
    fid = fopen(msldem_imgpath,"rb");
    
    ltu = 0;
    for(l=0;l<wright;l++){
        fread(&elevtmp[ltu][wleft],sz,S_dem,fid);
        ltu++;
        ltu = ltu % ave_wndw_size;
    }
    
    printf("%d \n",L_dem);
    
    
    // printf("%d,%d,%d\n",skip_l,msldemc_samples*s,skip_r);
    for(l=0;l<L_dem;l++){
        // printf("l=%d \n",l);
        l_elevtmp = l % ave_wndw_size;
        
        if(l<L_dem_mmrgn){
            fread(&elevtmp[ltu][wleft],sz,S_dem,fid);
            //printf("%d \n",l);
        } else {
            for(c=0;c<S_dem_pmrgn;c++)
                elevtmp[ltu][c] = data_ignore_value_float;
        }
        ltu++;
        ltu = ltu % ave_wndw_size;
        //printf("%d \n",ltu);
        
        for(c=0;c<S_dem;c++){
            //printf("c=%d \n",c);
            if(msldem_numeqneighbors[c][l]==-1){
                //printf("%d \n",msldem_numeqneighbors[c][l]);
                msldem_mod[c][l] = data_ignore_value_float;
            } else if(msldem_numeqneighbors[c][l]==0) {
                //printf("%d \n",msldem_numeqneighbors[c][l]);
                msldem_mod[c][l] = elevtmp[l_elevtmp][c+wleft];
            } else if(msldem_numeqneighbors[c][l]>0){
                /* Calculate the average value  */
                //printf("%d \n",msldem_numeqneighbors[c][l]);
                val = 0; N = 0;
                for(ll=0;ll<ave_wndw_size;ll++){
                    for(cc=c;cc<c+ave_wndw_size;cc++){
                        if(!(elevtmp[ll][cc] < data_ignore_value_float_p1)){
                            val += elevtmp[ll][cc];
                            N++;
                        }    
                    }
                }
                val = val / (float) N;
                msldem_mod[c][l] = val;
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
    float **msldem_mod;
    
    int32_T si,li;
    int32_T msldem_samples, msldem_lines;

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

    /* OUTPUT 0 msldem imFOV */
    plhs[0] = mxCreateNumericMatrix((mwSize) msldem_hdr.lines,(mwSize) msldem_hdr.samples,mxSINGLE_CLASS,mxREAL);
    msldem_mod = set_mxSingleMatrix(plhs[0]);
    
    // Initialize matrices
    msldem_samples = (int32_T) msldem_hdr.samples;
    msldem_lines = (int32_T) msldem_hdr.lines;
    for(si=0;si<msldem_samples;si++){
        for(li=0;li<msldem_lines;li++){
            // msldem_imFOVmask[si][li] = false;
            msldem_mod[si][li] = NAN;
        }
    }
    // printf("sim = %d\n",S_im);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    msldem_average_20x20_lowrespixels(msldem_imgpath,msldem_hdr,msldem_numeqneighbors,msldem_mod);
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldem_numeqneighbors);
    mxFree(msldem_mod);
    
}
