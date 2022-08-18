/* =====================================================================
 * get_num_equal_neighbors_mex.c
 * Evaluate if pixels in the MSL DEM image are potentially in 
 * 
 * INPUTS:
 * 0 msldem_imgpath        char*
 * 1 msldem_header         struct
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
        int8_T **msldem_numeqneighbors)
{
    int32_T c,l;
    int32_T cc,ll;
    
    float *elevlm1, *elevl ,*elevlp1;
    size_t nbytes_l;
    int32_T S_dem, L_dem;
    int32_T L_demm1;
    int32_T S_demp1;
    size_t sz = sizeof(float);
    //size_t sz1 = sizeof(size_t);
    FILE *fid;
    float elevcl;
    float data_ignore_value_float, data_ignore_value_float_p1;
    
    S_dem = (int32_T) msldem_hdr.samples;
    L_dem = (int32_T) msldem_hdr.lines;
    L_demm1 = L_dem - 1;
    S_demp1 = S_dem + 1;
    
    // printf("%d \n",L_dem);
    

    // printf("%d \n",L_dem);
    data_ignore_value_float = (float) msldem_hdr.data_ignore_value;
    data_ignore_value_float_p1 = data_ignore_value_float + 1.0;
    // printf("ignore: %f\n",data_ignore_value_float);
    
    // printf("%d \n",L_dem);

    
    // printf("size of size_t: %d\n",sz1);

    
    /* read the data */
    nbytes_l = sz * ((size_t) S_dem+2);
    elevlm1  = (float*) malloc(nbytes_l);
    elevl    = (float*) malloc(nbytes_l);
    elevlp1  = (float*) malloc(nbytes_l);
    
    // printf("%d \n",L_dem);
    
    for(c=0;c<S_dem+2;c++){
        elevlm1[c] = data_ignore_value_float;
        elevl[c]   = data_ignore_value_float;
        elevlp1[c] = data_ignore_value_float;
    }
    
    // printf("%d \n",L_dem);
    
    fid = fopen(msldem_imgpath,"rb");
    
    fread(&elevlp1[1],sz,S_dem,fid);
    for(c=1;c<S_demp1;c++){
        //if(elevlp1[c]<data_ignore_value_float)
        //    elevlp1[c] = NAN;
    }
    
    // printf("%d \n",L_dem);
    
    
    // printf("%d,%d,%d\n",skip_l,msldemc_samples*s,skip_r);
    for(l=0;l<L_dem;l++){
        // printf("l=%d \n",l);
        memcpy(elevlm1,elevl,nbytes_l);
        memcpy(elevl,elevlp1,nbytes_l);
        if(l<L_demm1){
            fread(&elevlp1[1],sz,S_dem,fid);
        } else {
            for(c=1;c<S_demp1;c++)
                elevlp1[c] = data_ignore_value_float;
        }
        
        for(c=0;c<S_dem;c++){
            elevcl = elevl[c+1];
            if(elevcl<data_ignore_value_float_p1){
                msldem_numeqneighbors[c][l] = -1;
            } else{
                /* perform only if the elevation value is valid */
                /* Evaluate the surrounding */
                for(cc=c;cc<c+3;cc++){
                    if(fabs(elevcl-elevlm1[cc]) < 1e-9)
                        msldem_numeqneighbors[c][l]++;
                    if(fabs(elevcl-elevlp1[cc]) < 1e-9)
                        msldem_numeqneighbors[c][l]++;
                }
                if(fabs(elevcl-elevl[c]) < 1e-9)
                    msldem_numeqneighbors[c][l]++;
                if(fabs(elevcl-elevl[c+2]) < 1e-9)
                    msldem_numeqneighbors[c][l]++;
            }
            
        }


    }
    free(elevlm1);
    free(elevl);
    free(elevlp1);
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
    msldem_hdr.samples = (int32_T) mxGetScalar(mxGetField(prhs[1],0,"samples"));
    msldem_hdr.lines = (int32_T) mxGetScalar(mxGetField(prhs[1],0,"lines"));
    msldem_hdr.bands = (int32_T) mxGetScalar(mxGetField(prhs[1],0,"bands"));
    msldem_hdr.data_type = (int32_T) mxGetScalar(mxGetField(prhs[1],0,"data_type"));
    msldem_hdr.byte_order = (int32_T) mxGetScalar(mxGetField(prhs[1],0,"byte_order"));
    msldem_hdr.header_offset = (int32_T) mxGetScalar(mxGetField(prhs[1],0,"header_offset"));
    msldem_hdr.data_ignore_value = mxGetScalar(mxGetField(prhs[1],0,"data_ignore_value"));

    /* OUTPUT 0 msldem imFOV */
    plhs[0] = mxCreateNumericMatrix((mwSize) msldem_hdr.lines,(mwSize) msldem_hdr.samples,mxINT8_CLASS,mxREAL);
    msldem_numeqneighbors = set_mxInt8Matrix(plhs[0]);
    
    // Initialize matrices
    msldem_samples = (int32_T) msldem_hdr.samples;
    msldem_lines = (int32_T) msldem_hdr.lines;
    for(si=0;si<msldem_samples;si++){
        for(li=0;li<msldem_lines;li++){
            // msldem_imFOVmask[si][li] = false;
            msldem_numeqneighbors[si][li] = 0;
        }
    }
    // printf("sim = %d\n",S_im);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    get_num_equal_neighbors(msldem_imgpath,msldem_hdr,msldem_numeqneighbors);
    
    /* free memories */
    mxFree(msldem_imgpath);
    // mxFree(msldem_imFOVmask);
    mxFree(msldem_numeqneighbors);
    
}
