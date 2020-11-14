/* =====================================================================
 * msldem_cleanup_numeqnb_step3p_mex.c
 * This function cleanup "msldem_numeqnb_cl2", which is the output of 
 * "msldem_cleanup_numeqnb_step2_mex". At this final step, low resolution 
 * pixels that only composed of two pixles are additionally detected. Such 
 * pixels have the value of 1 in "msldem_numeqnb" originally and not simple
 * to detect. Two criteria are used. First, such pixels are right next to 
 * low resolution pixels. Second, dem values of such pixels are distinct 
 * from the surrounding high resolution pixels. The second condition is 
 * evaluated if the dem values are more deviated than 4 sigma from the mean
 * of the surrounding high resolution pixels. 
 * 
 * 
 * INPUTS:
 * 0 msldem_imgpath        char*
 * 1 msldem_header         struct
 * 2 msldem_numeqneighbors int8 [L_dem x S_dem]
 * 3 msldem_numeqnb_cl2    int8 [L_dem x S_dem]
 * 
 * 
 * OUTPUTS:
 * 3 msldem_numeqnb_cl3    int8 [L_dem x S_dem]
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
void msldem_numeqnb_clean_3rdstep(char *msldem_imgpath, EnviHeader msldem_hdr,
        int8_T **msldem_numeqneighbors, int8_T **msldem_numeqnb_cl2,
        int8_T **msldem_numeqnb_cl3, int32_T wndw_size, int8_T eval_val)
{
    int32_T c,l,cp;
    int32_T cc,ll,ccp;
    int32_T ltu,l_elevtmp; // Line of Temporary array for update
    int32_T c_min,c_max,l_min,l_max;
    int32_T c_min3,c_max3,l_min3,l_max3;
    int32_T c_min5,c_max5,l_min5,l_max5;
    int32_T cc_min,ll_min;
    int32_T ccc,lll,cccp,lidx2;
    float **elevtmp;
    float *elevtmp_base;
    int32_T wleft, wright;
    int32_T S_dem, L_dem;
    int32_T L_dem_mmrgn;
    int32_T S_dem_pmrgn;
    size_t sz = sizeof(float);
    FILE *fid;
    float elevcl,elevccll;
    float data_ignore_value_float, data_ignore_value_float_p1;
    float m,v_min,v,vv,v_min_ccll,vv_ccll,v_max;
    int32_T N;
    int8_T flg;
    int32_T lidx,i;
    float v_min_lr;
    
    float *v_min_ccll_ar;
    
    // wndw_size = 20;
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
    
    fid = fopen(msldem_imgpath,"rb");
    
    ltu = 0;
    for(l=0;l<wright;l++){
        fread(&elevtmp[ltu][wleft],sz,S_dem,fid);
        ltu++;
        ltu = ltu % wndw_size;
    }
    
    // printf("%d \n",L_dem);
    
     v_min_ccll_ar = (float*) malloc(sizeof(float) * (size_t) wndw_size * (size_t) wndw_size);
    
    
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
        //printf("%d \n",ltu);
        
        l_min = (l-wleft>0)?(l-wleft):0;
        l_max = (l+wright+1<L_dem)?(l+wright+1):L_dem;
        
        l_min3 = (l-1>0)?(l-1):0;
        l_max3 = (l+2<L_dem)?(l+2):L_dem;
        
        l_min5 = (l-2>0)?(l-2):0;
        l_max5 = (l+3<L_dem)?(l+3):L_dem;
        
        for(c=0;c<S_dem;c++){
            cp = c+wleft;
            elevcl = elevtmp[l_elevtmp][cp];
            //printf("c=%d \n",c);
            /* if pixel (c,l) is lageled as high resolution pixels */
            if(msldem_numeqneighbors[c][l]==eval_val){
                c_min3 = (c-1>0)?(c-1):0;
                c_max3 = (c+2<S_dem)?(c+2):S_dem;
                
                /* Evaluate surrounding pixel difference */
                flg = 0;
                for(cc=c_min3;cc<c_max3;cc++){
                    if(msldem_numeqnb_cl2[cc][l]>0)
                        flg=1;
                }
                for(ll=l_min3;ll<l_max3;ll++){
                    if(msldem_numeqnb_cl2[c][ll]>0)
                        flg=1;
                }
                
                /* if the pixel contact a low resolution pixel */
                if(flg){
                    c_min = (c-wleft>0)?(c-wleft):0;
                    c_max = (c+wright+1<S_dem)?(c+wright+1):S_dem;
                    c_min5 = (c-2>0)?(c-2):0;
                    c_max5 = (c+3<S_dem)?(c+3):S_dem;
                    v_min = INFINITY;
                    v_min_lr = INFINITY;
                    for(ll=l_min5;ll<l_max5;ll++){
                        lidx = ll%wndw_size;
                        for(cc=c_min5;cc<c_max5;cc++){
                            if(cc!=c || ll!=l){
                                ccp = cc+wleft;
                                if(msldem_numeqnb_cl2[cc][ll]==0){
                                    vv = fabsf(elevtmp[lidx][ccp]-elevcl);
                                    if( (v_min > vv) && (vv > 1e-9) ){
                                        v_min = vv; cc_min = cc; ll_min = ll;
                                    }
                                }else if(msldem_numeqnb_cl2[cc][ll]>0){
                                    vv = fabsf(elevtmp[lidx][ccp]-elevcl);
                                    if( (v_min_lr > vv) && (vv > 1e-9) ){
                                        v_min_lr = vv; // cc_min = cc; ll_min = ll;
                                    }
                                }
                            }
                        }
                    }
                    
                    
                    
                    
                    /* local standard deviation */
                    m = 0;
                    N = 0;
                   
                    for(ll=l_min;ll<l_max;ll++){
                        lidx = ll%wndw_size;
                        for(cc=c_min;cc<c_max;cc++){
                            if(msldem_numeqnb_cl2[cc][ll]==0){
                                ccp = cc+wleft;
                                elevccll = elevtmp[lidx][ccp];
                                if(fabsf(elevcl-elevccll)>1e-9){
                                    v_min_ccll = INFINITY;
                                    // vv_ccll = 
                                    for(lll=l_min5;lll<l_max5;lll++){
                                        lidx2 = lll%wndw_size;
                                        for(ccc=c_min5;ccc<c_max5;ccc++){
                                            if(msldem_numeqnb_cl2[ccc][lll]==0 && (ccc!=cc || lll!=ll) && (ccc!=cc || lll!=ll)){
                                                cccp = ccc+wleft;
                                                vv_ccll = fabsf(elevtmp[lidx2][cccp]-elevccll);
                                                if( (v_min_ccll > vv_ccll) && (vv_ccll > 1e-9) ){
                                                    v_min_ccll = vv_ccll; // cc_min = cc; ll_min = ll;
                                                }
                                            }
                                        }
                                    }
                                    v_min_ccll_ar[N] = v_min_ccll;
                                    N++;
                                }
                            }
                        }
                    }
                    
                    m = 0;
                    for(i=1;i<N;i++){
                        m = m + v_min_ccll_ar[i];
                    }
                    m = m / (float) N;
                    
                    v_max = 0;
                    for(i=1;i<N;i++){
                        if(v_max < v_min_ccll_ar[i]){
                            v_max =  v_min_ccll_ar[i];
                        }
                        v = v + (v_min_ccll_ar[i]-m)*(v_min_ccll_ar[i]-m);
                    }
                    v = v / (float) N;
                    
                    // if((l==30679 && c==19760) || (l==30678 && c==19760) || (l==49260 && c==9259) || (l==49261 && c==9259)){
                    if((c==3680 && l==23459) || (l==23198 && c==3647)){
                        printf("l%d,c%d,%f,%f,%f\n",l,c,elevcl,v_min,v_min_lr);
                        printf("l%d,c%d,%f,%f,%f,%10f,%10f,%d,%d,%f\n",l,c,elevcl,v_min,v_max,(v_min-m)*(v_min-m),16*v,ll_min,cc_min,elevtmp[ll_min%wndw_size][cc_min+wleft]);
                    }
                    
                    if((v_min>0.1) && (v_min-m)*(v_min-m) > 16*v && v_min_lr<v_min){
                    // if((v_min-m)*(v_min-m) > 16*v){
                        // msldem_numeqnb_cl3[c][l] = msldem_numeqneighbors[c][l];
                        msldem_numeqnb_cl3[c][l] = 1;
                    }
                }
            }
        }
    }
    
    free(elevtmp);
    free(elevtmp_base);
    free(v_min_ccll_ar);
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
    int8_T **msldem_numeqnb_cl2;
    int8_T **msldem_numeqnb_cl3;
    
    int32_T si,li;
    int32_T msldem_samples, msldem_lines;
    
    int32_T wndw_size;
    int8_T eval_val;

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
    msldem_numeqnb_cl2 = set_mxInt8Matrix(prhs[3]);
    
    /* INPUT 4 wndw_size */
    wndw_size = (int32_T) mxGetScalar(prhs[4]);
    
    /* INPUT 4 wndw_size */
    eval_val = (int8_T) mxGetScalar(prhs[5]);

    /* OUTPUT 0 msldem imFOV */
    plhs[0] = mxCreateNumericMatrix((mwSize) msldem_hdr.lines,(mwSize) msldem_hdr.samples,mxINT8_CLASS,mxREAL);
    msldem_numeqnb_cl3 = set_mxInt8Matrix(plhs[0]);
    
    // Initialize matrices
    msldem_samples = (int32_T) msldem_hdr.samples;
    msldem_lines = (int32_T) msldem_hdr.lines;
    for(si=0;si<msldem_samples;si++){
        for(li=0;li<msldem_lines;li++)
            msldem_numeqnb_cl3[si][li] = msldem_numeqnb_cl2[si][li];
    }
    // printf("sim = %d\n",S_im);
    /* -----------------------------------------------------------------
     * CALL MAIN COMPUTATION ROUTINE
     * ----------------------------------------------------------------- */
    msldem_numeqnb_clean_3rdstep(msldem_imgpath,msldem_hdr,
            msldem_numeqneighbors,msldem_numeqnb_cl2,msldem_numeqnb_cl3,wndw_size,eval_val);
    
    /* free memories */
    mxFree(msldem_imgpath);
    mxFree(msldem_numeqneighbors);
    mxFree(msldem_numeqnb_cl2);
    mxFree(msldem_numeqnb_cl3);
    
}
