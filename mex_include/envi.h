/* envi.h */
#include <stdint.h>
#include "mex.h"


typedef enum EnviHeaderInterleave {
    BSQ,BIP,BIL
} EnviHeaderInterleave ;

typedef struct EnviHeader {
    int32_T samples;
    int32_T lines;
    int32_T bands;
    int32_T data_type;
    int32_T byte_order;
    int32_T header_offset;
    EnviHeaderInterleave interleave;
    char* file_type;
    double data_ignore_value;
} EnviHeader ;

EnviHeader mxGetEnviHeader(const mxArray *pm){
    EnviHeader msldem_hdr;
    
    if(mxGetField(pm,0,"samples")!=NULL){
        msldem_hdr.samples = (int32_T) mxGetScalar(mxGetField(pm,0,"samples"));
    }else{
        mexErrMsgIdAndTxt("envi:mexGetEnviHeader","Struct is not an envi header (sample)");
    }
    if(mxGetField(pm,0,"lines")!=NULL){
        msldem_hdr.lines = (int32_T) mxGetScalar(mxGetField(pm,0,"lines"));
    }else{
        mexErrMsgIdAndTxt("envi:mexGetEnviHeader","Struct is not an envi header");
    }
    if(mxGetField(pm,0,"bands")!=NULL){
        msldem_hdr.bands = (int32_T) mxGetScalar(mxGetField(pm,0,"bands"));
    }else{
        mexErrMsgIdAndTxt("envi:mexGetEnviHeader","Struct is not an envi header");
    }
    if(mxGetField(pm,0,"data_type")!=NULL){
        msldem_hdr.data_type = (int32_T) mxGetScalar(mxGetField(pm,0,"data_type"));
    }else{
        mexErrMsgIdAndTxt("envi:mexGetEnviHeader","Struct is not an envi header");
    }
    if(mxGetField(pm,0,"byte_order")!=NULL){
        msldem_hdr.byte_order = (int32_T) mxGetScalar(mxGetField(pm,0,"byte_order"));
    }else{
        mexErrMsgIdAndTxt("envi:mexGetEnviHeader","Struct is not an envi header");
    }
    if(mxGetField(pm,0,"header_offset")!=NULL){
        msldem_hdr.header_offset = (int32_T) mxGetScalar(mxGetField(pm,0,"header_offset"));
    }else{
        mexErrMsgIdAndTxt("envi:mexGetEnviHeader","Struct is not an envi header");
    }
    if(mxGetField(pm,0,"data_ignore_value")!=NULL){
        msldem_hdr.data_ignore_value = mxGetScalar(mxGetField(pm,0,"data_ignore_value"));
    }else{
        mexErrMsgIdAndTxt("envi:mexGetEnviHeader","Struct is not an envi header");
    }
    
    return msldem_hdr;
    
}