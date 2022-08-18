/* msldem_util.h */
/* Definitions 
 * 
 */
#ifndef MSLDEM_UTIL_H
#define MSLDEM_UTIL_H

#include <stdio.h>
#include <stdint.h>
#include "matrix.h"

/* struct for storing the offset and size of the sub image of MSLDEM */
typedef struct MSLDEMC_HEADER {
    int32_t sample_offset;
    int32_t line_offset  ;
    int32_t samples      ;
    int32_t lines        ;
} MSLDEMC_HEADER ;

/* --------------------------------------------------------------------- *
 * int mxGet_MSLDEMC_HEADER(const mxArray *pm, MSLDEMC_HEADER *msldemc_hdr)
 * 
 * Convert matlab struct of msldemc_hdr to C-format struct
 * Input Parameters
 *  const mxArray *pm
 *     pointer to the matlab struct of msldemc_hdr
 *  MSLDEMC_HEADER *msldemc_hdr
 *     pointer to C struct of msldemc_hdr
 * Return
 *  int err
 *     0 if no error, 1 otherwise
 * --------------------------------------------------------------------- */
int mxGet_MSLDEMC_HEADER(const mxArray *pm, MSLDEMC_HEADER *msldemc_hdr){
    mxArray *fldval;
    int err=0;
    
    fldval=mxGetField(pm,0,"sample_offset");
    if(fldval!=NULL){
        msldemc_hdr->sample_offset = (int32_t) mxGetScalar(fldval);
    } else {
        err=1;
        fprintf(stderr, "Error:mxGet_MSLDEMC_HEADER: Field sample_offset does not exist.");
    }
    
    fldval=mxGetField(pm,0,"line_offset");
    if(fldval!=NULL){
        msldemc_hdr->line_offset = (int32_t) mxGetScalar(fldval);
    } else {
        err=1;
        fprintf(stderr, "Error:mxGet_MSLDEMC_HEADER: Field line_offset does not exist.");
    }
    
    fldval=mxGetField(pm,0,"samples");
    if(fldval!=NULL){
        msldemc_hdr->samples = (int32_t) mxGetScalar(fldval);
    } else {
        err=1;
        fprintf(stderr, "Error:mxGet_MSLDEMC_HEADER: Field samples does not exist.");
    }
    
    fldval=mxGetField(pm,0,"lines");
    if(fldval!=NULL){
        msldemc_hdr->lines = (int32_t) mxGetScalar(fldval);
    } else {
        err=1;
        fprintf(stderr, "Error:mxGet_MSLDEMC_HEADER: Field lines does not exist.");
    }
    
    return err;
}


/* --------------------------------------------------------------------- *
 * mxArray* mxCreate_MSLDEMC_HDR_StructmxArray(
 *      mwSize sample_offset, 
 *      mwSize line_offset,
 *      mwSize samples, 
 *      mwSize lines)
 * Create Struct mxArray of msldemc_hdr
 * Return 
 *   msldemc_hdr_struct: sample_offset, line_offset, samples, lines
 * --------------------------------------------------------------------- */
mxArray* mxCreate_MSLDEMC_HDR_StructmxArray(MSLDEMC_HEADER *msldemc_hdr){
    const char *field_names[] = {"sample_offset","line_offset","samples","lines"};
    mxArray *msldemc_hdr_struct;
    mwSize dims[2] = {1, 1};
    
    msldemc_hdr_struct = mxCreateStructArray(2,dims,4,field_names);
    
    mxSetField(msldemc_hdr_struct,1,"sample_offset",mxCreateDoubleMatrix(1,1,mxREAL));
    *mxGetDoubles(mxGetField(msldemc_hdr_struct,1,"sample_offset")) = (double) msldemc_hdr->sample_offset;
    
    mxSetField(msldemc_hdr_struct,1,"line_offset",mxCreateDoubleMatrix(1,1,mxREAL));
    *mxGetDoubles(mxGetField(msldemc_hdr_struct,1,"line_offset")) = (double) msldemc_hdr->line_offset;
    
    mxSetField(msldemc_hdr_struct,1,"samples",mxCreateDoubleMatrix(1,1,mxREAL));
    *mxGetDoubles(mxGetField(msldemc_hdr_struct,1,"samples")) = (double) msldemc_hdr->samples;
    
    mxSetField(msldemc_hdr_struct,1,"samples",mxCreateDoubleMatrix(1,1,mxREAL));
    *mxGetDoubles(mxGetField(msldemc_hdr_struct,1,"lines")) = (double) msldemc_hdr->lines;
    
    return msldemc_hdr_struct;
    
}

#endif