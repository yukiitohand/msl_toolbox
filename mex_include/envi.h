/* envi.h */

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