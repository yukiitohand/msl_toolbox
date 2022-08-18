function [mst_filter_info] = mastcam_get_filter_info(msteye)
mst_filter_info = [];
switch upper(msteye)
    case 'L'
        mst_filter_info.L0 = [];
        mst_filter_info.L0.wv_ctr = [640 554 495];
        mst_filter_info.L0.HWHM   = [ 44  38  37];
        mst_filter_info.L1 = [];
        mst_filter_info.L1.wv_ctr = 527;
        mst_filter_info.L1.HWHM   = 7;
        mst_filter_info.L2 = [];
        mst_filter_info.L2.wv_ctr = 445;
        mst_filter_info.L2.HWHM   = 10;
        mst_filter_info.L3 = [];
        mst_filter_info.L3.wv_ctr = 751;
        mst_filter_info.L3.HWHM   = 10;
        mst_filter_info.L4 = [];
        mst_filter_info.L4.wv_ctr = 676;
        mst_filter_info.L4.HWHM   = 10;
        mst_filter_info.L5 = [];
        mst_filter_info.L5.wv_ctr = 867;
        mst_filter_info.L5.HWHM   = 10;
        mst_filter_info.L6 = [];
        mst_filter_info.L6.wv_ctr = 1012;
        mst_filter_info.L6.HWHM   = 21;
    case 'R'
        mst_filter_info.R0 = [];
        mst_filter_info.R0.wv_ctr = [638 551 493];
        mst_filter_info.R0.HWHM   = [ 44  39  38];
        mst_filter_info.R1 = [];
        mst_filter_info.R1.wv_ctr = 527;
        mst_filter_info.R1.HWHM = 7;
        mst_filter_info.R2 = [];
        mst_filter_info.R2.wv_ctr = 447;
        mst_filter_info.R2.HWHM = 10;
        mst_filter_info.R3 = [];
        mst_filter_info.R3.wv_ctr = 805;
        mst_filter_info.R3.HWHM = 10;
        mst_filter_info.R4 = [];
        mst_filter_info.R4.wv_ctr = 908;
        mst_filter_info.R4.HWHM = 11;
        mst_filter_info.R5 = [];
        mst_filter_info.R5.wv_ctr = 937;
        mst_filter_info.R5.HWHM = 11;
        mst_filter_info.R6 = [];
        mst_filter_info.R6.wv_ctr = 1013;
        mst_filter_info.R6.HWHM = 21;
    otherwise
        error('Undefined MASTCAM Eye %s',msteye);
end

end