function [mask_mastcam_cellar,mask_msldemc_UFOVmask_cellar,...
    msldemc_mask_all,msldemc_mask_all_hdr] = mstproj_combine_mask(MSTproj,x_east,y_north)

% get the size of the combined mask
s1 = inf; send = 1;
l1 = inf; lend = 1;
for i=1:length(MSTproj)
    s1 = min(MSTproj(i).msldemc_imUFOVhdr.sample_offset+1,s1);
    l1 = min(MSTproj(i).msldemc_imUFOVhdr.line_offset+1,l1);
    send = max(MSTproj(i).msldemc_imUFOVhdr.sample_offset+MSTproj(i).msldemc_imUFOVhdr.samples,send);
    lend = max(MSTproj(i).msldemc_imUFOVhdr.line_offset+MSTproj(i).msldemc_imUFOVhdr.lines,lend);
end

smplofst_all = s1-1; lnofst_all = l1-1;
smpls_all = send-s1+1; lns_all = lend-l1+1;

% store header information of the combined mask
msldemc_mask_all_hdr = [];
msldemc_mask_all_hdr.lines = lns_all;
msldemc_mask_all_hdr.samples = smpls_all;
msldemc_mask_all_hdr.sample_offset = smplofst_all;
msldemc_mask_all_hdr.line_offset = lnofst_all;


msldemc_mask_all = false(lns_all,smpls_all);

mask_mastcam_cellar = cell(1,length(MSTproj));
mask_msldemc_UFOVmask_cellar = cell(1,length(MSTproj));
for i=1:length(MSTproj)
    % =============================================================
    % First get the selected point/area in the MASTCAM image at 
    % MSLDEM resolution.
    [x_mst,y_mst,x_dem,y_dem,mask_mastcam,mask_msldemc_UFOVmask]...
        = MSTproj(i).get_SelectMask_fromNorthEast(x_east,y_north);

    mask_mastcam_cellar{i} = mask_mastcam;
    mask_msldemc_UFOVmask_cellar{i} = mask_msldemc_UFOVmask;
    
    si1 = MSTproj(i).msldemc_imUFOVhdr.sample_offset-msldemc_mask_all_hdr.sample_offset+1;
    li1 = MSTproj(i).msldemc_imUFOVhdr.line_offset-msldemc_mask_all_hdr.line_offset+1;
    siend = si1+MSTproj(i).msldemc_imUFOVhdr.samples-1;
    liend = li1+MSTproj(i).msldemc_imUFOVhdr.lines-1;
    
    if ~isempty(mask_msldemc_UFOVmask)
        msldemc_mask_all(li1:liend,si1:siend) = or(msldemc_mask_all(li1:liend,si1:siend),mask_msldemc_UFOVmask);
    end

end

% mask_mastcam_all
% for i=1:length(MSTproj)
%     si1 = MSTproj(i).msldemc_imUFOVhdr.sample_offset-msldemc_mask_all_hdr.sample_offset+1;
%     li1 = MSTproj(i).msldemc_imUFOVhdr.line_offset-msldemc_mask_all_hdr.line_offset+1;
%     siend = si1+mstproj.msldemc_imUFOVhdr.samples-1;
%     liend = li1+mstproj.msldemc_imUFOVhdr.lines-1;
%     msldemc_mask_all_i = msldemc_mask_all(li1:liend,si1:siend);
% end

end