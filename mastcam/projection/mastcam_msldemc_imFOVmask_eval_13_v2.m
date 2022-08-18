function [msldemc_imFOVmask_new] = mastcam_msldemc_imFOVmask_eval_13_v2( ...
    msldemc_imFOVmask_obj,mastcam_ref_msldem_obj)
% [msldemc_imFOVmask_new] = mastcam_msldemc_imFOVmask_eval_13_v2( ...
%     msldemc_imFOVmask_obj,mastcam_ref_msldem_obj)
% msldemc(,) = 3 is only expected at the close range pixels, so if the
% pixels does not observed in the mastcam_msldemc_ref, then the pixel won't
% be in the FOV.
%  INPUTS
%   msldemc_imFOVmask_obj  : 
%   mastcam_ref_msldem_obj : 

% Validate the input variables.
validateattributes(msldemc_imFOVmask_obj, ...
    {'ENVIRasterSingleLayerMSLDEMCProj','ENVIRasterMultBandMSLDEMCProj'}, ...
    {},mfilename,'msldemc_imFOVmask_obj');

validateattributes(mastcam_ref_msldem_obj,{'ENVIRasterMultBand'}, ...
    {},mfilename,'mastcam_ref_msldem_obj');

% Read imFOVmask
if isempty(msldemc_imFOVmask_obj.img)
    msldemc_imFOVmask_img = msldemc_imFOVmask_obj.readimg('precision','raw');
else
    msldemc_imFOVmask_img = msldemc_imFOVmask_obj.img;
end

% Read ref_msldem and modify indices so that they align with msldemc_imFOVmask
if isempty(mastcam_ref_msldem_obj.img)
    mastcam_ref_msldem_img = mastcam_ref_msldem_obj.readimg('precision','raw');
else
    mastcam_ref_msldem_img = mastcam_ref_msldem_obj.img;
end
ref_valid = mastcam_ref_msldem_img(:,:,1)>-1;
ref_msldemc_img1 = mastcam_ref_msldem_img(:,:,1);
ref_msldemc_img1(ref_valid) = ref_msldemc_img1(ref_valid) - msldemc_imFOVmask_obj.chdr.sample_offset;
ref_msldemc_img2 = mastcam_ref_msldem_img(:,:,2);
ref_msldemc_img2(ref_valid) = ref_msldemc_img2(ref_valid) - msldemc_imFOVmask_obj.chdr.line_offset;
mastcam_ref_msldemc_img = cat(3,ref_msldemc_img1,ref_msldemc_img2, ...
    mastcam_ref_msldem_img(:,:,3));

% Main computation
[row3,col3] = find(msldemc_imFOVmask_img==3);
[row1,col1] = find(msldemc_imFOVmask_img==1);

msldemc_imFOVmask_new = msldemc_imFOVmask_img;

if ~isempty(row3)
    
    %
    mastcam_refc = mastcam_ref_msldemc_img(:,:,1);
    mastcam_refl = mastcam_ref_msldemc_img(:,:,2);
    
    N = length(row3);
    for n=1:N
        c = col3(n); l = row3(n);
        c_min = max(1,c-1);
        l_min = max(1,l-1);
        cl_is_inFOV = false;
        for ll=l_min:l
            for cc=c_min:c
                if any(and(mastcam_refc==cc,mastcam_refl==ll),'all')
                    cl_is_inFOV = true;
                end
            end
        end
        if cl_is_inFOV
            msldemc_imFOVmask_new(l,c) = 4;
        else
            msldemc_imFOVmask_new(l,c) = 2;
        end
        
    end
    
end

if ~isempty(row1)
    
    % mastcam_ref's values are indexes starting at 1 (MATLAB style)
    mastcam_refc = mastcam_ref_msldemc_img(:,:,1);
    mastcam_refl = mastcam_ref_msldemc_img(:,:,2);
    
    N = length(row1);
    for n=1:N
        c = col1(n); l = row1(n);
        c_min = max(1,c-1);
        l_min = max(1,l-1);
        cl_is_inFOV = false;
        for ll=l_min:l
            for cc=c_min:c
                if any(and(mastcam_refc==cc,mastcam_refl==ll),'all')
                    cl_is_inFOV = true;
                end
            end
        end
        if cl_is_inFOV
            msldemc_imFOVmask_new(l,c) = 1;
        else
            msldemc_imFOVmask_new(l,c) = 0;
        end
        
    end
    
end

end