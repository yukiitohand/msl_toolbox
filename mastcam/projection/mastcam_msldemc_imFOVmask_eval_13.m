function [msldemc_imFOVmask_new] = mastcam_msldemc_imFOVmask_eval_13(msldemc_imFOVmask,mastcam_msldemc_ref)
% [msldemc_imFOVmask_new] = mastcam_msldemc_imFOVmask_eval_3(msldemc_imFOVmask,mastcam_msldemc_ref)
% msldemc(,) = 3 is only expected at the close range pixels, so if the
% pixels does not observed in the mastcam_msldemc_ref, then the pixel won't
% be in the FOV.

if isempty(msldemc_imFOVmask)
    error('non empty msldemc_imFOVmask is required.');
end

if isempty(mastcam_msldemc_ref)
    error('non empty mastcam_msldemc_ref is required.');
end

[row3,col3] = find(msldemc_imFOVmask==3);
[row1,col1] = find(msldemc_imFOVmask==1);

msldemc_imFOVmask_new = msldemc_imFOVmask;

if ~isempty(row3)
    
    % mastcam_ref's values are indexes starting at 0 for the first row and column.
    % bring back to MATLAB indexes starting at 1.
    mastcam_refc = mastcam_msldemc_ref(:,:,1) + 1;
    mastcam_refl = mastcam_msldemc_ref(:,:,2) + 1;
    
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
    
    % mastcam_ref's values are indexes starting at 0 for the first row and column.
    % bring back to MATLAB indexes starting at 1.
    mastcam_refc = mastcam_msldemc_ref(:,:,1) + 1;
    mastcam_refl = mastcam_msldemc_ref(:,:,2) + 1;
    
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