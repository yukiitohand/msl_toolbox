function [specResponse] = mastcam_getSpectralResponse(mstqe_obj,mstbp_obj)
% [specResponse] = mastcam_getSpectralResponse(mstqe_obj,mstbp_obj)
%   Multiplying multiple band pass filters and responsitivities to get
%   total spectral response function.
%  INPUTS
%   mstqe_obj: class object of MSTQuantumEfficiency.
%   mstbp_obj: class object of MSTBandpass. Can be multiple.
%  OUTPUTS
%   specResponse: struct having two fields, "wavelenth" and "response"

wv_step = 1;

N_mstbp = length(mstbp_obj); % the number of bandpass objects.

%==========================================================================
% Right now band pass functions are supposed to be sampled with integer
% wavelength.
%==========================================================================
for i = 1:N_mstbp
    if any(floor(mstbp_obj(i).wavelength)~=mstbp_obj(i).wavelength)
        error('Wavelength values are not integers');
    end
end

%==========================================================================
% Get the minimal intersection of wavelength samples
%==========================================================================
wv_strt = -inf; wv_end = inf;
for i=1:length(mstbp_obj)
    wvi_min = min(mstbp_obj(i).wavelength); 
    wvi_max = max(mstbp_obj(i).wavelength);
    wv_strt = max(wv_strt,wvi_min);
    wv_end  = min(wv_end,wvi_max);
end

wv_strt = max(wv_strt,min(mstqe_obj.wavelength));
wv_end  = min(wv_end, max(mstqe_obj.wavelength));

wv_new = wv_strt:wv_step:wv_end;

%==========================================================================
% Multiplying the transmissions
%==========================================================================
for i=1:length(mstbp_obj)
    idxi = arrayfun(@(x) find(x==mstbp_obj(i).wavelength),wv_new);
    if i==1
        trs_new = mstbp_obj(i).transmission(idxi);
    else
        trs_new = trs_new .* mstbp_obj(i).transmission(idxi);
    end
end

idx_qe = arrayfun(@(x) find(x==mstqe_obj.wavelength),wv_new);
trs_new = trs_new .* mstqe_obj.responsitivity(idx_qe);

%==========================================================================
% Summary
%==========================================================================
specResponse = [];
specResponse.wavelength = wv_new;
specResponse.response = trs_new;

end