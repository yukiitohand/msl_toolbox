function [north,east,elev] = mastcam_get_FOVenvelope(mastcam_NEE)
%[north,east,elev] = mastcam_get_FOVenvelope(mastcam_NEE)
%   GET an approximate enveloped of UFOV mask as a surrounding polygon from
%   mastcam_NEE without actually perfoming the costly caluculation of UFOV 
%   mask.
%  INPUTS
%   mastcam_NEE: [L_im x S_im x 3]
%    North-East-Elevation
%  OUTPUTS
%   north : north-coordinates
%   east  : east-coordinates
%   elev  : elevation-values

[L_im,S_im,~] = size(mastcam_NEE);

bw_mastcam_valid = ~isinf(mastcam_NEE);

imenvl_xy = bwboundaries(bw_mastcam_valid);

if length(imenvl_xy) == 1
    imenvl_xy = imenvl_xy{1};
else
    error('Multiple region is selected');
end

imenvl_xy1d = imenvl_xy(:,1) + L_im*(imenvl_xy(:,2)-1);

mastcam_NEE_1d = reshape(mastcam_NEE,[L_im*S_im,3]);

nee_poly = mastcam_NEE_1d(imenvl_xy1d,:);
north = nee_poly(:,1); east = nee_poly(:,2); elev = nee_poly(:,3);

end

