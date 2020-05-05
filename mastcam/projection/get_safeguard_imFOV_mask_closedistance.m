function [safeguard_mask] = get_safeguard_imFOV_mask_closedistance(imFOV_mask,right_dir)
% [safeguard_mask] = get_safeguard_imFOV_mask_closedistance(imFOV_mask,right_dir)
%  Get a mask of pixels in the image that are oppsoite in the direction of 
%  line of sight but may contribute to the camera FOV due to closeness to 
%  the camera center.
%  INPUTS:
%   imFOV_mask: Boolean, two dimensional mask
%   right_dir: Boolean, same size as imFOV_mask. true if the pixel is in
%              the same direction as the line of sight.
%  OUTPUTS:
%   safeguard_mask: Boolean, same size as imFOV_mask

% safeguarding for close distance (some pixels located in the
% opposite direction also needs 
% ddr_imFOV_mask_xyd = reshape(imFOV_mask,[L_ddr,S_ddr]);
imFOV_mask_xyd_pd = padarray(imFOV_mask,[1 1],'replicate','both');

imFOV_mask_nbr = false(size(imFOV_mask));
for i=1:3
    for j=1:3
        imFOV_mask_nbr = or(imFOV_mask_nbr,imFOV_mask_xyd_pd(i:end-3+i,j:end-3+j));
    end
end
safeguard_mask = and(~right_dir,imFOV_mask_nbr);

end