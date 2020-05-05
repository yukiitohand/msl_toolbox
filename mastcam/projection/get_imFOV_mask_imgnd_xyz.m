function [imgndFOV_mask,imgnd_imxy] = get_imFOV_mask_imgnd_xyz(cmmdl,camim_size,imgnd_xyz,varargin)
% [imgndFOV_mask,imgnd_imxy] = get_imFOV_mask_imgnd_xyz(cmmdl,camim_size,imgnd_xyz,varargin)
%  Obtain image FOV for given camera model and image_xyz and give image
%  size.
% INPUTS:
%   cmmdl: obj of CAHVOR_MODEL class
%   camim_size: size of the camera image [L_im, S_im]
%   imgnd_xyz: (L_g x S_g x 3)
%     coordinate in XYZ ground reference coordinate for image pixels.
%     The first, second, and third page correspond to coordinates of X, Y,
%     and Z.
%  OUTPUTS:
%   imgndFOV_mask: boolean (L_g x S_g), 
%    true if a pixel is potentially in FOV. 
%   imgnd_imxy: (L_g x S_g x 2)
%    representing the (x,y) coordinate in the image for each pixel in the
%    imgnd_xyz. The first page is x coordinate and the second is y.
%  OPTIONAL PARAMETERS
%   'RESOLUTION_FACTOR'
%    factor for amplifying the estimated resolution.
%    (default) 1.05
%
% TIPS
%  Make sure cmmdl and imgnd_xyz are defined in the same coordinte system.

resol_amp_factor = 1.05;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'RESOLUTION_FACTOR'
                resol_amp_factor = varargin{i+1};
            otherwise
                error('Unrecognized option: %s',varargin{i});
        end
    end
end

L_g = size(imgnd_xyz,1); S_g = size(imgnd_xyz,2);
L_im = camim_size(1); S_im = camim_size(2);
cmmdl.get_image_plane_unit_vectors();

imgnd_xyz_vec2d_mC = reshape(imgnd_xyz,[L_g*S_g,3])' - cmmdl.C';
dst2plane_inprod_vec = cmmdl.A*imgnd_xyz_vec2d_mC;

imgnd_imxy_vec2d = ([cmmdl.H;cmmdl.V] * imgnd_xyz_vec2d_mC)./dst2plane_inprod_vec;

% get resolution
[imgnd_xyz_resol] = get_instantaneous_resolution(imgnd_xyz);
resol_HVdash_vec2d = [abs(cmmdl.Hdash);abs(cmmdl.Vdash)] ...
    * reshape(imgnd_xyz_resol,[L_g*S_g,3])'.*resol_amp_factor;
% resol_HVdash = reshape(resol_HVdash_vec2d',[L_g,S_g,2]);

% get margin
hsvs_d_vec2d = [cmmdl.hs;cmmdl.vs].*resol_HVdash_vec2d./abs(dst2plane_inprod_vec);
% hsvs_d = reshape(hsvs_d_vec2d',[L_g,S_g,2]);

% test
imgnd_FOV_mask_xy_vec = all(...
    and( (-0.5-hsvs_d_vec2d) < imgnd_imxy_vec2d,...
         imgnd_imxy_vec2d < ([S_im;L_im]-0.5+hsvs_d_vec2d) )...
    ,1);

right_dir_vec = (dst2plane_inprod_vec > 0);

imgnd_FOV_mask_xyd = and(right_dir_vec,imgnd_FOV_mask_xy_vec);
imgnd_imxy_vec2d(:,imgnd_FOV_mask_xyd==0) = nan;
imgnd_FOV_mask_xyd = reshape(imgnd_FOV_mask_xyd,[L_g,S_g]);

% safeguarding for close distance (some pixels located in the
% opposite direction also needs 
[safeguard_mask] = get_safeguard_imFOV_mask_closedistance(...
    imgnd_FOV_mask_xyd,reshape(right_dir_vec,[L_g,S_g]));

imgndFOV_mask = or(imgnd_FOV_mask_xyd,safeguard_mask);

imgnd_imxy = reshape(imgnd_imxy_vec2d',[L_g,S_g,2]);

end