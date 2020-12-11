function [xyz_resol] = get_instantaneous_resolution(xyz)
% Evaluate approximate resolution in the first two dimensions at each pixel
%  INPUTS
%    xyz: (L x S x B) image, storing geographical information.
%      First two dimensions corresopond to spatial dimensions.
%  OUTPUTS
%   xyz_resol: (L x S x B) image, storing resolution at each pixel for each
%   geographical component.
xyz_pd = padarray(xyz,[1,1,0],nan,'both');
% get resolution
xyz_resol = nan(size(xyz));
for i=1:3
    for j=1:3
        xyz_resol = max(...
                xyz_resol, ...
                abs( xyz - xyz_pd(i:end-3+i,j:end-3+j,:) )...
            );
    end
end

end