function [NEEpolygon_cell] = mastcam_get_NEEpolygons_from_impolygons(mastcam_NEE,impolygoncell)
% 
[L_im,S_im,~] = size(mastcam_NEE);
Np = length(impolygoncell);
mastcam_NEE_1d = reshape(mastcam_NEE,[L_im*S_im,3]);

NEEpolygon_cell = cell(1,Np);
for i=1:Np
    impolygon_i = impolygoncell{i};
    impolygon_i_1d = impolygon_i(:,1) + (impolygon_i(:,2)-1) * L_im;
    neepolygon_i = mastcam_NEE_1d(impolygon_i_1d,:);
    NEEpolygon_cell{i} = neepolygon_i;
    
end
    

end