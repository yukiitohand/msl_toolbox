function [imprj] = proj_CRISMimage2MASTCAM(img,mastcam_prj,...
    crism_DDRdata,mastcamdata_obj)
% [imprj] = proj_CRISMimage2MASTCAM(img,mastcam_prj,...
%     crism_DDRdata,mastcamdata_obj)
%  project a CRISM sampled image onto a MASTCAM image space.

switch class(mastcamdata_obj)
    case 'MASTCAMdata'
        L_im = mastcamdata_obj.hdr.lines; S_im = mastcamdata_obj.samples;
    case 'MASTCAMgroup_eye'
        L_im = mastcamdata_obj.L_im; S_im = mastcamdata_obj.S_im;
    otherwise
        error('second input needs to be either MASTCAMdata or MASTCAMgroup_eye class.');
end

L_ddr = crism_DDRdata.hdr.lines; S_ddr = crism_DDRdata.hdr.samples;
nB = size(img,3);

% nearest neighbor based projection.
idx_nn = (mastcam_prj.DDR_nn(:,:,1)-1)*L_ddr + mastcam_prj.DDR_nn(:,:,2);
idx_nn = idx_nn(:);
idx_nn_1d_nisnan = ~isnan(idx_nn);
idx_nn_1d_ok = idx_nn(idx_nn_1d_nisnan);

imprj = nan(L_im*S_im,nB);
img = reshape(img,[L_ddr*S_ddr,nB]);

imprj(idx_nn_1d_nisnan,:) = img(idx_nn_1d_ok,:);

imprj = reshape(imprj,[L_im,S_im,nB]);


end