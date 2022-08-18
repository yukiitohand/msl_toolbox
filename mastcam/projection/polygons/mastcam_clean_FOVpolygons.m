function [impolid_new2] = mastcam_clean_FOVpolygons(impolid,mastcam_range)

im_sky = isinf(mastcam_range);

% reorder impolid
[impolid_new,Nid_new] = mastcam_reorder_FOVpolygonIDs(impolid,im_sky);

% fill the polygons
polid_rm_list = [];
for i=1:Nid_new
    impoli = (impolid_new==i);
    impoli_fill = imfill(impoli,'holes');
    holesi = (impoli ~= impoli_fill);
    holesid = impolid_new(holesi);
    holesid = unique(holesid(:));
    polid_rm_list = union(polid_rm_list,holesid);
end

polid_keep_list = setdiff(1:Nid_new,polid_rm_list);

impolid_new2 = zeros(size(impolid));
impolid_new2(im_sky) = -1;
for i=1:length(polid_keep_list)
    impoli = (impolid_new==polid_keep_list(i));
    impoli_fill = imfill(impoli,'holes');
    impolid_new2(impoli_fill) = i;
end

% reorder impolid
[impolid_new2,Nid_new2] = mastcam_reorder_FOVpolygonIDs(impolid_new2,im_sky);



end

