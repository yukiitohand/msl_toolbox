function [polygoncell] = mastcam_get_FOVpolygons_Vectors(impolid)

impolid_unique = unique(impolid(:));
id_pos_idx = (impolid_unique > 0);
impolid_unique = impolid_unique(id_pos_idx);
Nid = length(impolid_unique);

polygoncell = cell(1,Nid);
for i=1:Nid
    impoli = (impolid==impolid_unique(i));
    impoli_fill = bwboundaries(impoli,4,'noholes');
    polygoncell{i} = impoli_fill{1};
end

end
