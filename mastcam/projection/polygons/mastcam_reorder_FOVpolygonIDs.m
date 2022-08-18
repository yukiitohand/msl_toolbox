function [impolid_new,Nid_new] = mastcam_reorder_FOVpolygonIDs(impolid,im_sky)

impolid(im_sky) = -1;

unqid = unique(impolid(:),'sorted'); Nid = length(unqid);

Pid = zeros(1,Nid); idvalid_idx_strt = 1;
for i=1:Nid
    if unqid(i)<0
        idvalid_idx_strt = idvalid_idx_strt+1;
    end
    Pid(i) = sum(impolid(:)==unqid(i));
end
unqid = unqid(idvalid_idx_strt:end);
Pid = Pid(idvalid_idx_strt:end);
Nid_new = length(unqid);

[Pid_srtd,iPid_srtd] = sort(Pid,'descend');
unqid_srtd = unqid(iPid_srtd);

impolid_new = zeros(size(impolid));
impolid_new(im_sky) = -1;
for i=1:Nid_new
    impolid_new(impolid==unqid_srtd(i)) = i;
end

end