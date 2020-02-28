function [is_hidden,hidden_idx] = find_hidden_location_all(xyz_roted,lpvo,xy_im_mask)

[~,L,S] = size(xyz_roted);
is_hidden = false(L,S);
hidden_idx = nan(L,S,2);
parfor lt = 1:L
    [is_hidden(lt,:),hidden_idx(lt,:,:)] = wrapper_find_hidden_location(...
        lt,xy_im_mask,xyz_roted,lpvo);
end

end

function [is_hidden,hidden_idx] = wrapper_find_hidden_location(...
    lt,xy_im_mask,xyz_roted,lpvo)

[~,L,S] = size(xyz_roted);
is_hidden = false(1,S);
hidden_idx = nan(1,S,2);
for st = 1:S
    if xy_im_mask(lt,st)
        [is_hidden(st),hidden_idx(:,st,:)] = find_hidden_location(...
            xyz_roted,lpvo,lt,st,xy_im_mask);
    end
end

end