function [is_hidden,hidden_ls] = find_hidden_location(xyz_roted,lpvo,lt,st,xy_im_mask)
% [is_hidden,hidden_ls] = find_hidden_location(xyz_roted,lpvo,lpv1)
%   evaluate the the specified location lpv1 is hidden or not from lpvo.
%  INPUTS
%    xyz_roted: [3 x lines x samples]
%    lpvo     : [3 x 1], origin position vector
%    lt,st    : [3 x 1] target position vector location
%    xy_im_mask    : [lines x samples] field of view.
%  OUTPUTS
%    is_hidden: boolean, lpv1 is hidden or not
%    hidden_ls: return idx [l,s], that is determied to block the ray from
%    lpv1.

lpv1 = xyz_roted(:,lt,st); % line position vector
[~,L,S] = size(xyz_roted);
is_hidden = false;
hidden_ls = [nan,nan];
flg = 0;

idx_evaluated = find(xy_im_mask(:));
% l_idx_evaluated = mod(idx_evaluated,L);
% l_idx_evaluated(l_idx_evaluated==0) = L;
s_idx_evaluated = ceil(idx_evaluated/L);
l_idx_evaluated = idx_evaluated - L * (s_idx_evaluated-1);

for j=1:length(idx_evaluated)
    l = l_idx_evaluated(j);
    s = s_idx_evaluated(j);
    
    if l==L, continue; end
    
    % select three points
    if mod(s,2)==1
        ppv1 = xyz_roted(:,l,s); % plane position vector
        ppv2 = xyz_roted(:,l,s+1);
        ppv3 = xyz_roted(:,l+1,s);
        idx_vert = [l,s;l,s+1;l+1,s];
    else
        ppv1 = xyz_roted(:,l,s);
        ppv2 = xyz_roted(:,l+1,s);
        ppv3 = xyz_roted(:,l+1,s-1);
        idx_vert = [l,s;l+1,s;l+1,s-1];
    end
    if any(all(idx_vert==[lt,st],2))
        % skip one of the vertices is a target vector.
        continue;
    end
    % test line segment intersect with the plane determined by the
    % three points
    [line_param,is_intersect] = line_plane_intersect(lpvo,lpv1,ppv1,ppv2,ppv3);

    % if the line segment intersect with the plane, then test if the
    % intersection is within the triangle or not.
    if is_intersect
        pipv = lpvo + (lpv1-lpvo)*line_param; % plane intersection position vector
        [plane_param,is_in_face] = get_plane_param_coefficient(ppv1,ppv2,ppv3,pipv);
    else
        is_in_face = false;
    end
    if is_in_face
        is_hidden = true;
        hidden_ls = [l,s];
        flg = true;
    end
    if flg
        break;
    end

end


end