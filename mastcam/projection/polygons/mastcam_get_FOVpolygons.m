function [impolid] = mastcam_get_FOVpolygons(mastcam_NEE,ifovV,ifovH)
%[north,east,elev] = mastcam_get_FOVenvelope(mastcam_NEE)
%   GET an approximate polygons of UFOV mask without actually perfoming the
%   costly caluculation of UFOV mask.
%  INPUTS
%   mastcam_NEE: [L_im x S_im x 3]
%    North-East-Elevation
%  OUTPUTS
%   polygons_nee: cell array, each element has [ * x 3 ] components of the 


[L_im,S_im,~] = size(mastcam_NEE);

% PrepareIMage for POLygon group IDs.
impolid = nan(L_im,S_im);
impolid(1,1) = 1;
maxpolid = 1;

% first column
for l=2:L_im
    % if l==237
    %     keyboard;
    % end
    if sqrt(sum((mastcam_NEE(l,1,:)-mastcam_NEE(l-1,1,:)).^2,'all')) <= max(ifovV(l,1),ifovV(l-1,1))
        impolid(l,1) = impolid(l-1,1);
    elseif isnan(mastcam_NEE(l,1,1)) && isnan(mastcam_NEE(l-1,1,1))
        impolid(l,1) = impolid(l-1,1);
    else
        maxpolid = maxpolid+1;
        impolid(l,1) = maxpolid;
    end
end

% first line
for s=2:S_im    
    if sqrt(sum((mastcam_NEE(1,s,:)-mastcam_NEE(1,s-1,:)).^2,'all')) <= max(ifovH(1,s),ifovH(1,s-1))
        impolid(1,s) = impolid(1,s-1);
    elseif isnan(mastcam_NEE(1,s,1)) && isnan(mastcam_NEE(1,s-1,1))
        impolid(1,s) = impolid(1,s-1);
    else
        maxpolid = maxpolid+1;
        impolid(1,s) = maxpolid;
    end
end

% a = impolid(220,1);

figure;
for s=2:S_im
    for l=2:L_im
        if s==78 && l==195
            % keyboard;
        end
        if sqrt(sum((mastcam_NEE(l,s,:)-mastcam_NEE(l-1,s,:)).^2,'all')) <= max(ifovV(l,s),ifovV(l-1,s))
            impolid(l,s) = impolid(l-1,s);
            if sqrt(sum((mastcam_NEE(l,s,:)-mastcam_NEE(l,s-1,:)).^2,'all')) <= max(ifovH(l,s),ifovH(l,s-1)) ...
                    && (impolid(l,s-1) ~= impolid(l,s))
                grp_id_rm = impolid(l,s);
                grp_id_keep = impolid(l,s-1);
                % impolid(impolid==grp_id_rm) = grp_id_keep;
                % ll=l;
                % while ll>0 && (impolid(ll,s) == grp_id_rm)
                %     impolid(ll,s) = grp_id_keep;
                %     ll = ll-1;
                % end
                ss = s;
                while ss>0 && (any(impolid(:,ss)==grp_id_rm))
                    impolid(impolid(:,ss)==grp_id_rm,ss) = grp_id_keep;
                    ss = ss-1;
                end
            end
        elseif sqrt(sum((mastcam_NEE(l,s,:)-mastcam_NEE(l,s-1,:)).^2,'all')) < max(ifovH(l,s),ifovH(l,s-1))
            impolid(l,s) = impolid(l,s-1);
        elseif isnan(mastcam_NEE(l,s,1)) && isnan(mastcam_NEE(l,s-1,1))
            impolid(l,s) = impolid(l,s-1);
        elseif isnan(mastcam_NEE(l,s,1)) && isnan(mastcam_NEE(l-1,s,1))
            impolid(l,s) = impolid(l-1,s);
        else
            maxpolid = maxpolid+1;
            impolid(l,s) = maxpolid;
        end
    end
    % imagesc(impolid); set(gca,'DataAspectRatio',[1,1,1]); drawnow;
    % pause(0.1);
end



end