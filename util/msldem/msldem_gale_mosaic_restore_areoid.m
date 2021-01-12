function msldem_gale_mosaic_restore_areoid(MSLDEMdata,molaAdata,varargin)
% msldem_gale_mosaic_restore_areoid(MSLDEMdata,molaAdata,varargin)
%  Restore areoid to get apparent topography of the MSLDEM
% INPUTS
%   MSLDEMdata: MSLGaleDEMMosaic_v3 obj
%   molaAdata : MOLA_MEGTRdata obj, areoid
% OPTIONAL Parameters
%   "COMPENSATE_OFFSET": whether or not to compensate the offset present in
%   the molaAdata (33960km)

% 

S = MSLDEMdata.hdr.samples;
L   = MSLDEMdata.hdr.lines;

% create a big matrix
out_dirpath  = MSLDEMdata.dirpath;
out_basename = [ MSLDEMdata.basename '_ra'];
out_imgpath  = joinPath(out_dirpath,[out_basename '.img']);

msldem_restore_areoid_create_image_mex( out_imgpath,MSLDEMdata.hdr);


% divide into blocks.
sz_blck_s = 8000;
sz_blck_l = 6000;
iter_line = ceil(L/sz_blck_l);
iter_sample = ceil(S/sz_blck_s);
for l=1:iter_line
    line_strt = 1 + sz_blck_l * (l-1);
    line_end = min(sz_blck_l*l ,L);
    for s=1:iter_sample
        smpl_strt = 1 + sz_blck_s * (s-1);
        smpl_end = min(sz_blck_s*s, S);
        
        xrange = [smpl_strt smpl_end];
        yrange = [line_strt line_end];
        
        [imgA_ua,img_dem,xmsldem,ymsldem] = mola_megdr_upsample_areoid(...
            molaAdata,MSLDEMdata,xrange,yrange,'Pixels');
        
        % compensate the areoid only on 
        img_new = img_dem + imgA_ua;
        
        % save the block to the file
        sample_offset = smpl_strt-1; 
        line_offset = line_strt-1;
        samplesc = smpl_end - smpl_strt + 1;
        linesc   = line_end - line_strt + 1;
        msldem_restore_areoid_fill_imageblock_mex(out_imgpath,...
            MSLDEMdata.hdr,sample_offset,line_offset,samplesc,linesc,...
            img_new);
        
        
    end
end

% 
