MSLOrthodata = MSLGaleOrthoMosaic_v3('MSL_Gale_Orthophoto_Mosaic_25cm_v3','/Volumes/LaCie/data');

% outpath = joinPath('/Volumes/LaCie/data','MSL_Gale_Orthophoto_Mosaic_25cm_v3_ave10.img');
% mslortho_mosaic_downsample_mex(MSLOrthodata.imgpath,MSLOrthodata.hdr,outpath,10);

outpath = joinPath('/Volumes/LaCie/data','MSL_Gale_Orthophoto_Mosaic_25cm_v3_ave80.img');
mslortho_mosaic_downsample_mex(MSLOrthodata.imgpath,MSLOrthodata.hdr,outpath,80);