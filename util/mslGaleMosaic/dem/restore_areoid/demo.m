pds_geosciences_node_setup
basenameA = 'MEGA90N000EB';
molaAdata = MOLA_MEGTRdata(basenameA,'');
MSLDEMdata = MSLGaleDEMMosaic_v3('MSL_Gale_DEM_Mosaic_1m_v3_dave',...
    '/Users/yukiitoh/data/');


msldem_gale_mosaic_restore_areoid(MSLDEMdata,molaAdata);