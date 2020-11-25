function [crism_imUFOVmask] = convert_UFOVmask_msldem2crism_wrapper(...
    msldemc_imUFOVmask,msldemc_imUFOVmask_hdr,GLTdata)


msldemc_northing = msldemc_imUFOVmask_hdr.y;
msldemc_easting = msldemc_imUFOVmask_hdr.x;

crism_dn = GLTdata.hdr.map_info.dy;
crism_de = GLTdata.hdr.map_info.dx;
Lcrism   = GLTdata.hdr.lines;
Scrism   = GLTdata.hdr.samples;

crism_n0 = GLTdata.hdr.map_info.mapy - crism_dn * (1-GLTdata.hdr.map_info.image_coords(1));
crism_e0 = GLTdata.hdr.map_info.mapx + crism_de * (1-GLTdata.hdr.map_info.image_coords(2));

crism_imUFOVmask = convert_UFOVmask_msldem2crism_mex(msldemc_imUFOVmask,...
    msldemc_northing,msldemc_easting,crism_n0,crism_e0,crism_dn,crism_de,...
    Lcrism,Scrism);

end