classdef MASTCAMCameraProjectionMSLDEM < handle
    %MASTCAMCameraProjectionMSLDEM
    %  class handling camera projection.
    %  Properties
    %    MASTCAMdata: MASTCAMdata or MASTCAMgroup_eye class object
    %    MSLDEMdata: HSI class obj, 
    %      MSLDEMdata at
    %      https://astrogeology.usgs.gov/search/map/Mars/MarsScienceLaboratory/
    %      Mosaics/MSL_Gale_DEM_Mosaic_10m
    %    msldemc_imFOVmask: boolean image, [L_demc x S_demc x 1]
    %      true if in the FOV, false otherwise.
    %      FOV should be true if the linear transformation shows the value
    %      between -200 and 200+image_edge. The values outside of this are
    %      considered to be inaccurate or invalid. In addition, pixels within
    %      the distance 50m are considered to be in FOV, just in case.
    %    msldemc_imxy: [L_demc x L_demc x 2] 
    %      page 1 is the x coordinate in the image frame,
    %      page 2 is the y coordinate in the image frame.
    %      x and y coordinate values outside of the range between -200 and 
    %      200+image_edge are replaced with nans because the values outside of
    %      this is far away from the calibration range, therefore accuracy is
    %      not guaranteed. Size depend on FOV. The rectangle that minimally
    %      encloses the FOV.
    %    mastcam_NEE: [L_im x S_im x 3]
    %      pages 1,2,3 are northing, easting, and elevation. 
    %    mastcam_msldemc_ref: [L_im x S_im x 3]
    %      indicate which triangle in the MSLDEMdata, the pixel is located.
    %      pages 1,2,3 are sample, line, and j. j is either 1 or 2, indicating
    %      which triangle at the (sample, line).
    %     mastcam_range: [L_im x S_im]
    %      range for each pixel.
    %     mastcam_msldemc_nn: [L_im x S_im x 2]
    %      nearest neighbor pixels in DDR image. The first page is the sample
    %      indices and the second is the line indexes.
    %     msldemc_imFOVmask_nh: [L_demc x S_demc x 1]
    %       Boolean, imFOV_mask with hidden points are removed.
    
    properties
        MASTCAMdata
        MSLDEMdata
        msldemc_hdr_imxy
        msldemc_imxy
        msldemc_imFOVmask
        msldemc_imFOVmask_nh
        mastcam_NEE
        mastcam_range
        mastcam_msldemc_ref
        mastcam_msldemc_nn
        msldemc_imFOVmash_nh2
    end
    
    methods
        function obj = MASTCAMCameraProjectionMSLDEM(mstdata_obj,MSLDEMdata_obj)
            obj.MASTCAMdata = mstdata_obj;
            obj.MSLDEMdata  = MSLDEMdata_obj;
        end
        
        function proj_MSLDEM2mastcam_old(obj)
            [MSLDEMprj] = proj_MSLDEM2mastcam_v2(obj.MSLDEMdata,obj.MASTCAMdata);
            l1 = MSLDEMprj.hdr_imxy.line_offset+1;
            lend = MSLDEMprj.hdr_imxy.line_offset+MSLDEMprj.hdr_imxy.lines;
            s1 = MSLDEMprj.hdr_imxy.sample_offset+1;
            send = MSLDEMprj.hdr_imxy.sample_offset+MSLDEMprj.hdr_imxy.samples;
            dem_imFOV_mask_crop = MSLDEMprj.imFOV_mask(l1:lend,s1:send);
            obj.msldemc_imxy = MSLDEMprj.imxy;
            obj.msldemc_imFOVmask = dem_imFOV_mask_crop;
            obj.msldemc_hdr_imxy = MSLDEMprj.hdr_imxy;
        end
        
        function proj_MSLDEM2mastcam(obj)
            [obj.msldemc_imFOVmask,obj.msldemc_imxy,obj.msldemc_hdr_imxy] = proj_MSLDEM2mastcam_v3(obj.MSLDEMdata,obj.MASTCAMdata);
            
        end
        
        function proj_mastcam2MSLDEM(obj)
            proj_mastcam2MSLDEM_v5_mexw(obj.MASTCAMdata,obj.MSLDEMdata,obj);
        end
        
        function assess_occlusion(obj)
            msldemc_noccFOV = mastcam_assess_occlusion_wMSLDEM_mexw(obj);
            obj.msldemc_imFOVmash_nh2 = msldemc_noccFOV;
        end
        
        
    end
end

