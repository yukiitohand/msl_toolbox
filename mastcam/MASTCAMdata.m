classdef MASTCAMdata < HSI
    % MASTCAMdata class
    %  MASTCAM data
    
    properties
        lblpath
        lbl
        prop
        PRODUCT_ID
        RMC
        CAM_MDL
        ROVER_NAV
        CAM_MDL_GEO
    end
    
    methods
        function obj = MASTCAMdata(basename,dirpath,varargin)
            obj@HSI(basename,dirpath,varargin{:});
            obj.lblpath = joinPath(dirpath,[basename '.lbl']);
            obj.lbl = pds3lblread(obj.lblpath);
            obj.hdr = mastcam_extract_imghdr_from_lbl(obj.lbl);
            obj.prop = getProp_basenameMASTCAM(basename);
            obj.CAM_MDL = get_cammera_model(obj);
            obj.RMC = get_rmc(obj);
            obj.PRODUCT_ID = obj.lbl.PRODUCT_ID;
            obj.ROVER_NAV  = obj.get_rover_nav();

        end
        function [cahvor_mdl] = get_cammera_model(obj)
            [cahvor_mdl] = mastcam_get_cahvor_model(obj.lbl);
        end
        function [rmc_mdl] = get_rmc(obj)
            [rmc_mdl] = mastcam_get_RMC(obj.lbl);
        end
        function [rover_nav] = get_rover_nav(obj)
            [rover_nav] = get_ROVER_NAV_from_telemetry(...
                obj.RMC.SITE,obj.RMC.DRIVE,obj.RMC.POSE);
        end
        
        function [] = get_CAM_MDL_GEO(obj)
            [imxy_direc_rov] = get_3d_pointing_from_CAHV(...
                [obj.hdr.lines,obj.hdr.samples],obj.CAM_MDL,'gpu',0);
            cmmdl_A_rov0 = obj.ROVER_NAV.rot_mat * obj.CAM_MDL.A';
            cmmdl_C_rov0 = obj.ROVER_NAV.rot_mat * obj.CAM_MDL.C';
            imxy_direc_rov_2d = reshape(imxy_direc_rov,[obj.hdr.lines*obj.hdr.samples,3])';
            imxy_direc_rov0_2d = obj.ROVER_NAV.rot_mat * imxy_direc_rov_2d;
            imxy_direc_rov0 = reshape(imxy_direc_rov0_2d',[obj.hdr.lines,obj.hdr.samples,3]);
            cmmdl_C_geo = cmmdl_C_rov0 + [obj.ROVER_NAV.NORTHING; 
                              obj.ROVER_NAV.EASTING;
                              -obj.ROVER_NAV.ELEVATION];
                          
            obj.CAM_MDL_GEO.A_rov0 = cmmdl_A_rov0;
            obj.CAM_MDL_GEO.C_rov0 = cmmdl_C_rov0;
            obj.CAM_MDL_GEO.C_geo  = cmmdl_C_geo;
            obj.CAM_MDL_GEO.imxy_direc_rov = imxy_direc_rov;
            obj.CAM_MDL_GEO.imxy_direc_rov0 = imxy_direc_rov0;
            
        end
        
            
    end
end

