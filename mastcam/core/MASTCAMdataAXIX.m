classdef MASTCAMdataAXIX < HSI
    % MASTCAMAIdata class
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
        FILTER_NUMBER
        INSTRUMENT_ID
        L_im
        S_im
        MASTCAMdata_ref % referenced mastcam data
        RADIANCE_FACTOR
        RADIANCE_OFFEST
    end
    
    methods
        function obj = MASTCAMdataAXIX(basename,dirpath,varargin)
            
            global msl_env_vars
            dirpath_AXII = msl_env_vars.dirpath_AXI1;
            
            rover_nav_ver = 'localized_interp';
            rover_nav_mstcode = '';
            mstdata_ref = [];
            varargin_rmidx = [];
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case {'ROVER_NAV_VERSION','ROVER_NAV_VER'}
                            rover_nav_ver = varargin{i+1};
                            varargin_rmidx = [varargin_rmidx i i+1];
                        case 'ROVER_NAV_MSTCAM_CODE'
                            rover_nav_mstcode = varargin{i+1};
                            varargin_rmidx = [varargin_rmidx i i+1];
                        case {'MASTCAMDATA_REF'}
                            mstdata_ref = varargin{i+1};
                            varargin_rmidx = [varargin_rmidx i i+1];
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            
            if isempty(dirpath)
                dirpath = dirpath_AXII;
            end
            
            varargin_HSI = varargin(setdiff(1:length(varargin),varargin_rmidx));
            obj@HSI(basename,dirpath,varargin_HSI{:});
            obj.lblpath = joinPath(dirpath,[basename '.IMG']);
            [obj.lbl,posend] = pds3lblread(obj.lblpath);
                
            obj.hdr = mastcam_extract_imghdr_from_lbl(obj.lbl);
            lbytes = obj.hdr.samples * get_bytes_envi_data_type(obj.hdr.data_type);
            obj.hdr.header_offset = ceil(posend/lbytes)*lbytes;
            obj.L_im = obj.hdr.lines;
            obj.S_im = obj.hdr.samples;
            
            if ~isempty(mstdata_ref)
                % obj.lblpath = joinPath(dirpath,[basename '.lbl']);
                % obj.lbl = pds3lblread(obj.lblpath);
                obj.prop = getProp_basenameMASTCAM(basename);
                [issamescene] = obj.compare_prop_ref(mstdata_ref);
                if ~issamescene
                    error('Input reference cube is not an ancestor of %s',basename);
                end
                % obj.hdr = mstdata_ref.hdr;
                obj.CAM_MDL = mstdata_ref.CAM_MDL;
                obj.RMC = mstdata_ref.RMC;
                % obj.PRODUCT_ID = obj.lbl.PRODUCT_ID;
                obj.ROVER_NAV  = mstdata_ref.ROVER_NAV;
                obj.FILTER_NUMBER = mstdata_ref.FILTER_NUMBER;
                obj.INSTRUMENT_ID = mstdata_ref.INSTRUMENT_ID;
                % obj.get_radiance_factor();
                % obj.L_im = mstdata_ref.L_im;
                % obj.S_im = mstdata_ref.S_im;
                
                
                obj.MASTCAMdata_ref = mstdata_ref;
            end
            
            
        end
        
        function [issamescene] = compare_prop_ref(obj,mst_obj)
            if (obj.prop.sol==mst_obj.prop.sol)...
                    && strcmpi(obj.prop.cam_code,mst_obj.prop.cam_code)...
                    && (obj.prop.seq_id==mst_obj.prop.seq_id)...
                    && (obj.prop.command_num==mst_obj.prop.command_num)...
                    && (obj.prop.cdpid_counter==mst_obj.prop.cdpid_counter)...
                    && (obj.prop.unique_cdpid==mst_obj.prop.unique_cdpid)...
                    && strcmpi(obj.prop.product_type,mst_obj.prop.product_type)...
                    && strcmpi(obj.prop.gop,mst_obj.prop.gop)...
                    && (obj.prop.version==mst_obj.prop.version)
                issamescene = true;
            else
                issamescene = false;
            end
                
                
        end
        
        function [cahvor_mdl] = get_cammera_model(obj)
            [cahvor_mdl] = obj.MASTCAMdata_ref.get_cammera_model();
        end
        
        function [rmc_mdl] = get_rmc(obj)
            [rmc_mdl] = obj.MASTCAMdata_ref.get_rmc();
        end
        
        function [rover_nav] = get_rover_nav(obj,varargin)
            [rover_nav] = obj.MASTCAMdata_ref.get_rover_nav();
            % [rover_nav] = get_ROVER_NAV_from_mslplc_view_csv(...
            %     obj.RMC.SITE,obj.RMC.DRIVE,obj.RMC.POSE,varargin{:});
        end
        
%         function [] = get_CAM_MDL_GEO(obj)
%             [imxy_direc_rov] = get_3d_pointing_from_CAHV_v2(...
%                 [obj.hdr.lines,obj.hdr.samples],obj.CAM_MDL);
%             cmmdl_A_rov0 = obj.ROVER_NAV.rot_mat * obj.CAM_MDL.A';
%             cmmdl_C_rov0 = obj.ROVER_NAV.rot_mat * obj.CAM_MDL.C';
%             imxy_direc_rov_2d = reshape(imxy_direc_rov,[obj.hdr.lines*obj.hdr.samples,3])';
%             imxy_direc_rov0_2d = obj.ROVER_NAV.rot_mat * imxy_direc_rov_2d;
%             imxy_direc_rov0 = reshape(imxy_direc_rov0_2d',[obj.hdr.lines,obj.hdr.samples,3]);
%             cmmdl_C_geo = cmmdl_C_rov0 + [obj.ROVER_NAV.NORTHING; 
%                               obj.ROVER_NAV.EASTING;
%                               -obj.ROVER_NAV.ELEVATION];
%                           
%             obj.CAM_MDL_GEO.A_rov0 = cmmdl_A_rov0;
%             obj.CAM_MDL_GEO.C_rov0 = cmmdl_C_rov0;
%             obj.CAM_MDL_GEO.C_geo  = cmmdl_C_geo;
%             obj.CAM_MDL_GEO.imxy_direc_rov = imxy_direc_rov;
%             obj.CAM_MDL_GEO.imxy_direc_rov0 = imxy_direc_rov0;
%             
%         end
        
        function get_instrument_id(obj)
            obj.INSTRUMENT_ID = obj.MASTCAMdata_ref.get_instrument_id();
        end
        
        function get_filter_number(obj)
            obj.FILTER_NUMBER = obj.MASTCAMdata_ref.get_filter_number();
        end
        
        
        function get_radiance_factor(obj)
            obj.RADIANCE_FACTOR = obj.lbl.OBJECT_IMAGE.SCALING_FACTOR;
            obj.RADIANCE_OFFSET = obj.lbl.OBJECT_IMAGE.OFFSET;
        end
        
        function [img_iof] = get_IoF(obj)
            if isempty(obj.img)
                obj.readimg();
            end
            switch obj.FILTER_NUMBER
                case 0
                    img_iof = reshape(obj.RADIANCE_FACTOR,[1,1,3]) .* obj.img + reshape(obj.RADIANCE_OFFEST,[1,1,3]);
                otherwise
                    img_iof = obj.RADIANCE_FACTOR .* obj.img + obj.RADIANCE_OFFEST;
            end
        end
        
        function delete(obj)
            delete(obj.CAM_MDL);
            delete(obj.CAM_MDL_GEO);
            delete(obj.ROVER_NAV);
            delete(obj.RMC);
        end   
    end
end