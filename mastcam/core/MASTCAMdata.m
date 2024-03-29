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
        FILTER_NUMBER
        INSTRUMENT_ID
        L_im
        S_im
        EYE
        Linearization
    end
    
    methods
        function obj = MASTCAMdata(basename,dirpath,varargin)
            
            rover_nav_ver = 'localized_interp';
            rover_nav_mstcode = '';
            rover_nav_lr = [];
            varargin_rmidx = [];
            lr = [];
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
                        case 'ROVER_NAV_LINEARIZATION'
                            rover_nav_lr = varargin{i+1};
                            varargin_rmidx = [varargin_rmidx i i+1];
                        case 'LINEARIZATION'
                            lr = varargin{i+1};
                            varargin_rmidx = [varargin_rmidx i i+1];
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            varargin_HSI = varargin(setdiff(1:length(varargin),varargin_rmidx));
            obj@HSI(basename,dirpath,varargin_HSI{:});
            obj.lblpath = joinPath(dirpath,[basename '.LBL']);
            obj.lbl = pds3lblread(obj.lblpath);
            obj.hdr = mastcam_extract_imghdr_from_lbl(obj.lbl);
            obj.prop = getProp_basenameMASTCAM(basename);
            obj.Linearization = lr;
            obj.CAM_MDL = get_cammera_model(obj);
            obj.RMC = get_rmc(obj);
            obj.PRODUCT_ID = obj.lbl.PRODUCT_ID;
            obj.ROVER_NAV  = obj.get_rover_nav('VERSION',rover_nav_ver,...
                'MSTCAM_CODE',rover_nav_mstcode,'LINEARIZATION',rover_nav_lr);
            obj.get_filter_number();
            obj.get_instrument_id();
            obj.L_im = obj.hdr.lines;
            obj.S_im = obj.hdr.samples;
        end
        function [cahvor_mdl] = get_cammera_model(obj)
            [cahvor_mdl] = mastcam_get_cahvor_model(obj.lbl);
        end
        function [rmc_mdl] = get_rmc(obj)
            [rmc_mdl] = mastcam_get_RMC(obj.lbl);
        end
        function [rover_nav] = get_rover_nav(obj,varargin)
            [rover_nav,option] = msl_get_ROVER_NAV(obj,varargin{:});
            % [rover_nav] = get_ROVER_NAV_from_mslplc_view_csv(...
            %     obj.RMC.SITE,obj.RMC.DRIVE,obj.RMC.POSE,varargin{:});
        end
        
        function get_instrument_id(obj)
            obj.INSTRUMENT_ID = obj.lbl.INSTRUMENT_ID;
        end
        
        function get_filter_number(obj)
            obj.FILTER_NUMBER = obj.lbl.GROUP_IMAGE_REQUEST_PARMS.FILTER_NUMBER;
        end
        
        function update_ROVER_NAV(obj,rover_nav_new)
           obj.ROVER_NAV = rover_nav_new;
        end
        function update_ROVER_NAV_DEM(obj,MSLDEMdata)
           obj.ROVER_NAV.update_DEM(MSLDEMdata);
        end
        function update_ROVER_NAV_MAP(obj,MSLOrthodata)
           obj.ROVER_NAV.update_DEM(MSLOrthodata);
        end
        
        function delete(obj)
            % delete(obj.CAM_MDL);
            % delete(obj.CAM_MDL_GEO);
            % delete(obj.ROVER_NAV);
            % delete(obj.RMC);
        end   
    end
end

