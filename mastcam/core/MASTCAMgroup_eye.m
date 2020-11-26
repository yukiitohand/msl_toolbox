classdef MASTCAMgroup_eye < dynamicprops
    % MASTCAMgroup_eye
    %   class handling a group of MASTCAMdata that are taken from the same
    %   location and the same camera (Left of Right).
    %   Images are stored under the tree of two layers. The first layer is
    %   categorized by their PRODUCT_TYPE, which is one of the characters 
    %   (A-Z), and the second layer is by their DATA_PROCESSING_CODE, which
    %   consists of four characters among (D,R,C,L,X), depending on the
    %   processing levels. This is implemented using a dynamic assignment
    %   of properties to a class object.
    %   Some common properties are stored in the most upper layer.
    %    eye: {'L' or 'R'}
    %    PROCUCT_ID:, cell array of the product ids of the images loaded
    %    RMC: MSL_RMC class object, storing Rover Motion Counter
    %    information.
    %    CAM_MDL: CAHVOR_MODEL class object, storing camera models.
    %
    %   The edit of this class is only accessible through append method.
    
    properties
        eye         % eye of the group
        PRODUCT_ID  % PRODUCT ID list
        RMC         % ROVER MOTION COUNTER
        CAM_MDL     % CAMERA MODEL
        ROVER_NAV   % ROVER NAVIGATION FRAME
        CAM_MDL_GEO % CAMERA MODEL in a geographic coordinate system
        L_im
        S_im
        DRXX
        DRCX
        DRLX
        DRCL
        AXI1
        addedProps
    end
    
    methods
        function obj = MASTCAMgroup_eye()
            obj.DRXX = MASTCAMgroupDRXX();
            obj.DRCX = MASTCAMgroupDRCX();
            obj.DRLX = MASTCAMgroupDRLX();
            obj.DRCL = MASTCAMgroupDRCL();
            obj.AXI1 = MASTCAMgroupAXIX(1);
            
            % obj.PRODUCT_ID = struct('DRXX',[],'DRCX',[],'DRLX',[],...
            %     'DRCL',[],'AXI1',[]);
        end
        
        function append(obj,mst_obj)
            % input check
            if ~isa(mst_obj,'MASTCAMdata') && ~isa(mst_obj,'MASTCAMdataAXIX')
                error('INPUT must be an object of MASTCAMdata class or MASTCAMdataAXIX');
            end
            % check the data is appropriate or not first
            if isempty(obj.PRODUCT_ID)
                obj.RMC = mst_obj.RMC;
                obj.ROVER_NAV = mst_obj.ROVER_NAV;
                obj.CAM_MDL  = mst_obj.CAM_MDL;
                obj.PRODUCT_ID = {mst_obj.PRODUCT_ID};
                obj.L_im = mst_obj.hdr.lines;
                obj.S_im = mst_obj.hdr.samples;
                switch upper(mst_obj.prop.cam_code)
                    case 'ML'
                        obj.eye = 'L';
                    case 'MR'
                        obj.eye = 'R';
                end  
            else
                if obj.RMC ~= mst_obj.RMC %% || obj.CAM_MDL ~= mst_obj.CAM_MDL
                    error('INPUT is not in this group');
                end
                obj.PRODUCT_ID = [obj.PRODUCT_ID {mst_obj.PRODUCT_ID}];
            end
            
            % append the data to an appropirate branch
            % add a dynamic property if the 
            data_proc_code = mst_obj.prop.data_proc_code;
            if isprop(obj,data_proc_code)
                obj.(data_proc_code).append(mst_obj);
            else
                error('DATA_PROC_CODE %s is undefined.',data_proc_code);
            end

        end
        
        function load_AXIX(obj,varargin)
            vr_AXIX = 1;
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case 'VEERSION'
                            vr_AXIX = varargin{i+1};
                    end
                end
            end
            if isnumeric(vr_AXIX), vr_AXIX = num2str(vr_AXIX,'%1d'); end
            propname = sprintf('AXI%s',vr_AXIX);
            if ~isempty(obj.DRXX)
                dtype_list = {'E','D','C'};
                for j=1:length(dtype_list)
                    dtypej = dtype_list{j};
                    if isprop(obj.DRXX,dtypej) && ~isempty(obj.DRXX.(dtypej))
                        for i=1:length(obj.DRXX.(dtypej))
                            [mstaxixdata] = mastcam_get_AXIXdata_from_ref(obj.DRXX.(dtypej)(i),varargin{:});
                            obj.(propname).append(mstaxixdata);
                        end
                    end
                end
            end
        end

        function delete(obj)
            if ~isempty(obj.RMC)
                delete(obj.RMC);
            end
            if ~isempty(obj.ROVER_NAV)
                delete(obj.ROVER_NAV);
            end
            if ~isempty(obj.CAM_MDL)
                delete(obj.CAM_MDL);
            end
            delete(obj.DRXX);
            delete(obj.DRCX);
            delete(obj.DRLX);
            delete(obj.DRCL);
            delete(obj.AXI1);
            for i=1:length(obj.addedProps)
                propi = obj.addedProps{i};
                delete(obj.(propi));
            end
            
        end
        
        function [mstmsi] = MASTCAMMSIConstructor(obj,data_proc_code,varargin)
            priority_dtype_list0 = {'E','C'};
            switch upper(data_proc_code)
                case {'DRXX','AXI1'}
                    mstdata_ColorCor = find_MASTCAMdata_Filterk(obj.DRCX,0,priority_dtype_list0);
                case {'DRLX'}
                    mstdata_ColorCor = find_MASTCAMdata_Filterk(obj.DRCL,0,priority_dtype_list0);
                otherwise
                    mstdata_ColorCor = [];
            end
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case {'MASTCAMCOLORCOR'}
                            mstdata_ColorCor = varargin{i+1};
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            
            switch upper(data_proc_code)
                case {'DRXX','DRLX','AXI1'}
                    mstmsi = MASTCAMMSI(obj.(data_proc_code),...
                        'MASTCAMColorCor',mstdata_ColorCor);
                case {'DRCL','DRCX'}
                    error('Does not make sense to contruct MSI for color corrected images.');
                otherwise
                    error('Undefined DATA_PROC_CODE %s.',data_proc_code);
            end
            
        end
        
%         function [] = get_CAM_MDL_GEO(obj)
%             [imxy_direc_rov] = get_3d_pointing_from_CAHV_v2(...
%                 [obj.L_im,obj.S_im],obj.CAM_MDL);
%             cmmdl_A_rov0 = obj.ROVER_NAV.rot_mat * obj.CAM_MDL.A';
%             cmmdl_C_rov0 = obj.ROVER_NAV.rot_mat * obj.CAM_MDL.C';
%             imxy_direc_rov_2d = reshape(imxy_direc_rov,[obj.L_im*obj.S_im,3])';
%             imxy_direc_rov0_2d = obj.ROVER_NAV.rot_mat * imxy_direc_rov_2d;
%             imxy_direc_rov0 = reshape(imxy_direc_rov0_2d',[obj.L_im,obj.S_im,3]);
%             cmmdl_C_geo = cmmdl_C_rov0 + [obj.ROVER_NAV.NORTHING; 
%                               obj.ROVER_NAV.EASTING;
%                               -obj.ROVER_NAV.ELEVATION];
%                           
%             obj.CAM_MDL_GEO.A_rov0 = cmmdl_A_rov0';
%             obj.CAM_MDL_GEO.C_rov0 = cmmdl_C_rov0';
%             obj.CAM_MDL_GEO.C_geo  = cmmdl_C_geo';
%             obj.CAM_MDL_GEO.imxy_direc_rov = imxy_direc_rov;
%             obj.CAM_MDL_GEO.imxy_direc_rov0 = imxy_direc_rov0;
%             
%         end
        
%         function append_AXIX(obj,basename,varargin)
%             mstaidata_obj = mastcam_get_AXIXdata(basename,obj,varargin{:});
%             obj.append(mstaidata_obj);
%             % fpath = joinPath(dirpath,[basename '.IMG']);
%             % propA6I1 = get_basenameMASTCAM_fromProp(basename);
% 
%         end
        
    end
end