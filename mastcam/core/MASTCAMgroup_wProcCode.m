classdef MASTCAMgroup_wProcCode < dynamicprops
    % MASTCAMgroup_wProcCode
    %   class handling a group of MASTCAMdata that are taken from the same
    %   location and the same camera (Left of Right) and shares the same
    %   DATA PROCESSING CODE.
    %   Images are stored categorized by its data_type ('C','D','E',...) 
    %   This is implemented using a dynamic assignment of properties to a 
    %   class object.
    %   Some common properties are stored in the most upper layer.
    %    eye: {'L' or 'R'}
    %    PROCUCT_ID: struct, the product ids of the images loaded,
    %    organized by data_type
    %    RMC: MSL_RMC class object, storing Rover Motion Counter
    %    information.
    %    CAM_MDL: CAHVOR_MODEL class object, storing camera models.
    %
    %   The edit of this class is only accessible through append method.
    
    properties
        eye         % eye of the group
        PRODUCT_ID  % PRODUCT ID list
        DATA_PROC_CODE % DATA PROCESSING CODE
        RMC         % ROVER MOTION COUNTER
        CAM_MDL     % CAMERA MODEL
        ROVER_NAV   % ROVER NAVIGATION FRAME
        CAM_MDL_GEO % CAMERA MODEL in a geographic coordinate system
        L_im
        S_im
        Linearization
        addedProps
        member_class_name
        C
        D
        E
    end
    
    methods
        function obj = MASTCAMgroup_wProcCode(data_proc_code,varargin)
            obj.DATA_PROC_CODE    = data_proc_code;
            obj.member_class_name = {'MASTCAMdata','MASTCAMdataAXIX'};
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case 'MEMBER_CLASS_NAME'
                            obj.member_class_name = varargin{i+1};
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            
            
            
        end
        
        function append(obj,mst_obj)
            % input check
            if ~obj.mstdata_ismember(mst_obj)
                if iscell(obj.member_class_name)
                    clnm = strjoin(obj.member_class_name,',');
                else
                    clnm = obj.member_class_name;
                end
                error('INPUT must be an object of class %s',clnm);
            end
            % check the data is appropriate or not first
            if ~obj.compare_DATA_PROC_CODE(mst_obj)
                error('Processing Code should be %s. (your input has %s)',...
                    obj.DATA_PROC_CODE,mst_obj.prop.data_proc_code);
            end
            if isempty(obj.PRODUCT_ID)
                
                obj.RMC = mst_obj.RMC;
                obj.ROVER_NAV = mst_obj.ROVER_NAV;
                obj.CAM_MDL  = mst_obj.CAM_MDL;
                obj.PRODUCT_ID = {mst_obj.PRODUCT_ID};
                obj.L_im = mst_obj.hdr.lines;
                obj.S_im = mst_obj.hdr.samples;
                obj.DATA_PROC_CODE = mst_obj.prop.data_proc_code;
                switch upper(mst_obj.prop.cam_code)
                    case 'ML'
                        obj.eye = 'L';
                    case 'MR'
                        obj.eye = 'R';
                end  
            else
                if obj.RMC ~= mst_obj.RMC
                    error('INPUT is not in this group');
                end
                obj.PRODUCT_ID = [obj.PRODUCT_ID {mst_obj.PRODUCT_ID}];
            end
            
            % append the data to an appropirate branch
            % add a dynamic property if the 
            product_type = mst_obj.prop.product_type;
            if ~isprop(obj,product_type) 
                addprop(obj,product_type);
                obj.addedProps = [obj.addedProps {product_type}];
            end
                
            if isempty(obj.(product_type))
                obj.(product_type) = mst_obj;
                % obj.(product_type) = mst_obj.PRODUCT_ID;
                % obj.PRODUCT_ID.(product_type) = mst_obj.PRODUCT_ID;
            else
                obj.(product_type) = [obj.(product_type) mst_obj];
                % obj.PRODUCT_ID.(product_type) = ...
                %         [obj.PRODUCT_ID.(product_type) {mst_obj.PRODUCT_ID}];
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
        
        function [tf] = compare_DATA_PROC_CODE(obj,mst_obj)
            tf = strcmpi(obj.DATA_PROC_CODE,mst_obj.prop.data_proc_code);
        end
        
        function [tf] = mstdata_ismember(obj,mst_obj)
            if iscell(obj.member_class_name)
                tf = any(cellfun(@(x) isa(mst_obj,x), obj.member_class_name));
            else
                tf = isa(mst_obj,obj.member_class_name);
            end
        end
        
        function [tf] = isSameScene(mstgrp_wpc)
            if strcmpi(obj.eye,mstgrp_wpc.eye) && obj.RMC == mstgrp_wpc.RMC ...
                && obj.Linearization == mstgrp_wpc.Linearization
                tf = true;
            else
                tf = false;
            end
        end
        
        function update_ROVER_NAV(obj,varargin)
            if isempty(varargin) 
            elseif length(varargin)==1
                obj.ROVER_NAV = varargin{1};
            else
                error('Input is invalid');
            end
            if ~isempty(obj.C)
                for k=1:length(obj.C)
                    obj.C(k).update_ROVER_NAV(obj.ROVER_NAV);
                end
            end
            if ~isempty(obj.D)
                for k=1:length(obj.D)
                    obj.D(k).update_ROVER_NAV(obj.ROVER_NAV);
                end
            end
            if ~isempty(obj.E)
                for k=1:length(obj.E)
                    obj.E(k).update_ROVER_NAV(obj.ROVER_NAV);
                end
            end
            for i=1:length(obj.addedProps)
                propi = obj.addedProps{i};
                for k=1:length(obj.(propi))
                    obj.(propi)(k).update_ROVER_NAV(obj.ROVER_NAV);
                end
            end
        end
        
        function update_ROVER_NAV_DEM(obj,MSLDEMdata)
            obj.ROVER_NAV.update_DEM(MSLDEMdata);
            obj.update_ROVER_NAV();
        end
        
        function update_ROVER_NAV_MAP(obj,MSLOrthodata)
            obj.ROVER_NAV.update_DEM(MSLOrthodata);
            obj.update_ROVER_NAV();
        end
        
        function delete(obj)
            if ~isempty(obj.RMC)
                % delete(obj.RMC);
            end
            if ~isempty(obj.ROVER_NAV)
                % delete(obj.ROVER_NAV);
            end
            if ~isempty(obj.CAM_MDL)
                % delete(obj.CAM_MDL);
            end
            % if ~isempty(obj.CAM_MDL_GEO) && isvalid(obj.CAM_MDL_GEO)
            %     delete(obj.CAM_MDL_GEO);
            % end
            if ~isempty(obj.C)
                for k=1:length(obj.C)
                    delete(obj.C(k));
                end
            end
            if ~isempty(obj.D)
                for k=1:length(obj.D)
                    delete(obj.D(k));
                end
            end
            if ~isempty(obj.E)
                for k=1:length(obj.E)
                    delete(obj.E(k));
                end
            end
            for i=1:length(obj.addedProps)
                propi = obj.addedProps{i};
                for k=1:length(obj.(propi))
                    delete(obj.(propi)(k));
                end
            end
            
        end
        
    end
end