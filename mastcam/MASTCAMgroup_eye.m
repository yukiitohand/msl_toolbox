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
        addedProps
    end
    
    methods
        function obj = MASTCAMgroup_eye()
            
        end
        
        function append(obj,mst_obj)
            % input check
            if ~isa(mst_obj,'MASTCAMdata')
                error('INPUT must be an object of MASTCAMdata class');
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
            if ~isprop(obj,mst_obj.prop.product_type)
                addprop(obj,mst_obj.prop.product_type);
                obj.(mst_obj.prop.product_type).(mst_obj.prop.data_proc_code) = mst_obj;
                obj.addedProps = [obj.addedProps {mst_obj.prop.product_type}];
            else
                if isfield(obj.(mst_obj.prop.product_type),mst_obj.prop.data_proc_code)
                    obj.(mst_obj.prop.product_type).(mst_obj.prop.data_proc_code)...
                        = [obj.(mst_obj.prop.product_type).(mst_obj.prop.data_proc_code) mst_obj];
                else
                    obj.(mst_obj.prop.product_type).(mst_obj.prop.data_proc_code) = mst_obj;
                end
            end

        end
        function [] = get_CAM_MDL_GEO(obj)
            [imxy_direc_rov] = get_3d_pointing_from_CAHV_v2(...
                [obj.L_im,obj.S_im],obj.CAM_MDL);
            cmmdl_A_rov0 = obj.ROVER_NAV.rot_mat * obj.CAM_MDL.A';
            cmmdl_C_rov0 = obj.ROVER_NAV.rot_mat * obj.CAM_MDL.C';
            imxy_direc_rov_2d = reshape(imxy_direc_rov,[obj.L_im*obj.S_im,3])';
            imxy_direc_rov0_2d = obj.ROVER_NAV.rot_mat * imxy_direc_rov_2d;
            imxy_direc_rov0 = reshape(imxy_direc_rov0_2d',[obj.L_im,obj.S_im,3]);
            cmmdl_C_geo = cmmdl_C_rov0 + [obj.ROVER_NAV.NORTHING; 
                              obj.ROVER_NAV.EASTING;
                              -obj.ROVER_NAV.ELEVATION];
                          
            obj.CAM_MDL_GEO.A_rov0 = cmmdl_A_rov0';
            obj.CAM_MDL_GEO.C_rov0 = cmmdl_C_rov0';
            obj.CAM_MDL_GEO.C_geo  = cmmdl_C_geo';
            obj.CAM_MDL_GEO.imxy_direc_rov = imxy_direc_rov;
            obj.CAM_MDL_GEO.imxy_direc_rov0 = imxy_direc_rov0;
            
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
            % if ~isempty(obj.CAM_MDL_GEO) && isvalid(obj.CAM_MDL_GEO)
            %     delete(obj.CAM_MDL_GEO);
            % end
            for i=1:length(obj.addedProps)
                propi = obj.addedProps{i};
                fldnms = fieldnames(obj.(propi));
                for j=1:length(fldnms)
                    fldnm = fldnms{j};
                    for k=1:length(obj.(propi).(fldnm))
                        delete(obj.(propi).(fldnm)(k));
                    end
                end
            end
            
        end
        
    end
end