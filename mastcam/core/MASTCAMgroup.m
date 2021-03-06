classdef MASTCAMgroup < handle
    % MASTCAMgroup
    %   class handling a group of MASTCAMdata that are taken from the same
    %   location. 
    %   The edit of this class is only accessible through append method.
    %   Properties
    %    L: class object of MASTCAMgroup_eye
    %    R: class object of MASTCAMgroup_eye
    %    PROCUCT_ID:, cell array of the product ids of the images loaded
    %    RMC: MSL_RMC class object, storing Rover Motion Counter
    %    information.
    
    properties
        L          % = MASTCAMgroup_eye(); % LEFT CAMERA
        R          % = MASTCAMgroup_eye(); % RIGHT CAMERA
        PRODUCT_ID                       % PRODUCT ID list
        RMC                              % ROVER MOTION COUNTER
        ROVER_NAV
    end
    
    methods
        function obj = MASTCAMgroup()
            obj.L = MASTCAMgroup_eye();
            obj.R = MASTCAMgroup_eye();
            
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
                 obj.PRODUCT_ID = {mst_obj.PRODUCT_ID};
            else
                if obj.RMC ~= mst_obj.RMC
                    error('INPUT is not in this group');
                end
                obj.PRODUCT_ID = [obj.PRODUCT_ID {mst_obj.PRODUCT_ID}];
            end
            
            % append the data to an appropirate branch
            switch upper(mst_obj.prop.cam_code)
                case 'ML'
                    mst_eye = 'L';
                case 'MR'
                    mst_eye = 'R';
            end
            obj.(mst_eye).append(mst_obj);
        end
        
        function append_AXIX(obj,basename,varargin)
            propAXIX = get_basenameMASTCAM_fromProp(basename);
            switch upper(propAXIX.cam_code)
                case 'ML'
                    obj.L.add_AXIX(basename,varargin{:});
                case 'MR'
                    obj.R.add_AXIX(basename,varargin{:});
                otherwise
                    error('Something wrong with basename %s.',basename);
            end
        end
        
        function update_ROVER_NAV(obj,varargin)
            if isempty(varargin) 
            elseif length(varargin)==1
                obj.ROVER_NAV = varargin{1};
            else
                error('Input is invalid');
            end
            obj.L.update_ROVER_NAV(obj.ROVER_NAV);
            obj.R.update_ROVER_NAV(obj.ROVER_NAV);
        end
        function update_ROVER_NAV_DEM(obj,MSLDEMdata)
            obj.ROVER_NAV.update_DEM(MSLDEMdata);
            obj.update_ROVER_NAV();
        end
        function update_ROVER_NAV_MAP(obj,MSLOrthodata)
            obj.ROVER_NAV.update_MAP(MSLOrthodata);
            obj.update_ROVER_NAV();
        end
        
        function delete(obj)
            if ~isempty(obj.L)
                % delete(obj.L);
            end
            if ~isempty(obj.R)
                % delete(obj.R);
            end
            if ~isempty(obj.RMC)
                % delete(obj.RMC);
            end
            if ~isempty(obj.ROVER_NAV)
                % delete(obj.ROVER_NAV);
            end
        end
    end
end