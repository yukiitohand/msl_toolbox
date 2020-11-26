classdef MASTCAMgroupDRLX < MASTCAMgroup_wProcCode
    % MASTCAMgroupDRLX
    %   class handling a group of MASTCAMdata that are taken from the same
    %   location and the same camera (Left of Right) and shares the same
    %   DATA PROCESSING CODE ('DRLX').
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
        
    end
    
    methods
        function obj = MASTCAMgroupDRLX()
            obj@MASTCAMgroup_wProcCode(...
                'DRLX','MEMBER_CLASS_NAME','MASTCAMdataDRLX');
        end
        
        function [imrgb] = get_rgb(obj,varargin)
            tol = 0;
            priority_dtype_list0 = {'E','C'};
            mstdata0 = find_MASTCAMdata_Filterk(obj,0,priority_dtype_list0);
            imrgb = mstdata0.readimg('datatype','IoF');
            imrgb = RGBImage(imrgb,'Tol',tol);
        end
        
    end
end