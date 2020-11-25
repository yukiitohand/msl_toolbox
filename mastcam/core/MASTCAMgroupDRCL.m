classdef MASTCAMgroupDRCL < MASTCAMgroup_wProcCode
    % MASTCAMgroupDRCL
    %   class handling a group of MASTCAMdata that are taken from the same
    %   location and the same camera (Left of Right) and shares the same
    %   DATA PROCESSING CODE ('DRCL').
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
        function obj = MASTCAMgroupDRCL()
            obj@MASTCAMgroup_wProcCode(...
                'DRCL','MEMBER_CLASS_NAME','MASTCAMdataDRCL');
        end
        
        function [imrgb] = get_rgb(obj)
            imrgb = obj.E.readimg();
            imrgb = uint8(imrgb);
        end
        
    end
end