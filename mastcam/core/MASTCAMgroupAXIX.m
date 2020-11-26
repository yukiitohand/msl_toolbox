classdef MASTCAMgroupAXIX < MASTCAMgroup_wProcCode
    % MASTCAMgroupAXIX
    %   class handling a group of MASTCAMdata that are taken from the same
    %   location and the same camera (Left of Right) and shares the same
    %   DATA PROCESSING CODE ('AXIX').
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
        DATA_PROC_CODE_PTRN % DATA PROCESSING CODE
    end
    
    methods
        function obj = MASTCAMgroupAXIX(vr)
            if isnumeric(vr)
                vr = num2str(vr,'%1d');
            end
            data_proc_code = sprintf('AXI%s',vr);
            data_proc_code_ptrn = sprintf('A\\dI%s',vr);
            obj@MASTCAMgroup_wProcCode(...
                data_proc_code,'MEMBER_CLASS_NAME','MASTCAMdataAXIX');
            obj.DATA_PROC_CODE_PTRN = data_proc_code_ptrn;
        end
        
        function [tf] = compare_DATA_PROC_CODE(obj,mst_obj)
            tf = ~isempty(regexpi(mst_obj.prop.data_proc_code,obj.DATA_PROC_CODE_PTRN,'once'));
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