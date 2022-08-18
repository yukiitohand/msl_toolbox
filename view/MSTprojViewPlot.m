classdef MSTprojViewPlot < handle
    % MSTprojViewPlot class
    %   Internal class for handling the plot. Mainly used for linking the
    %   spectral plot with image cursors in the ImageStackView
    
    properties
        cursor_obj
        line_obj
        im_obj
        imObj_msldemc_UFOVmask_hsiresol
        imObj_msldemc_UFOVmask_msldemresol
        imObj_mastcam_hsiresol
        imObj_mastcam_msldemresol
        cursor_msldemc_derived
        cursor_mastcam_derived
        lnObj_spc_mastcam_hsiresol
        lnObj_spc_mastcam_msldemresol
        lnObj_spc_hsi_mastcam_resol
    end
    
    methods
        function obj = MSTprojViewPlot()
        end
        
        function add_imobj(obj,imobj)
            if isempty(obj.im_obj)
                obj.im_obj = imobj;
            else
                obj.im_obj = [obj.im_obj imobj];
            end
        end
        
        function add_lineobj(obj,lineobj)
            if isempty(obj.line_obj)
                obj.line_obj = lineobj;
            else
                obj.line_obj = [obj.line_obj lineobj];
            end
        end
        
        function delete(obj)
            for i=1:length(obj.line_obj)
                delete(obj.line_obj(i));
            end
            for i=1:length(obj.im_obj)
                delete(obj.im_obj(i));
            end
            delete(obj.imObj_msldemc_UFOVmask_hsiresol);
            delete(obj.imObj_msldemc_UFOVmask_msldemresol);
            delete(obj.imObj_mastcam_hsiresol);
            delete(obj.imObj_mastcam_msldemresol);
            delete(obj.cursor_msldemc_derived);
            delete(obj.cursor_mastcam_derived);
            delete(obj.lnObj_spc_mastcam_hsiresol);
            delete(obj.lnObj_spc_mastcam_msldemresol);
            delete(obj.lnObj_spc_hsi_mastcam_resol);
        end
    end
end
