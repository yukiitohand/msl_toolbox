classdef MASTCAMMSIview < handle
    %MASTCAMview
    %   Spectral Viewer for MASTCAM images
    properties
        obj_HSIview
        MASTCAMdataseq_eye
        MSTMSI
    end
    
    methods
        function obj = MASTCAMMSIview(mstmsi,varargin)
            %Constructor
            obj.MSTMSI = mstmsi;
            [imRGB] = mstmsi.get_rgb();
            obj.obj_HSIview = HSIview({imRGB},...
            {{obj.MSTMSI}},...
            ...'SPC_XLIM',[300 1200],...
            'varargin_ImageStackView',{'Ydir','reverse','XY_COORDINATE_SYSTEM','IMAGEPIXELS'});
        end
        
        
        
    end
end