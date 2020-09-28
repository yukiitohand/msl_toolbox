classdef MASTCAMview < handle
    %MASTCAMview
    %   Spectral Viewer for MASTCAM images
    properties
        obj_HSIview
        MASTCAMdataseq_eye
        MSTMSI
    end
    
    methods
        function obj = MASTCAMview(mstdataseq_eye,varargin)
            %Constructor
            DATA_PROC_CODE = 'DRLX';
            obj.MSTMSI = MASTCAMMSI(mstdataseq_eye,DATA_PROC_CODE);
            obj.MASTCAMdataseq_eye = mstdataseq_eye;
            [imRGB] = obj.get_rgb();
            obj.obj_HSIview = HSIview({imRGB},...
            {{obj.MSTMSI}},...
            ...'SPC_XLIM',[300 1200],...
            'varargin_ImageStackView',{'Ydir','reverse','XY_COORDINATE_SYSTEM','IMAGEPIXELS'});
        end
        
        function [imRGB] = get_rgb(obj)
            if isprop(obj.MASTCAMdataseq_eye,'E') && isfield(obj.MASTCAMdataseq_eye.E,'DRCL')
                imRGB = uint8(obj.MASTCAMdataseq_eye.E.DRCL.readimg());
            elseif isprop(obj.MASTCAMdataseq_eye,'C') && isfield(obj.MASTCAMdataseq_eye.C,'DRCL')
                imRGB = uint8(obj.MASTCAMdataseq_eye.C.DRCL.readimg());
            else
                error('No C or E images are present');
            end
        end
        
        
    end
end