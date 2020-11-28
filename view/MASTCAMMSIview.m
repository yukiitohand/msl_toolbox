classdef MASTCAMMSIview < handle
    %MASTCAMview
    %   Spectral Viewer for MASTCAM images
    properties
        obj_HSIview
        MASTCAMdataseq_eye
        MSTMSI
    end
    
    methods
        function obj = MASTCAMMSIview(mstmsi_input,varargin)
            objSpecView = [];
            % varargin_rmidx = [];
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case {'SPECVIEW'}
                            objSpecView = varargin{i+1};
                            % varargin_rmidx = [varargin_rmidx i i+1];
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            % varargin_HSI = varargin(setdiff(1:length(varargin),varargin_rmidx));
            %Constructor
            % obj.MSTMSI = mstmsi;
            
            
            % for i=1:length(hsiar_input)
            [rgbList] = obj.process_mstmsi_input(mstmsi_input);
            if ~iscell(rgbList),rgbList = {rgbList}; end
            obj.obj_HSIview = HSIview(rgbList,...
            mstmsi_input,...
            ...'SPC_XLIM',[300 1200],...
            'varargin_ImageStackView',{'Ydir','reverse','XY_COORDINATE_SYSTEM','IMAGEPIXELS'},...
            'SpecView',objSpecView);
            end
        
            function [rgbList] = process_mstmsi_input(obj,mstmsi_input)
                
                if isempty(mstmsi_input)
                    rgbList = [];
                    % skip if input hsiar is empty.
                elseif isa(mstmsi_input,'MASTCAMMSI')
                    rgbList = mstmsi_input.get_rgb();
                elseif iscell(mstmsi_input)
                    if length(mstmsi_input)==1
                        if isa(mstmsi_input{1},'MASTCAMMSI')
                            rgbList = mstmsi_input{1}.get_rgb();
                        elseif iscell(mstmsi_input{1}) && isa(mstmsi_input{1}{1},'MASTCAMMSI')
                            rgbList = mstmsi_input{1}{1}.get_rgb();
                        else
                            error('Input is not proper.');
                        end
                    elseif length(mstmsi_input)>1
                        if isa(mstmsi_input{1},'MASTCAMMSI') && ~(isa(mstmsi_input{2},'MASTCAMMSI') || iscell(mstmsi_input{2}))
                            rgbList = mstmsi_input{1}.get_rgb();
                        else
                            nelem = length(mstmsi_input);
                            rgbList = cell(1,nelem);
                            for i=1:nelem
                                if iscell(mstmsi_input{i})
                                    rgbList{i} = mstmsi_input{i}{1}.get_rgb();
                                elseif isa(mstmsi_input{i},'MASTCAMMSI')
                                    rgbList{i} = mstmsi_input{i}.get_rgb();
                                else
                                    error('Input mstmsi is not proper.');
                                end
                            end
                        end
                    end
                end
            end
        
        
        
        
        
    end
end