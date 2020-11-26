classdef MASTCAMMSI < dynamicprops
    % MASTCAM Multispectral Image class
    %   
    
    properties
        MASTCAMgroup_wpc
        % filternames
        filter_info
        wavelength
        wavelength_filternames
        wavelength2filter
        img
        hdr
        MASTCAMColorCor
    end
    
    methods
        function obj = MASTCAMMSI(mstgrp_wpc,varargin)
            
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case {'MASTCAMCOLORCOR'}
                            obj.MASTCAMColorCor = varargin{i+1};
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            
            % decide which products are used for constructing the
            % multispectral image cube.
            
            % find filter 0
            priority_dtype_list0 = {'E','C'};
            mstdata0 = find_MASTCAMdata_Filterk(mstgrp_wpc,0,priority_dtype_list0);
            if ~isempty(mstdata0)
                propName0 = [mstgrp_wpc.eye,'0'];
                addprop(obj,propName0);
                obj.(propName0) = mstdata0;
                % propNameList = [propNameList {propName0}];
            end
            
            % find filter 1-6
            priority_dtype_list = {'D','C'};
            for k=1:6
                mstdatak = find_MASTCAMdata_Filterk(mstgrp_wpc,k,priority_dtype_list);
                if ~isempty(mstdatak)
                    propNamek = [mstgrp_wpc.eye,num2str(k,'%1d')];
                    % fprintf('%s,',propNamek);
                    addprop(obj,propNamek);
                    obj.(propNamek) = mstdatak;
                    % propNameList = [propNameList {propNamek}];
                end
            end
            
            % create a filter info struct
            [filter_info_all] = mastcam_get_filter_info(mstgrp_wpc.eye);
            for k=0:6
                propNamek = [mstgrp_wpc.eye,num2str(k,'%1d')];
                if isprop(obj,propNamek) && ~isempty(obj.(propNamek))
                    obj.filter_info.(propNamek) = filter_info_all.(propNamek);
                end
            end
            
            obj.set_wavelength();
            obj.MASTCAMgroup_wpc = mstgrp_wpc;
            obj.set_header();
            obj.readimg();

        end
        
        function set_wavelength(obj)
            % filternames = fieldnames(obj.filter_info);
            wv = [];
            fltids = fieldnames(obj.filter_info);
            nFlt = length(fltids);
            wv2filter = []; wv2filind = [];
            cumIdx = 1;
            for i=1:nFlt
                fltid = fltids{i};
                ni = length(obj.filter_info.(fltid).wv_ctr);
                wv = [wv obj.filter_info.(fltid).wv_ctr];
                wv2filter = [wv2filter repmat({fltid},[1 ni])];
                wv2filind = [wv2filind 1:ni];
                cumIdx_i = cumIdx:(cumIdx+ni-1);
                obj.filter_info.(fltid).cumIdx = cumIdx_i;
                % update cumIdx
                cumIdx = cumIdx+ni;
            end
            B = length(wv);
            [wv_sort,wv_sort_index] = sort(wv,'ascend');
            [~,wv_sort_index_r] = sort(wv_sort_index,'ascend');
            for i=1:nFlt
                fltid = fltids{i};
                obj.filter_info.(fltid).sortIdx = wv_sort_index_r(obj.filter_info.(fltid).cumIdx);
            end
            
            % mapping to filter 
            wv_sort2filter = wv2filter(wv_sort_index);
            wv_sort2filind = wv2filind(wv_sort_index);
            
            % construct filter names
            obj.wavelength_filternames = cell(1,B);
            for i=1:B
                if strcmpi(wv_sort2filter{i}(2),'0')
                    if wv_sort2filind(i)==1
                        fltname = [wv_sort2filter{i} 'Red'];
                    elseif wv_sort2filind(i)==2
                        fltname = [wv_sort2filter{i} 'Green'];
                    elseif wv_sort2filind(i)==3
                        fltname = [wv_sort2filter{i} 'Blue'];
                    end
                else
                    fltname = wv_sort2filter{i};
                end
                obj.wavelength_filternames{i} = fltname;
            end
            
            obj.wavelength = wv_sort;
            
            obj.wavelength2filter = struct(...
                'FILTER_ID', wv_sort2filter,...
                'FILTER_SUBIND', num2cell(wv_sort2filind));
            
        end
        
        function set_header(obj)
            obj.hdr = [];
            obj.hdr.samples = obj.MASTCAMgroup_wpc.S_im;
            obj.hdr.lines = obj.MASTCAMgroup_wpc.L_im;
            obj.hdr.bands = length(obj.wavelength);
        end
        
        function readimg(obj)
            % load images from mstdat_seq
            obj.img = nan(obj.hdr.lines,obj.hdr.samples,obj.hdr.bands);
            fltids = fieldnames(obj.filter_info);
            nFlt = length(fltids);
            for i=1:nFlt
                fltid = fltids{i};
                img_flt = obj.(fltid).readimg();
                obj.img(:,:,obj.filter_info.(fltid).sortIdx) = img_flt; 
            end
        end
        
        function [spc,wv,bdxes] = get_spectrum(obj,s,l,varargin)
            spc = squeeze(obj.img(l,s,:));
            wv = obj.wavelength;
            bdxes = 1:obj.hdr.bands;
        end
        
        function [imrgb] = get_rgb(obj,varargin)
            if isa(obj.MASTCAMColorCor,'MASTCAMdataDRCX') || isa(obj.MASTCAMColorCor,'MASTCAMdataDRCL')
                imrgb = obj.MASTCAMColorCor.readimg('datatype','uint8');
            elseif isa(obj.MASTCAMColorCor,'MASTCAMdataDRXX') || isa(obj.MASTCAMColorCor,'MASTCAMdataDRLX') || isa(obj.MASTCAMColorCor,'MASTCAMdataAXIX')
                imrgb = obj.MASTCAMColorCor.readimg('datatype','IoF');
            else
                imrgb = [];
            end
        end
        
    end
end
