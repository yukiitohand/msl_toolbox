classdef MASTCAMMSI < dynamicprops
    % MASTCAM Multispectral Image class
    %   
    
    properties
        MASTCAMgroup_eye
        filternames
        filter_info
        wavelength_filternames
        wavelength
        DATA_PROC_CODE
        img
        img_iof
        eye
        radiance_factor
        hdr
    end
    
    methods
        function obj = MASTCAMMSI(mstdata_seq,DATA_PROC_CODE)
            switch mstdata_seq.eye
                case 'L'
                    obj.eye = 'L';
                    [obj.filter_info] = mastcamL_get_filter_info();
                    obj.filternames = {'L0Red','L0Green','L0Blue','L1','L2','L3','L4','L5','L6'};
                case 'R'
                    obj.eye = 'R';
                    [obj.filter_info] = mastcamR_get_filter_info();
                     obj.filternames = {'R0Red','R0Green','R0Blue','R1','R2','R3','R4','R5','R6'};
            end
            addprop(obj,[obj.eye,'0']);
            for i=1:length(obj.filternames)
                addprop(obj,obj.filternames{i});
            end
            
            obj.set_wavelength();

            obj.DATA_PROC_CODE = DATA_PROC_CODE;
            obj.set_img(mstdata_seq);
            if any(strcmpi(DATA_PROC_CODE,{'DRLX','DRXX'}))
                obj.set_imgIoF(mstdata_seq);
            end
            obj.img = obj.img_iof;
            obj.MASTCAMgroup_eye = mstdata_seq;
            obj.set_header();

        end
        
        function set_wavelength(obj)
            % filternames = fieldnames(obj.filter_info);
            filterwavelength = cellfun(@(x) obj.filter_info.(x).wv_ctr,obj.filternames);
            [wv_sort,wv_sort_index] = sort(filterwavelength,'ascend');
            filternames_sorted = obj.filternames(wv_sort_index);
            obj.wavelength = wv_sort;
            obj.wavelength_filternames = filternames_sorted;
        end
        
        function set_header(obj)
            obj.hdr = [];
            if isprop(obj.MASTCAMgroup_eye,'E')
                obj.hdr.samples = obj.MASTCAMgroup_eye.E.DRCL.hdr.samples;
                obj.hdr.lines = obj.MASTCAMgroup_eye.E.DRCL.hdr.lines;
            elseif isprop(obj.MASTCAMgroup_eye,'C')
                obj.hdr.samples = obj.MASTCAMgroup_eye.C.DRCL.hdr.samples;
                obj.hdr.lines = obj.MASTCAMgroup_eye.C.DRCL.hdr.lines;
            end
            obj.hdr.bands = size(obj.img,3);
        end
        
        function set_img(obj,mstdata_seq)
            % load images from mstdat_seq
            if isprop(mstdata_seq,'E')
                imRGB = mstdata_seq.E.(obj.DATA_PROC_CODE).readimg();
            elseif isprop(mstdata_seq,'C')
                imRGB = mstdata_seq.C.(obj.DATA_PROC_CODE).readimg();
            end
            obj.([obj.eye '0']) = imRGB;
            obj.([obj.eye '0Red']) = imRGB(:,:,1);
            obj.([obj.eye '0Green']) = imRGB(:,:,2);
            obj.([obj.eye '0Blue']) = imRGB(:,:,3);
            
            for i=1:length(mstdata_seq.D.(obj.DATA_PROC_CODE))
                filter_number = mstdata_seq.D.(obj.DATA_PROC_CODE)(i).FILTER_NUMBER;
                filter_id = sprintf('%1s%1d',obj.eye,filter_number);
                obj.(filter_id) = mstdata_seq.D.(obj.DATA_PROC_CODE)(i).readimg();
            end
            
            % order the image by their filter wavelengths
            for i=1:length(obj.wavelength_filternames)
                if i==1
                    obj.img = obj.(obj.wavelength_filternames{i});
                else
                    obj.img = cat(3,obj.img,obj.(obj.wavelength_filternames{i}));
                end
            end
        end
        
        function set_imgIoF(obj,mstdata_seq)
            obj.radiance_factor = [];
            iof = [];
            if isprop(mstdata_seq,'E')
                imRGB_iof = mstdata_seq.E.(obj.DATA_PROC_CODE).get_IoF();
            elseif isprop(mstdata_seq,'C')
                imRGB_iof = mstdata_seq.C.(obj.DATA_PROC_CODE).get_IoF();
            end
            iof.([obj.eye '0Red']) = imRGB_iof(:,:,1);
            iof.([obj.eye '0Green']) = imRGB_iof(:,:,2);
            iof.([obj.eye '0Blue']) = imRGB_iof(:,:,3);
            if isprop(mstdata_seq,'E')
                obj.radiance_factor.([obj.eye '0']) = mstdata_seq.E.(obj.DATA_PROC_CODE).RADIANCE_FACTOR;
            elseif isprop(mstdata_seq,'C')
                obj.radiance_factor.([obj.eye '0']) = mstdata_seq.C.(obj.DATA_PROC_CODE).RADIANCE_FACTOR;
            end
            
            for i=1:length(mstdata_seq.D.(obj.DATA_PROC_CODE))
                filter_number = mstdata_seq.D.(obj.DATA_PROC_CODE)(i).FILTER_NUMBER;
                filter_id = sprintf('%1s%1d',obj.eye,filter_number);
                iof.(filter_id) = mstdata_seq.D.(obj.DATA_PROC_CODE)(i).get_IoF();
                obj.radiance_factor.(filter_id) = mstdata_seq.D.(obj.DATA_PROC_CODE)(i).RADIANCE_FACTOR;
            end
            
            % order the image by their filter wavelengths
            for i=1:length(obj.wavelength_filternames)
                if i==1
                    obj.img_iof = iof.(obj.wavelength_filternames{i});
                else
                    obj.img_iof = cat(3,obj.img_iof,iof.(obj.wavelength_filternames{i}));
                end
            end
        end
        
        function [spc,wv,bdxes] = get_spectrum(obj,s,l,varargin)
            spc = squeeze(obj.img(l,s,:));
            wv = obj.wavelength;
            bdxes = 1:obj.hdr.bands;
        end
        
    end
end