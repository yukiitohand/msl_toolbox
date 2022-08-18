classdef MASTCAMImages < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        filternames
        filternamesL
        filternamesR
        filter_info
        wavelength_filternames
        wavelength
        img
        eye_tags
        isL
        isR
        L0Red
        L0Green
        L0Blue
        L1
        L2
        L3
        L4
        L5
        L6
        L5Red
        L5Green
        L5Blue
        L6Red
        L6Green
        L6Blue
        R0Red
        R0Green
        R0Blue
        R1
        R2
        R3
        R4
        R5
        R6
        R4Red
        R4Green
        R4Blue
        R5Red
        R5Green
        R5Blue
        R6Red
        R6Green
        R6Blue
    end
    
    methods
        function obj = MASTCAMImages()
            [obj.filter_info] = mastcam_get_filter_info();
            obj.filternamesL = {'L0Red','L0Green','L0Blue','L1','L2','L3','L4','L5','L6'};
            obj.filternamesR = {'R0Red','R0Green','R0Blue','R1','R2','R3','R4','R5','R6'};
            obj.filternames = [obj.filternamesL,obj.filternamesR];
            obj.set_wavelength();
            % obj.set_img();
        end
        function [] = set_wavelength(obj)
            filternames = fieldnames(obj.filter_info);
            filterwavelength = cellfun(@(x) obj.filter_info.(x).wv_ctr,filternames);
            [wv_sort,wv_sort_index] = sort(filterwavelength,'ascend');
            filternames_sorted = filternames(wv_sort_index);
            obj.wavelength = wv_sort;
            obj.wavelength_filternames = filternames_sorted;
            
            obj.isL = cellfun(@(x) any(strcmpi(x,obj.filternamesL)),filternames_sorted);
            obj.isR = cellfun(@(x) any(strcmpi(x,obj.filternamesR)),filternames_sorted);
            obj.eye_tags = cell(length(filterwavelength),1);
            [obj.eye_tags{obj.isL}] = deal('L'); [obj.eye_tags{obj.isR}] = deal('R');
        end
            
        function [] = set_img(obj)
            for i=1:length(obj.wavelength_filternames)
                if i==1
                    obj.img = obj.(obj.wavelength_filternames{i});
                else
                    obj.img = cat(3,obj.img,obj.(obj.wavelength_filternames{i}));
                end
            end
        end
        
        function [mstcamif] = rd2if_crism(obj,mastcam_images_SF,CDRSFlbl)
            d_km = CDRSFlbl.SOLAR_DISTANCE.value;
            [ d_au ] = km2au( d_km );
            mstcamif = MASTCAMImages();
            for i=1:length(obj.filternames)
                filname = obj.filternames{i};
                [mstcamif.(filname)] = rd2if_general(obj.(filname),mastcam_images_SF.(filname),d_au);
            end 
            mstcamif.set_wavelength();
            mstcamif.set_img();
        end
        
        function [spc,wv] = get_spectrum(obj,s,l,varargin)
            spc = squeeze(obj.img(l,s,:));
            wv = obj.wavelength;
        end
        
        function [spc,wv] = get_spectrumL(obj,s,l,varargin)
            spc = squeeze(obj.img(l,s,obj.isL));
            wv = obj.wavelength(obj.isL);
        end
        
        function [spc,wv] = get_spectrumR(obj,s,l,varargin)
            spc = squeeze(obj.img(l,s,obj.isR));
            wv = obj.wavelength(obj.isR);
        end
        
    end
end