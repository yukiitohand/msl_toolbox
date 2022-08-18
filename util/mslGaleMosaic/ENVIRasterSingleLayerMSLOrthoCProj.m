classdef ENVIRasterSingleLayerMSLOrthoCProj < ENVIRasterSingleLayerEquirectProjRot0
    % ENVIRasterSingleLayerMSLOrthoCProj
    %  ENVI style (*.img & *.hdr) single-layered image defined on a 
    %  rectangular region cropped from MSLOrthoGaleMosaic_v3 data. 
    %  "C" means "Cropped".
    % 
    %  chdr: struct storing information on the cropped region, composed of
    %    of following fields
    %      sample_offset
    %      line_offset
    %      samples
    %      lines
    %  base: point to the obj of MSLDEMGaleMosaic_v3, base DEM data.
    properties
        chdr; % stroing crop information
        base;
    end
     methods
        function obj = ENVIRasterSingleLayerMSLOrthoCProj(basename,dirpath,varargin)
            
            obj@ENVIRasterSingleLayerEquirectProjRot0(...
                basename,dirpath,varargin{:});
            
            if ~isempty(obj.hdr)
                obj.get_proj_info();

                obj.chdr.samples = obj.hdr.samples;
                obj.chdr.lines   = obj.hdr.lines;

                if isfield(obj.hdr,'mslorthoc_sample_offset')
                    obj.chdr.sample_offset = obj.hdr.mslorthoc_sample_offset;
                else
                    error('mslorthoc_sample_offset is not defined in the header file.');
                end
                if isfield(obj.hdr,'mslorthoc_line_offset')
                    obj.chdr.line_offset   = obj.hdr.mslorthoc_line_offset;
                else
                    error('mslorthoc_line_offset is not defined in the header file.');
                end
                
                obj.base = MSLGaleOrthoMosaic_v3(obj.hdr.mslortho_basename,obj.hdr.mslortho_dirpath);
                % obj.get_proj_info_from_base();
                
            end
        end
        
        function [] = get_proj_info_from_base(obj)
            proj_info_from_base = SphereEquiRectangularProj( ...
                'Radius',obj.base.proj_info.radius                           , ...
                'STANDARD_PARALLEL',obj.base.proj_info.standard_parallel     , ...
                'CenterLongitude',obj.base.proj_info.center_longitude        , ...
                'Latitude_of_origin',obj.base.proj_info.latitude_of_origin   , ...
                'Longitude_of_origin',obj.base.proj_info.longitude_of_origin   ...
                );
            proj_info_from_base.rdlat = obj.base.proj_info.rdlat;
            proj_info_from_base.rdlon = obj.base.proj_info.rdlon;
            proj_info_from_base.map_scale_x = obj.base.proj_info.map_scale_x;
            proj_info_from_base.map_scale_y = obj.base.proj_info.map_scale_y;
            proj_info_from_base.set_lon1(obj.base.longitude(obj.get_x_base(1)));
            proj_info_from_base.set_lat1(obj.base.latitude(obj.get_y_base(1)));
            proj_info_from_base.longitude_range = obj.base.proj_info.get_lon_range(obj.get_xrange_base());
            proj_info_from_base.latitude_range  = obj.base.proj_info.get_lat_range(obj.get_yrange_base());
            obj.proj_info = proj_info_from_base;
        end
        
        function [x_base] = get_x_base(obj,x)
            x_base = x + obj.chdr.sample_offset;
        end
        
        function [y_base] = get_y_base(obj,y)
            y_base = y + obj.chdr.line_offset;
        end
        
        function [xy_base] = get_xy_base(obj,xy)
            % xy needs to be [* x 2] size. Each row is considered to be
            % [x,y]. If only one [x,y] input, the vector can be either of
            % column or row.
            if iscol(xy)
                xy = xy';
                xy_iscol = true;
            else
                xy_iscol = false;
            end
            xy_base = xy + [obj.chdr.sample_offset obj.chdr.line_offset];
            if xy_iscol
                xy = xy';
            end
            
        end
        
        function [xrange_base] = get_xrange_base(obj)
            xrange_base = obj.get_x_base([1 obj.hdr.samples]);
        end
        
        function [yrange_base] = get_yrange_base(obj)
            yrange_base = obj.get_y_base([1 obj.hdr.lines]);
        end
        
        function [xrange_base,yrange_base] = get_xyrange_base(obj)
            xrange_base = obj.get_x_base([1 obj.hdr.samples]);
            yrange_base = obj.get_y_base([1 obj.hdr.lines]);
        end
        
        
     end
end