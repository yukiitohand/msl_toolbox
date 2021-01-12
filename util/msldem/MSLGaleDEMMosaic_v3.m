classdef MSLGaleDEMMosaic_v3 < HSI
    
    properties
        lblpath;
        lbl;
        proj_info
    end
     methods
        function obj = MSLGaleDEMMosaic_v3(basename,dirpath,varargin)
            
            obj@HSI(basename,dirpath,varargin{:});
            obj.lblpath = joinPath(dirpath,[basename '.lbl']);
            obj.lbl = pds3lblread(obj.lblpath);
            obj.hdr = MSLGaleDEMMosaic_v3_lbl2hdr(obj.lbl);
            obj.get_proj_info();

        end
        
        function get_proj_info(obj)
            [obj.proj_info] = mslGaleMosaic_get_cylindrical_proj_info(obj.lbl);
        end
        
        function [subimg] = get_subimage_wPixelRange(obj,xrange,yrange,varargin)
            precision = 'double';
            mex_mode = 1;
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case 'PRECISION'
                            precision = varargin{i+1};
                        case 'MEX_MODE'
                            mex_mode = varargin{i+1};
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            sample_offset = xrange(1)-1; line_offset = yrange(1)-1;
            samplesc = xrange(2)-xrange(1)+1; linesc = yrange(2)-yrange(1)+1;
            [subimg] = msldem_lazyenvireadRect(obj,...
                sample_offset,line_offset,samplesc,linesc,...
                'Precision',precision,'MEX_MODE',mex_mode);
        end
        
        function [subimg,xrange,yrange] = get_subimage_wlatlon(obj,lon_range,lat_range,varargin)
            mrgn = 0;
            precision = 'double';
            mex_mode = 1;
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case 'MARGIN'
                            mrgn = varargin{i+1};
                        case 'PRECISION'
                            precision = varargin{i+1};
                        case 'MEX_MODE'
                            mex_mode = varargin{i+1};
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            [xrange] = obj.proj_info.get_xrange_wlon(lon_range,mrgn);
            [yrange] = obj.proj_info.get_yrange_wlat(lat_range,mrgn);
            s1 = max(1,xrange(1)); send = min(obj.hdr.samples,xrange(2));
            l1 = max(1,yrange(1)); lend = min(obj.hdr.lines,yrange(2));
%             s1 = max(1,...
%                 ceil( (lon_range(1)-obj.proj_info.lon_range(1))*obj.proj_info.rdlon ) - mrgn ...
%                 );
%             send = min(obj.hdr.samples,...
%                 ceil( (lon_range(2)-obj.cylindrical_proj_info.lon_range(1))*obj.cylindrical_proj_info.rdlon ) + mrgn ...
%                 );
%             l1 = max(1,...
%                 ceil( (obj.cylindrical_proj_info.lat_range(1)-lat_range(1))*obj.cylindrical_proj_info.rdlat ) - mrgn ...
%                 );
%             lend = min(obj.hdr.lines,...
%                 ceil( (obj.cylindrical_proj_info.lat_range(1)-lat_range(2))*obj.cylindrical_proj_info.rdlat ) + mrgn ...
%                 );
            sample_offset = s1-1; line_offset = l1-1;
            samplesc = send-s1+1; linesc = lend-l1+1;
            [subimg] = msldem_lazyenvireadRect(obj,...
                sample_offset,line_offset,samplesc,linesc,...
                'Precision',precision,'MEX_MODE',mex_mode);
            xrange = [s1 send];
            yrange = [l1 lend];
        end
        
     end
    
end