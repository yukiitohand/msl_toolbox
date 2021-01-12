classdef MSLGaleOrthoMosaic_v3 < HSI
    
    properties
        lblpath;
        lbl;
        proj_info
    end
     methods
        function obj = MSLGaleOrthoMosaic_v3(basename,dirpath,varargin)
            
            obj@HSI(basename,dirpath,varargin{:});
            obj.lblpath = joinPath(dirpath,[basename '_pds3.lbl']);
            obj.lbl = pds3lblread(obj.lblpath);
            obj.hdr = MSLGaleDEMMosaic_v3_lbl2hdr(obj.lbl);
            obj.get_proj_info();

        end
        
        function [img] = readimg(obj,varargin)
            [img] = msldem_lazyenvireadRect(obj,...
                0,0,obj.hdr.samples,obj.hdr.lines,'precision','uint8',varargin{:});
            if nargout<1
                obj.img = img;
            end
        end
        
        function get_proj_info(obj)
            [obj.proj_info] = mslGaleMosaic_get_cylindrical_proj_info(obj.lbl);
        end
        
        function [subimg] = get_subimage_wPixelRange(obj,xrange,yrange,varargin)
            precision = 'uint8';
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case 'PRECISION'
                            precision = varargin{i+1};
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            sample_offset = xrange(1)-1; line_offset = yrange(1)-1;
            samplesc = xrange(2)-xrange(1)+1; linesc = yrange(2)-yrange(1)+1;
            [subimg] = msldem_lazyenvireadRect(obj,...
                sample_offset,line_offset,samplesc,linesc,...
                'Precision',precision);
        end
        
        function [subimg,xrange,yrange] = get_subimage_wlatlon(obj,lon_range,lat_range,varargin)
            mrgn = 0;
            precision = 'uint8';
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case 'MARGIN'
                            mrgn = varargin{i+1};
                        case 'PRECISION'
                            precision = varargin{i+1};
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            [xrange] = obj.proj_info.get_xrange_wlon(lon_range,mrgn);
            [yrange] = obj.proj_info.get_yrange_wlat(lat_range,mrgn);
            s1 = max(1,xrange(1)); send = min(obj.hdr.samples,xrange(2));
            l1 = max(1,yrange(1)); lend = min(obj.hdr.lines,yrange(2));
            sample_offset = s1-1; line_offset = l1-1;
            samplesc = send-s1+1; linesc = lend-l1+1;
            [subimg] = msldem_lazyenvireadRect(obj,...
                sample_offset,line_offset,samplesc,linesc,...
                'Precision',precision);
            xrange = [s1 send];
            yrange = [l1 lend];
        end
        
%         function [subimg,xrange,yrange] = get_subimage_wlatlon(obj,lon_range,lat_range,varargin)
%             s1 = max(1,floor((lon_range(1)-obj.cylindrical_proj_info.lon_range(1))*obj.cylindrical_proj_info.rdlon));
%             send = min(obj.hdr.samples,ceil((lon_range(2)-obj.cylindrical_proj_info.lon_range(1))*obj.cylindrical_proj_info.rdlon)+1);
%             l1 = max(1,floor((obj.cylindrical_proj_info.lat_range(1)-lat_range(1))*obj.cylindrical_proj_info.rdlat));
%             lend = min(obj.hdr.lines,ceil((obj.cylindrical_proj_info.lat_range(1)-lat_range(2))*obj.cylindrical_proj_info.rdlat)+1);
%             sample_offset = s1-1;
%             line_offset = l1-1;
%             samplesc = send-s1+1;
%             linesc = lend-l1+1;
%             [subimg] = mslortho_lazyenvireadRect_mexw(obj,...
%                 sample_offset,line_offset,samplesc,linesc,varargin{:});
%             xrange = [s1 send];
%             yrange = [l1 lend];
%         end
     end
    
end