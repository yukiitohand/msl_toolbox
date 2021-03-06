classdef CAHV_MODEL < handle
    % CAHV CAMERA model class
    %  CAHVOR model information stored. You can only have CAHV. 
    %  Properties
    %   C, A, H, V: model components. 1x3 vector.
    
    properties
        type   %'CAHV'
        C
        A
        H
        V
        hc
        vc
        hs
        vs
        Hdash
        Vdash
        REFERENCE_COORD_SYSTEM
    end
    
    methods
        function obj = CAHV_MODEL(varargin)
            obj.C = []; obj.A = []; obj.H = []; obj.V = [];
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case {'C','CENTER'}
                            obj.C = reshape(varargin{i+1},1,[]);
                        case {'A','AXIS'}
                            obj.A = reshape(varargin{i+1},1,[]);
                        case {'H','HORIZONTAL'}
                            obj.H = reshape(varargin{i+1},1,[]);
                        case {'V','VERTICAL'}
                            obj.V = reshape(varargin{i+1},1,[]);
                        case {'REFERENCE_COORD_SYSTEM'}
                            obj.REFERENCE_COORD_SYSTEM = varargin{i+1};
                        otherwise
                            error('Parameter: %s', varargin{i});   
                    end
                end
            end
            obj.type = 'CAHV';
        end
        function tf = eq(obj1,obj2)
            if all(obj1.C == obj2.C) ...
                && all(obj1.A == obj2.A) ...
                && all(obj1.H == obj2.H) ...
                && all(obj1.V == obj2.V)
                tf = true;
            else
                tf = false;
            end
        end
        function tf = ne(obj1,obj2)
            if (all(obj1.C ~= obj2.C)) ...
                || (all(obj1.A ~= obj2.A)) ...
                || (all(obj1.H ~= obj2.H)) ...
                || (all(obj1.V ~= obj2.V))
                tf = true;
            else
                tf = false;
            end
        end
        
        function [hc,vc] = get_optical_center_pixel(obj)
            hc = (obj.A*obj.H'); %./norm(obj.A); (A should be normalized)
            vc = (obj.A*obj.V'); %./norm(obj.A); (A should be normalized)
            obj.hc = hc; obj.vc = vc;
        end
        
        function [hs,vs] = get_pixel_focal_length(obj)
            hs = norm(cross(obj.A,obj.H'));
            vs = norm(cross(obj.A,obj.V'));
            obj.hs = hs; obj.vs = vs;
        end
        
        function [Hdash,Vdash,hs,vs,hc,vc] = get_image_plane_unit_vectors(obj)
            [hc,vc] = obj.get_optical_center_pixel();
            [hs,vs] = obj.get_pixel_focal_length();
            Hdash = (obj.H-hc*obj.A) / hs;
            Vdash = (obj.V-vc*obj.A) / vs;
            obj.Hdash = Hdash;
            obj.Vdash = Vdash;
        end
        
        function [p_minus_c] = get_p_minus_c_from_xy(obj,xy)
            [p_minus_c] = cahv_get_p_minus_c_from_xy_v2(xy,obj);
        end
        
        function [xy] = get_xy_from_p_minus_c(obj,pmc)
            [xy] = cahv_get_xy_from_p_minus_c(pmc,obj);
        end
        
        function [pmc_h,pmc_v] = get_pmc_FOV(obj,h_range,v_range)
            obj.get_image_plane_unit_vectors();
            pmc_lrtb   = obj.get_p_minus_c_from_xy(...
                [ h_range(1) obj.vc    ;
                  h_range(2) obj.vc    ;
                  obj.hc     v_range(1);
                  obj.hc     v_range(2) ]');
            pmc_h = pmc_lrtb(:,1:2);
            pmc_v = pmc_lrtb(:,3:4);
            % pmc_l_vc   = obj.get_p_minus_c_from_xy([h_range(1) obj.vc]);
            % pmc_r_vc   = obj.get_p_minus_c_from_xy([h_range(2) obj.vc]);
            % pmc_top_hc = obj.get_p_minus_c_from_xy([obj.hc v_range(1)]);
            % pmc_btm_hc = obj.get_p_minus_c_from_xy([obj.hc v_range(2)]);
            % pmc_h = [pmc_l_vc' pmc_r_vc'];
            % pmc_v = [pmc_top_hc' pmc_btm_hc'];
        end
        
        function [fovh,fovv] = get_FOV(obj,h_range,v_range)
            [pmc_h,pmc_v] = obj.get_pmc_FOV(h_range,v_range);
            fovl = acos(obj.A*pmc_h(:,1));
            fovr = acos(obj.A*pmc_h(:,2));
            fovt = acos(obj.A*pmc_v(:,1));
            fovb = acos(obj.A*pmc_v(:,2));
           
            fovh = [fovl,fovr];
            fovv = [fovt;fovb];
        end

        function [pmc_ul,pmc_ll,pmc_lr,pmc_ur] = get_pmc_FOVvertex(obj,h_range,v_range)
            obj.get_image_plane_unit_vectors();
            pmc_ul_counter_clock   = obj.get_p_minus_c_from_xy(...
                [ h_range(1) v_range(1);
                  h_range(1) v_range(2);
                  h_range(2) v_range(1);
                  h_range(2) v_range(2) ]');
            pmc_ul = pmc_ul_counter_clock(:,1);
            pmc_ll = pmc_ul_counter_clock(:,2);
            pmc_lr = pmc_ul_counter_clock(:,3);
            pmc_ur = pmc_ul_counter_clock(:,4);
        end
        
        % function [fovh,fovv] = get_FOV(obj,h_range,v_range)
        %     obj.get_optical_center_pixel();
        %     pmc_l_vc = cahv_get_p_minus_c_from_xy([h_range(1) obj.vc],obj);
        %     pmc_r_vc = cahv_get_p_minus_c_from_xy([h_range(2) obj.vc],obj);
        %     pmc_top_hc = cahv_get_p_minus_c_from_xy([obj.hc v_range(1)],obj);
        %     pmc_btm_hc = cahv_get_p_minus_c_from_xy([obj.hc v_range(2)],obj);
        %     fovl = acos(pmc_l_vc*obj.A');
        %     fovr = acos(pmc_r_vc*obj.A');
        %     fovt = acos(pmc_top_hc*obj.A');
        %     fovb = acos(pmc_btm_hc*obj.A');
        %    
        %     fovh = [fovl,fovr];
        %     fovv = [fovt;fovb];
        %end

    end
end