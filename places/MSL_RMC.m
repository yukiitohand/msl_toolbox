classdef MSL_RMC < handle
    %MSL rover motion counter class
    %  Properties
    %   SITE,DRIVE,POSE,RSM,ARM,CHIMRA,DRILL,HGA,DRT,IC
    %   -1 indicates that the value for the id is not set.
    
    properties
        SITE  = -1 % Site ID
        DRIVE = -1 % Drive ID
        POSE  = -1 % Pose ID
        RSM   = -1 % Remote Sensing Mast Motion Counter
        ARM   = -1
        CHIMRA = -1
        DRILL = -1
        HGA = -1
        DRT = -1
        IC = -1
    end
    
    methods
        function obj = MSL_RMC(varargin)
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        case 'SITE'
                            obj.SITE = varargin{i+1};
                        case 'DRIVE'
                            obj.DRIVE = varargin{i+1};
                        case 'POSE'
                            obj.POSE = varargin{i+1};
                        case 'RSM'
                            obj.RSM = varargin{i+1};
                        case 'ARM'
                            obj.ARM = varargin{i+1};
                        case 'CHIMRA'
                            obj.CHIMRA = varargin{i+1};
                        case 'DRILL'
                            obj.DRILL = varargin{i+1};
                        case 'HGA'
                            obj.HGA = varargin{i+1};
                        case 'DRT'
                            obj.DRT = varargin{i+1};
                        case 'IC'
                            obj.IC = varargin{i+1};
                        otherwise
                            error('Parameter: %s', varargin{i});   
                    end
                end
            end
        end
        function tf = eq(obj1,obj2)
            if (obj1.SITE == obj2.SITE) ...
                && (obj1.DRIVE == obj2.DRIVE) ...
                && (obj1.POSE == obj2.POSE) ...
                && (obj1.RSM == obj2.RSM) ...
                ...&& (obj1.ARM == obj2.ARM) ...
                ...&& (obj1.CHIMRA == obj2.CHIMRA) ...
                ...&& (obj1.DRILL == obj2.DRILL) ...
                ...&& (obj1.HGA == obj2.HGA) ...
                ...&& ((obj1.DRT == obj2.DRT) || (isnan(obj1.DRT) && isnan(obj2.DRT))) ...
                ...&& ((obj1.IC == obj2.IC) || (isnan(obj1.IC) && isnan(obj2.IC)))
                tf = true;
            else
                tf = false;
            end
        end
        function tf = ne(obj1,obj2)
            if (obj1.SITE ~= obj2.SITE) ...
                || (obj1.DRIVE ~= obj2.DRIVE) ...
                || (obj1.POSE ~= obj2.POSE) ...
                || (obj1.RSM ~= obj2.RSM) ...
                ...|| (obj1.ARM ~= obj2.ARM) ...
                ...|| (obj1.CHIMRA ~= obj2.CHIMRA) ...
                ...|| (obj1.DRILL ~= obj2.DRILL) ...
                ...|| (obj1.HGA ~= obj2.HGA) ...
                ...|| ((obj1.DRT ~= obj2.DRT) && ~(isnan(obj1.DRT) && isnan(obj2.DRT)) ) ...
                ...|| ((obj1.IC ~= obj2.IC) && ~(isnan(obj1.IC) && isnan(obj2.IC)))
                tf = true;
            else
                tf = false;
            end
        end

    end
end

