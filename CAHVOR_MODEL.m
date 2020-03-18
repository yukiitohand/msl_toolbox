classdef CAHVOR_MODEL < handle
    % CAHVOR CAMERA model class
    %  CAHVOR model information stored. You can only have CAHV. 
    %  Properties
    %   C, A, H, V, O, R: model components. 1x3 vector.
    %   type: type of the model {'CAHV','CAHVOR'}
    
    properties
        type   % 'CAHVOR' or 'CAHV'
        C = []
        A = []
        H = []
        V = []
        O = []
        R = []
    end
    
    methods
        function obj = CAHVOR_MODEL(varargin)
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
                        case {'O','OPTICAL'}
                            if ~isempty(varargin{i+1})
                                obj.O = reshape(varargin{i+1},1,[]);
                            end
                        case {'R','RADIAL'}
                            if ~isempty(varargin{i+1})
                                obj.R = reshape(varargin{i+1},1,[]);
                            end
                        otherwise
                            error('Parameter: %s', varargin{i});   
                    end
                end
            end
            if isempty(obj.O) && isempty(obj.R)
                obj.type = 'CAHV';
            elseif ~isempty(obj.O) && ~isempty(obj.R)
                obj.type = 'CAHVOR';
            else
                error('Only either of O and R is given');
            end
        end
        function tf = eq(obj1,obj2)
            if all(obj1.C == obj2.C) ...
                && all(obj1.A == obj2.A) ...
                && all(obj1.H == obj2.H) ...
                && all(obj1.V == obj2.V)
                ...&& all(obj1.O == obj2.O) ...
                ...&& all(obj1.R == obj2.R)
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
                ...|| (all(obj1.O ~= obj2.O)) ...
                ...|| (all(obj1.R ~= obj2.R))
                tf = true;
            else
                tf = false;
            end
        end

    end
end