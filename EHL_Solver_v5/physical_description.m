classdef physical_description
    properties
        U
        L
        eta0
        lambda
        h0
        p0
    end

    methods
        function obj = physical_description(varargin)
            if nargin > 0
                for i = 1:2:numel(varargin)
                    propertyName = varargin{i};
                    if isprop(obj, propertyName)
                        obj.(propertyName) = varargin{i + 1};
                    else
                        error('Invalid property name: %s', propertyName);
                    end
                end
            end
        end
    end
end