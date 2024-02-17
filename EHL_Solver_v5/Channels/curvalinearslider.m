classdef curvalinearslider
    %Parabolic slider with controlled curvature property
    properties
        Rx   %radii of curvature
        R2   %the difference of the heights 
        mesh    %mesh object 
        height  %return property
        h0_initial %initial height
    end
    
    methods
         % Constructor method
        function obj = curvalinearslider(varargin)
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

            obj = set_initial_height(obj)
        end

        function obj = set_initial_height(obj)
            obj.height = obh.mesh.dimensional("x").^2 / (2 * obj.Rx) 
            obj.height = obj.height + obj.h0_initial
        end    
    end
end