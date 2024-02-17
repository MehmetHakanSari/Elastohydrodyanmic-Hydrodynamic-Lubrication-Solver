classdef parabolicslider
    %Parabolic slider with controlled curvature property
    properties
        alfa    %the curvature of the slider
        nodes    %mesh object 
        height  %return property
    end
    
    methods
         % Constructor method
        function obj = parabolicslider(varargin)
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

            obj = obj.set_initial_height();
        end


        function obj = set_initial_height(obj)
            x_non = linspace(0, 1, obj.nodes);
            obj.height = 4 * (1 - obj.alfa) * x_non.^2 + 4 * (obj.alfa - 1) * x_non + 1;
            obj.height =  (obj.height / max(obj.height));
        end


    end
end