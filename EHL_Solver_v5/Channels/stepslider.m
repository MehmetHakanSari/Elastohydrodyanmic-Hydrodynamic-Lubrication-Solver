classdef stepslider
    %Parabolic slider with controlled curvature property
    properties
        a   %the location of the step
        e   %the difference of the heights 
        mesh    %mesh object 
        height  %return property
        number_of_steps = 1 %default value of #_of_step
    end
    
    methods
         % Constructor method
        function obj = stepslider(varargin)
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
            x_non = linspace(0, 1, obj.mesh.N_x);
            h_1 = 1;
            h_0 = h_1 - obj.e;
            obj.height = (x_non <= obj.a) * h_1;
            obj.height = obj.height + (x_non > obj.a) * h_0;
            obj.height =  (obj.height / max(obj.height));
        end
    end
end