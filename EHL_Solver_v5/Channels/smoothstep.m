classdef smoothstep
    %Parabolic slider with controlled curvature property
    properties
        a   %the location of the step
        e   %the difference of the heights 
        c  % curvature of the step
        nodes    %the node number of height
        height  %return property
        number_of_steps = 1 %default value of #_of_step
    end
    
    methods
         % Constructor method
        function obj = smoothstep(varargin)
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
            h_1 = 1;
            k = h_1 -  abs(obj.e) / 2;
            obj.height = k - sign(obj.e) * (1 - k) * tanh(obj.c *(x_non - obj.a)) ;
            obj.height =  (obj.height / max(obj.height));            
        end
    end
end

         
    
