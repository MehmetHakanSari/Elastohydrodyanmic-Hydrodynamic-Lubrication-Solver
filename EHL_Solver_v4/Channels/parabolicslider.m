classdef parabolicslider < channel
    
    properties
        alfa = 0.5
        height
    end
    
    methods
  
        function [obj, h, dhdx, d2hdx2] = create(obj,x)
          % return a parabolic slider height
          % x:     length of the channel. Usually row array. [1 x N]
          % alpha: curvature parameter.   Double
            h = 4 * (1 - obj.alfa) * x.^2 + 4 * (obj.alfa - 1) * x + 1;
            dx = x(2) - x(1);               %If x is non-dimensional derivative is taken as non-dimensional.
            dhdx = OneDcentraldiff(h, dx);
            d2hdx2 = OneDcentraldiff(dhdx, dx);
            obj.height = h;
        end
          
    end
    
end