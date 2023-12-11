classdef pocket < channel
    
    properties
        type   %rear, center, front, multi
        a = 0.1;    %miniumum height of the channel
        c = 0.15;   
        b_front = 0.10;
        b_center = 0.425;
        b_rear = 0.75;
        height
        % h = h_min + (1-h_min) * tanh(100*(x-b)) / 2 - (1-h_min) * tanh(100*(x - b - c)) / 2;
    end
    
    methods
  
        function [obj, h, dhdx, d2hdx2] = create(obj,x)
          % return a parabolic slider height
          % x:     length of the channel. Usually row array. [1 x N]
            
            if obj.type == "front"
                h = obj.a + (1-obj.a) * tanh(100*(x-obj.b_front)) / 2 - (1-obj.a) * tanh(100*(x - obj.b_front - obj.c)) / 2;
            elseif obj.type == "center"
                h = obj.a + (1-obj.a) * tanh(100*(x-obj.b_center)) / 2 - (1-obj.a) * tanh(100*(x - obj.b_center - obj.c)) / 2;
            elseif obj.type == "rear"
                h = obj.a + (1-obj.a) * tanh(100*(x-obj.b_rear)) / 2 - (1-obj.a) * tanh(100*(x - obj.b_rear - obj.c)) / 2;
            elseif obj.type == "multi"
                %for non-diemnsional x
                for i = 1:length(x)
                if 0.05 < x(i) && x(i) <= 0.35
                    b = obj.b_front;
                    
                    h(i) = obj.a;
                    h(i) = h(i) + (1 - obj.a) * tanh(100 * (x(i) - b)) / 2;
                    h(i) = h(i) - (1 - obj.a) * tanh(100 * (x(i) - b - obj.c)) / 2;
                    
                elseif 0.37 < x(i) && x(i) <= 0.63
                    b = obj.b_center;
                    
                    h(i) = obj.a;
                    h(i) = h(i) + (1 - obj.a) * tanh(100 * (x(i) - b)) / 2;
                    h(i) = h(i) - (1 - obj.a) * tanh(100 * (x(i) - b - obj.c)) / 2;
                    
                elseif 0.64 < x(i) && x(i) <= 1
                    b = obj.b_rear;
                    
                    h(i) = obj.a;
                    h(i) = h(i) + (1 - obj.a) * tanh(100 * (x(i) - b)) / 2;
                    h(i) = h(i) - (1 - obj.a) * tanh(100 * (x(i) - b - obj.c)) / 2;
                    
                else
                    h(i) = obj.a;
                end
                end
            end
            
            dx = x(2) - x(1);               %If x is non-dimensional derivative is taken as non-dimensional.
            
            dhdx = OneDcentraldiff(h, dx);
            d2hdx2 = OneDcentraldiff(dhdx, dx);
            
            obj.height = h;
        end
          
    end
    
end