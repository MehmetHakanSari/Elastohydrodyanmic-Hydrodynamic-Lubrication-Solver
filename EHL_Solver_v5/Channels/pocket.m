classdef pocket 
    
    properties
        type   %rear, center, front, multi
        a = 0.1;    %miniumum height of the channel
        c = 0.15;   %the length of the pocket
        e   %the difference of the heights 
        b_front = 0.10;
        b_center = 0.425;
        b_rear = 0.75;
        curv  %curvature
        height
        nodes
        
    end
    
    methods
  
          function obj = pocket(varargin)
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
        
        function [obj] = set_initial_height(obj)
          % return a parabolic slider height
          % x:     length of the channel. Usually row array. [1 x N]
            x = linspace(0, 1, obj.nodes);
            h_1 = 1;
            k = h_1 -  abs(obj.e) / 2;
%             plot(k - sign(obj.e) * (1 - k) * tanh(obj.curv*(x - obj.b_center)) + sign(obj.e) * (1 - k) * tanh(obj.curv *(x - obj.b_center - obj.c)) + obj.e/2 );
            if obj.type == "front"
                h = obj.a + (1 - obj.a) * tanh(obj.curv * (x - obj.b_front)) / 2 - (1-obj.a) * tanh(obj.curv*(x - obj.b_front - obj.c)) / 2;
            elseif obj.type == "center"
                obj.b_center = (1 - obj.c) / 2 ; 
                h = k - sign(obj.e) * (1 - k) * tanh(obj.curv*(x - obj.b_center)) + sign(obj.e) * (1 - k) * tanh(obj.curv *(x - obj.b_center - obj.c)) + obj.e/2 ;
%                 h = h/max(h);
%                 h = obj.a + (1-obj.a) * tanh(obj.curv*(x-obj.b_center)) / 2 - (1-obj.a) * tanh(obj.curv*(x - obj.b_center - obj.c)) / 2;
            elseif obj.type == "rear"
                h = obj.a + (1-obj.a) * tanh(obj.curv*(x-obj.b_rear)) / 2 - (1-obj.a) * tanh(obj.curv*(x - obj.b_rear - obj.c)) / 2;
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
                        
            obj.height = h;
        end
          
    end
    
end