classdef channel 
    
    properties     %physical properties
       h_max       %maxiumum channel height  
       h_min       %minimum channel height
       L           %length of the channel
       U           %velocity of upper wall
       mu          %viscosity of the lubricant
    end
    
    methods
        
        function [epsilon] = epsilon_max(obj)
           epsilon = obj.h_max / obj.L; 
        end
        
        function [epsilon] = epsilon_min(obj, len)
            %len: user defined length scale for minimum epsilon
           epsilon = obj.h_min / len; 
        end
        
        function obj = set_properties(obj, h_max, L, U, mu)
           obj.h_max = h_max;
           obj.L = L;
           obj.U = U;
           obj.mu = mu;
           
           if length(obj.height) > 2
               %calculate minimum h 
           end
        end
        
        function plot_nd(obj)
           %plot non-dimensional channel height 
           %use with subclass (channel type)
           %I abonden it yani. Stupid matlab. Didnt understand how to
           %control axes inside a class. I did not understand axes, gca at
           %all. 
           
           N = length(obj.height);
           
           figure(1)
           plot(linspace(0,1,N), obj.height,"k","LineWidth",2);
          
%            plot(linspace(0,1,N), zeros(1,N),"k","LineWidth",1.5);
%            xlabel("\bf{x}","Interpreter","latex")
%            ylabel("\bf{h}","Interpreter","latex")
%            set(gca,"Fontsize",20);
%            set(gca,'TickLabelInterpreter','latex')
%            box on;
%            grid off;
%            set(gca, 'Position',[0 200 800 550])
%            
           
           
        end
        
    end
    
    
    
    
    
    
    
end