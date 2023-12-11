classdef domain
    
    properties
        R1 = 0.02;
        R2 = Inf;
        E1 = 10^7;
        E2 = Inf;
        nu1 = 0.5;
        nu2 = 0;
        mu = 0.1;
        Rx 
        E  
        
        N_x = 501;
        x
        dx
        bH
        pH
        K
        L
        
        h0_initial = 1e-6;
        h_max 
        height
        height_profile
        href
        epsilon
        domain_coeff
        
    end
    
    methods
        
        function [obj] = set_domain_coeff(obj, u)
           obj.domain_coeff = domain_coeff_function(u);
        end
          
        function [obj] = set_viscosity(obj, viscosity)
            obj.mu = viscosity;
        end
        
        function [obj] = set_eq_radii_curvature(obj)
            obj.Rx = 1 / (1 / obj.R1 + 1 / obj.R2); 
        end
        
        function [obj] = set_eq_elastic_modules(obj)
            obj.E = 1 / ( (1 - obj.nu1^2) / obj.E1 + (1 - obj.nu2^2) / obj.E2 ); 
        end
        
        function [obj] = set_bH(obj, F)
            obj.bH = (4 * F * obj.Rx / pi / obj.E) ^ 0.5;  
        end
        
        function [obj] = set_hmax(obj, h_max)
            obj.h_max = h_max;
         end
        
        function [obj] = set_ehl_length(obj, u)
            obj.L = 2 * obj.domain_coeff * obj.bH;
            obj.x = linspace(-obj.L / 2, obj.L / 2, obj.N_x);
            obj.dx = obj.x(2) - obj.x(1);
        end
        
        function [obj] = set_length(obj, given_length)
            obj.L = given_length;
            obj.x = linspace(0, obj.L, obj.N_x);
            obj.dx = obj.x(2) - obj.x(1);
        end
        
        function [obj] = set_hertizan_pressure(obj, force)
           obj.pH = 2 * force / (pi * obj.bH^2) .* ((obj.bH)^2 - obj.x.^2).^0.5; 
        end
        
        function obj = elastic_map(obj) 
            
            magnification = 10;                         
            fundemental_wavelength = magnification * obj.L;    
            npoints = magnification * obj.N_x;            % npoints for the L domain (calculated in order to have the same dx of the physical domain)

            [obj.K, Dlength,npuntiD ,GFcut] = ...
                ElasticKernel(obj.E1, obj.nu1, fundemental_wavelength, npoints, obj.dx, magnification);  
        end
        
        function obj = set_initial_height(obj, channel_type)
            if channel_type == "punch"
                obj.height = obj.x.^2 / (2 * obj.Rx) + obj.h0_initial;
                obj.height_profile = obj.x.^2 / (2 * obj.Rx);
            end
            
            if channel_type == "pocket"

            end
            
            if channel_type == "parabolicslider"
                alfa = 0.01;
                x_non = linspace(0, 1, obj.N_x);
                obj.height = 4 * (1 - alfa) * x_non.^2 + 4 * (alfa - 1) * x_non + 1;
                obj.height =  (obj.height / max(obj.height)) * obj.h_max;
            end
            
            if channel_type == "cosineslider"
                
            end

            if channel_type == "step"
                x_non = linspace(0, 1, obj.N_x);
                h_1 = 0.4;
                h_0 = 0.4 + 0.002;
                obj.height = (x_non <= 0.5) * h_1;
                obj.height = obj.height + (x_non > 0.5) * h_0;
            end

        end
        
        function obj = set_epsilon(obj, height_type, length_type)
            
            if isstring(height_type)
           
                if height_type == "max"           
                    h_ref = max(obj.height);
                elseif height_type == "center"
                    h_ref = obj.height(round(length(obj.height) / 2));
                elseif height_type == "min"
                    h_ref = min(obj.height);    
                elseif height_type == "curvature"
                    h_ref = obj.Rx;
                end               
            else
                h_ref = height_type; %double
            end
                
            obj.href = h_ref;
            
            if length_type == "max"
                L_ref = obj.L;
            elseif length_type == "center"
                 
            elseif length_type == "min"
                
            else    %predefined L
                L_ref = 2 * obj.bH;
            end
            obj.epsilon = h_ref / L_ref;
        end
        
        
    end
    
end


 