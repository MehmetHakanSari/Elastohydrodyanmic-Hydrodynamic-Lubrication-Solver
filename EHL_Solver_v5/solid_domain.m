classdef solid_domain
    
    properties
        R_top = 0.02;
        R2_bottom = Inf;
        E_top = 10^7;
        E_bottom = Inf;
        nu_top = 0.5;
        nu_bottom = 0;
        
        R 
        E  

        channel_type  %it is class that determines height of the channel 
        simulation    %it is the type of the simulation  
        mesh          %it is the mesh of the domain  
        
        bH
        p_hertzian
        deformation_kernel
        
        
        h0_initial = 1e-6;
        h_max 
        height
        height_profile
        href
        epsilon
        domain_coeff
        
    end
    
    methods
        function obj = solution(varargin)
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
      
        function [obj] = set_domain_coeff(obj, speed)
           obj.domain_coeff = domain_coeff_function(speed);
        end
          
        function [obj] = set_eq_radii_curvature(obj)
            obj.Rx = 1 / (1 / obj.R_top + 1 / obj.R2_bottom); 
        end
        
        function [obj] = set_eq_elastic_modules(obj)
            obj.E = 1 / ( (1 - obj.nu_top^2) / obj.E_top + (1 - obj.nu_bottom^2) / obj.E_bottom); 
        end
        
        function [obj] = set_bH(obj, force)
            obj.bH = (4 * force * obj.R / pi / obj.E) ^ 0.5;  
        end
        
        function [obj] = set_ehl_length(obj)
            obj.L = 2 * obj.domain_coeff * obj.bH;
            obj.x = linspace(-obj.L / 2, obj.L / 2, obj.N_x);
            obj.dx = obj.x(2) - obj.x(1);
        end
        
        function [obj] = set_hertizan_pressure(obj, force)
           obj.p_hertzian = 2 * force / (pi * obj.bH^2) .* ((obj.bH)^2 - obj.mesh.dimensional("x").^2).^0.5; 
        end
        
        function obj = elastic_map(obj) 
            
            magnification = 10;                         
            fundemental_wavelength = magnification * obj.mesh.length;    
            N_magnified = magnification * obj.mesh.N_x;            % npoints for the L domain (calculated in order to have the same dx of the physical domain)

            [obj.deformation_kernel, Dlength, npuntiD ,GFcut] = ...
                ElasticKernel(obj.E1, obj.nu1, fundemental_wavelength, N_magnified, obj.mesh.dimensional("dx"), magnification);  
        end
        
        function obj = set_channel(obj, channel)
            obj.channel = channel
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


 