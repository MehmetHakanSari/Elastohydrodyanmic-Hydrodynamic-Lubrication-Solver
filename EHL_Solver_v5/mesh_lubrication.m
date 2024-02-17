classdef mesh_lubrication
    
    properties
       
        channel  %it is class that determines height of the channel 
        %non-dimensional quantaties
        N_x
        N_y
        x
        dx
        X
        dy
        Y

        %dimensional quantaties
        h_max       %maximum height of interest
        L_max       %maximum length of interest
        
        %dimensional
        dimensional = containers.Map;
        h0_initial 
    
        href
        epsilon
        domain_coeff
        Jacobian
        
    end
    
    methods
        
        function obj = mesh_lubrication(varargin)
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

        function obj = ehl_mesh(obj, domain)
                %sets dimensional ehl mesh properties from given domain
                function [obj] = set_ehl_length(obj, u)
                    obj.L = 2 * obj.domain.ehl{domain_coeff} * obj.domain.ehl{bH}; %I ll translate them to hashmap
                    obj.x = linspace(-obj.L / 2, obj.L / 2, obj.N_x);
                    obj.dx = obj.x(2) - obj.x(1);
                end
                
          end 
        
        function obj = set_channel(obj, channel)
            obj.channel = channel;
        end

        function obj = consturct_mesh_matrix(obj)
            %GRIDING
            obj.X = zeros(obj.N_y, obj.N_x);
            obj.Y = zeros(obj.N_y, obj.N_x);
            obj.x = linspace(0, 1, obj.N_x);
            obj.dx = 1 / (obj.N_x - 1);
            for i = 1 : obj.N_y
                obj.X(i,:) = obj.x;
            end
            obj.dy = obj.channel.height / (obj.N_y - 1);
            for i = 1 : obj.N_y
                obj.Y(i,:) = obj.channel.height - obj.dy .* (i-1);
            end
        end

        function obj = set_Jacaboian(obj)
            obj.Jacobian = containers.Map; %hashmap

            dhdx = OneDcentraldiff(obj.channel.height,  obj.dx, "CD2");
            % dydx = -Y ./ h.^2 .* dhdx;          %Coordiante transform y' = y/h where the rest of the code uses y' as y.
            dydx = -obj.Y ./ obj.channel.height .* dhdx;             %Modified by multiplying "h". where it come from [d(.)/dy' = h*d(.)/dy]
            obj.Jacobian("dydx") = dydx;
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


 