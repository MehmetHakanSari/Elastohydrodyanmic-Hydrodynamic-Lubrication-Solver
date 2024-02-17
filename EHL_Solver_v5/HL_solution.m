classdef HL_solution
    
    properties
        simulation_proporties
        mesh 
        fluid
        h 
        
        %LIN
        theta_LIN = containers.Map;
        pressure_LIN = containers.Map;
        pressure_Re

        velocity_field
        stress_LIN
        vel_logic_map
         
        %VR
        pressure_VR  
        pressure_field 
        stress_VR 
        theta_VR
        
        load_dic
        force_dic 
        friction_dic
        flowrate_dic
        
    end
    
    methods
        
        function obj = HL_solution(varargin)
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
            obj = obj.initilizer();
        end

        function [obj] = initilizer(obj)
            N_x = obj.mesh.N_x; 
            obj.pressure_LIN("p") = ones(1, N_x);
            obj.pressure_LIN("p0")  = ones(1, N_x);
            obj.pressure_LIN("p1") = zeros(1, N_x);
            obj.pressure_Re = ones(1, N_x);
            obj.theta_LIN("theta0") = ones(1, N_x);
            obj.theta_LIN("theta") = ones(1, N_x);
            obj.theta_LIN("theta1") = zeros(1, N_x);
            obj.h = obj.mesh.channel.height;
        end
        
        function [obj] = Reynolds(obj)
            [obj.pressure_Re] = ...
            ReynoldsCondiation(obj.h, obj.pressure_Re, obj.mesh, obj.fluid, obj.simulation_proporties);   
        end
        
        function [obj] = LIN(obj, varargin)
           
             if nargin > 0
                    for i = 1:2:numel(varargin)
                        if varargin{i} == "cavitation"
                            cavitation_flag = varargin{i + 1};
                        end
                    end
             end
            
            if cavitation_flag == "on"
                [obj.pressure_LIN, obj.theta_LIN, obj.flowrate_dic] = ...
                    LIN_EHL_cavitation(obj.h,  obj.pressure_LIN, obj.theta_LIN, obj.fluid, obj.mesh, obj.simulation_proporties);
            elseif cavitation_flag == "off"
                [obj.pressure_LIN] = ...
                    LIN_EHL(obj.h, obj.fluid, obj.mesh);
            end
            
        end
        
        function [obj] = calculate_Friction(obj)
            [obj.friction_dic] ...
                = friction_LIN(obj, obj.fluid, obj.mesh);
        end

        function [obj] = FieldProperties(obj)
            [obj.velocity_field, obj.stress_LIN, obj.force_dic] ...
                = StressVelocity_LIN(obj, obj.fluid, obj.mesh); %how to add optinal settings. Check it.
                obj.vel_logic_map = containers.Map;
                obj.vel_logic_map("u_uw") = obj.velocity_field("u") > 0;
                obj.vel_logic_map("v_uw") = obj.velocity_field("v") > 0;
                obj.vel_logic_map("u_dw") = obj.velocity_field("u") <= 0;
                obj.vel_logic_map("v_dw") = obj.velocity_field("v") <= 0;
        end
        
        function [obj] = VR(obj, dt, Spatial_Scheme, Temporal_Scheme, varargin)
            
             if nargin > 0
                    for i = 1:2:numel(varargin)
                        if varargin{i} == "cavitation"
                            cavitation_flag = varargin{i + 1};
                        end
                    end
             end
            
            [obj.pressure_VR, obj.theta_VR, obj.stress_VR] = ...
                VRCavitation(dt, Spatial_Scheme, Temporal_Scheme, obj, obj.mesh, obj.fluid, cavitation_flag);
        end
       
        function obj = Calculate_Load(obj)
                obj.load_dic = calculate_load(obj);
        end

    end
end