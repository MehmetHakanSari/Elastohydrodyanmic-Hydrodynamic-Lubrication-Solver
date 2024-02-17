classdef EHL_solution
    
    properties
        simulation_proporties
        mesh 
        fluid
        solid_domain 

        pressure_Re
        
        ehl_properties  %h pressure theta 
        
        force 

        %LIN
        LIN_data

        %VR
        VR_data
        
        h0
        h
        displacement %dimensional displacement
        r            %relexation_parameter

        load_dic

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

            obj = obj.initilizer();
        end

        function [obj] = initilizer(obj)
            N_x = obj.mesh.N_x;   
            obj.pressure = ones(1, N_x);
            obj.pressure_Re = ones(1, N_x);
            obj.displacement = zeros(1, N_x);
            obj.r = zeros(1, N_x);
            obj.theta0 = ones(1, N_x);
            obj.theta = ones(1, N_x);
            obj.theta1 = zeros(1, N_x);
            obj.pressure0 = ones(1, N_x);
            obj.pressure1= zeros(1, N_x);
            obj.h = obj.mesh.channel.height;   
        end
       
        function [obj, loadflag] = loadsolver(obj, F) 
            %balance the load for a punch for unfixed film thickness
            %F: force (N)
            [obj, loadflag] = LoadSolver(obj, F);                  
        end
        
        function [obj, dspflag] = displacementsolver(obj)
            [obj, dspflag] = DisplacementSolver(obj);
        end
           
        function [obj] = Reynolds(obj)
            [obj.pressure_Re] = ..
            ReynoldsCondiation(obj.h, obj.pressure_Re, obj.mesh, obj.fluid, obj.simulation_proporties)   
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

                [obj.pressure, obj.pressure0, obj.pressure1, obj.thetaT, obj.thetaD, obj.theta, obj.flowrate_dic] = ...
                    LIN_EHL_cavitation(obj.h, obj.pressure, obj.pressure0, obj.pressure1, obj.thetaT, obj.thetaD, obj.theta, obj.fluid_domain, obj.mesh);

            elseif cavitation_flag == "off"

                [obj.pressure, obj.pressure0, obj.pressure1] = ...
                    LIN_EHL(obj.h, obj.fluid, obj.mesh);
            end
        end
        
        function [obj] = calculate_Friction(obj)
            [obj.friction] ...
                = friction_LIN(obj, obj.fluid, obj.mesh);
        end

        function [obj] = FieldProperties(obj)
            [obj.velocity_field, obj.stress_LIN, obj.force] ...
                = StressVelocity_LIN(obj, obj.fluid, obj.mesh); %how to add optinal settings. Check it.
            obj = obj.velocity_logic_map();

            function obj = velocity_logic_map(obj)
                obj.vel_logic_map = containers.Map;
                obj.vel_logic_map("u_uw") = obj.velocity_field("u") > 0;
                obj.vel_logic_map("v_uw") = obj.velocity_field("v") > 0;
                obj.vel_logic_map("u_dw") = obj.velocity_field("u") <= 0;
                obj.vel_logic_map("v_dw") = obj.velocity_field("v") <= 0;
            end 
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
            if isempty(obj.pressure_VR)
                W_N= trapz(obj.mesh.x, obj.pressure_Newtonain);
                W_NN = trapz(obj.mesh.x, obj.pressure);
            else
                W_N= trapz(obj.mesh.x, obj.pressure_Newtonain);
                W_NN = trapz(obj.mesh.x, obj.pressure);
                W_VR = trapz(obj.mesh.x, obj.pressure_VR("p"));
            end
            obj.load_dic = containers.Map({char('W_N'), char('W_NN'), char('W_VR')}, {W_N, W_NN, W_VR});
        end
        
    end
    
end