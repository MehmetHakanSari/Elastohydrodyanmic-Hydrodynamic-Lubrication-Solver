classdef solution
    
    properties
        velocity
        applied_load
        viscocity_ratio
        wiessenberg_Number 
        deborah_Number
        hersey_Number 
        
        p_cavitation = 0
        p_out = 0
        p_in = 0
        
        relaxp = 1.95; % Sovrarelaxation rate on Reynolds pressure
        relaxD = 0.1; % Underrelaxation rate on EHL dispalacement (initial value)
        relaxL = 0.5; % Underrelaxation rate on load control

        relaxDmax = 1e-2; % maximum underrelaxation rate on EHL dispalacement for Aitk. acc
        relaxDmin = 1e-3; % minimum underrelaxation rate on EHL dispalacement for Aitk. acc

        tollp = 1e-9; % Error on on Reynolds pressure
        tollD = 1e-6; % Error on EHL dispalacement
        tollL = 1e-3; % Error on EHL load control
        
        errD
        errL
        
        h0 = 1e-6;
        domain
        hertzian_pressure 
        
        friction
        stress_field
        velocity_field
        pressure_field
        force 
         
%         p_old 
%         dsp_old
%         r_old 
        
        theta_old 
        thetaT_old 
        thetaD_old 
        pressure0_old 
        pressure1_old 
            
        h 
        pressure
        theta
        displacement 
        r
        
        pressure_VR  
        stress_VR 
        theta_VR
       
    end
    
    methods
        
        function [obj] = calculate_He(obj)
            obj.hersey_Number = obj.domain.mu * obj.domain.L * ...
                obj.velocity / obj.applied_load;
        end

        function [obj] = set_solution_parameters(obj, Wi, beta, u, F)
            obj.velocity = u;
            obj.viscocity_ratio = beta;
            obj.applied_load = F;
            obj.wiessenberg_Number = Wi;
            obj.deborah_Number = Wi * obj.domain.epsilon;
            
            if Wi < 0.5 % Deborah Number Approximation
                obj.wiessenberg_Number =  Wi / obj.domain.epsilon; 
                obj.deborah_Number =  Wi; 
            end

        end
        
        function [obj] = initilizer(obj)
            
            N_x = obj.domain.N_x; 
            
            obj.stress_field = {zeros(N_x, N_x), zeros(N_x, N_x), zeros(N_x, N_x)};
            
            obj.pressure = ones(1, N_x);
            obj.displacement = zeros(1, N_x);
            obj.r = zeros(1, N_x);

            obj.theta_old = ones(1, N_x);
            obj.thetaT_old = ones(1, N_x);
            obj.thetaD_old = zeros(1, N_x);
            obj.pressure0_old = ones(1, N_x);
            obj.pressure1_old = zeros(1, N_x);
     
            obj.h = obj.domain.height;
            
            if isempty(obj.deborah_Number)
                obj.deborah_Number = obj.wiessenberg_Number * obj.domain.epsilon;
            end
            if isempty(obj.wiessenberg_Number)
                obj.wiessenberg_Number = obj.deborah_Number / obj.domain.epsilon;
            end
            
        end
        
        function [obj, loadflag] = loadsolver(obj, F) 
            %F: given force
            %pressure: pressure output taken from the finished displacement
            %solver.
            
            domain = obj.domain;
            
            F_new = real(trapz(domain.x, obj.pressure));
            
            % error on Load
            obj.errL = abs(F_new - F) / abs(F);            
            
            if obj.errL < obj.tollL
                loadflag = false;
            else % Load Control on gentral gap
                
                h00 = obj.h0 - obj.displacement(floor(domain.N_x / 2 + 1));  %displacement is column vector
                h00n = h00 * (F_new / F) ^ obj.relaxL;
                h0n = h00n + obj.displacement(floor(domain.N_x / 2 + 1));
                
                obj.h0 = h0n;
                obj.h = obj.h0 + domain.height_profile - obj.displacement;  
                loadflag = true;
            end
                  
        end
        
        function [obj, dspflag] = displacementsolver(obj)
            
            domain = obj.domain;
            
            Dsp_old = obj.displacement;
                                
            % Elastic displacement
            Dsp = (domain.K * (obj.pressure)')';            % (Nx1) = (NxN) (Nx1) changing Dps from (Nx1) to (1xN)
            
            % Aitken acceleration
            r_old = (obj.r)';                           %transpoze due to store preferences. 
            
            r = ( abs( (Dsp - Dsp_old) ./ Dsp) )';
            
            
            relaxD = obj.relaxD * (r_old' * (r_old - r)) ./ ((r - r_old)' * (r - r_old));       
            relaxD = min(max(relaxD, obj.relaxDmin), obj.relaxDmax);
             
            % error on displacament
            obj.errD = max(abs((Dsp - Dsp_old) ./ Dsp));
            % The displacement is underelaxated
            Dsp = (1 - relaxD) * Dsp_old + relaxD * Dsp;
            h_new = obj.h0 + domain.height_profile - Dsp;      %the h0 changes when each load balance occurs. 
            
            obj.r = r';
            obj.displacement = Dsp;
            
            if obj.errD < obj.tollD
                dspflag = false;
            else
                obj.h = h_new;
                dspflag = true;
            end
            
        end
            
        function [obj] = set_hertizan_pressure(obj, F)
            obj.hertzian_pressure = 2 * F / (pi * obj.domain.bH^2) .* ((obj.domain.bH)^2 - obj.domain.x.^2).^0.5;
        end
        
        function [obj] = Reynolds(obj)
            
            domain = obj.domain;
              
            [dhdx, a, c, u, cn, cs] = ...
                thicknessDer(obj.h, domain.dx);
            
            [pressure] = FluidPressISOSliding(domain.N_x, domain.mu, obj.h, dhdx, cn, cs, u, ...
                obj.p_in, obj.pressure, obj.p_out, obj.p_cavitation, obj.velocity , obj.relaxp, obj.tollp, 1000000);
            
            obj.pressure = pressure;
            
        end
        
        function [obj] = LIN(obj, cavitation_sting, cavitation_flag, De_fake)
            
            domain = obj.domain; 
            De = domain.epsilon * obj.wiessenberg_Number; 
            
            if De_fake ~= 0
                De = De_fake;
            end
            
            if cavitation_sting == "cavitation"
                if cavitation_flag == "on"
                    
                    
                    [pressure, flowrateT, flowrate, flowrateD, thetaT, theta, thetaD, pressure0, pressure1] = ...
                        LIN_EHL_cavitation(domain.N_x, domain.mu, obj.h, obj.p_in, obj.pressure, obj.p_out, obj.p_cavitation, obj.velocity, obj.tollp, domain.L, 10000000, ...
                        De, obj.viscocity_ratio, obj.pressure0_old, obj.pressure1_old, obj.theta_old, obj.thetaT_old, obj.thetaD_old, domain.href);
                    
                    obj.pressure = pressure;
                    obj.pressure0_old = pressure0;
                    obj.pressure1_old = pressure1;
                    
                    obj.theta_old = theta;
                    obj.thetaT_old = thetaT;
                    obj.thetaD_old = thetaD;
                   
                    
                elseif cavitation_flag == "off"
                    
                    [pressure, pressure0, pressure1] = ...
                        LIN_EHL(domain.N_x, domain.mu, obj.h, obj.p_in, obj.pressure, obj.p_out, obj.p_cavitation, obj.velocity, obj.tollp, domain.L, 1000000, ...
                        De, obj.viscocity_ratio, domain.href);
                    
                    obj.pressure = pressure;
                    obj.pressure0_old = pressure0;
                    obj.pressure1_old = pressure1;
                end
            end
            


        end

        function [obj] = FieldProperties(obj, N_x, N_y)
            [obj.friction, obj.velocity_field, obj.stress_field, obj.force] = StressVelocity_LIN(obj, N_x, N_y);
        end
        
        
        
        function [obj] = VR(obj, dt, N_x, N_y, Spatial_Scheme, Temporal_Scheme, De_fake, cavitation_flag)
            
            domain = obj.domain;
            De = domain.epsilon * obj.wiessenberg_Number;
            
            if De_fake ~= 0
                De = De_fake
            end
            
            [J_temp, N_temp] = size(obj.velocity_field{1,1});
            
            if J_temp ~= N_y || N_temp ~= N_x
                [obj.friction, obj.velocity_field, obj.stress_field] = StressVelocity_LIN(obj, N_x, N_y);
            end
            
%             stress = {zeros(N_y, N_x), zeros(N_y, N_x), zeros(N_y, N_x)};
            
            [pressure_new, theta_new, stress_new, obj.pressure_field] = ...
                VRCavitation(dt, N_x, N_y, Spatial_Scheme, Temporal_Scheme, obj.h, ...
            domain.mu, obj.p_in, obj.pressure, obj.p_cavitation, obj.velocity, obj.tollp, domain.L, 1000000,...
            De, obj.viscocity_ratio, obj.stress_field, obj.velocity_field, obj.pressure0_old, obj.pressure1_old, 1, ...
            domain.href, cavitation_flag);
            
            obj.pressure_VR = pressure_new;
            obj.stress_VR = stress_new;
            obj.theta_VR = theta_new;
            
        end
        
    end
    
end