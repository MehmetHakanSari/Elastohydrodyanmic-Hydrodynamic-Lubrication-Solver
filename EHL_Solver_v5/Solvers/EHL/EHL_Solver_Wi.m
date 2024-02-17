function EHL_Solver_Wi(velocity_list, wi_list, beta_list, load_list, savenote)

ehl_channel(length(load_list), length(velocity_list)) = domain;

for i = 1:length(load_list)
    for j = 1:length(velocity_list)
        ehl_channel(i, j) = ehl_channel(i, j).set_domain_coeff(velocity_list(j));   
        ehl_channel(i, j) = ehl_channel(i, j).set_eq_radii_curvature();
        ehl_channel(i, j) = ehl_channel(i, j).set_eq_elastic_modules(); 
        ehl_channel(i, j) = ehl_channel(i, j).set_bH(load_list(i));
        ehl_channel(i, j) = ehl_channel(i, j).set_ehl_length(velocity_list(j));
        ehl_channel(i, j) = ehl_channel(i, j).elastic_map();
        ehl_channel(i, j) = ehl_channel(i, j).set_initial_height("punch");
        ehl_channel(i, j) = ehl_channel(i, j).set_epsilon(h_cent(j), "None");       
    end
end
         
ehl_solution(length(load_list), length(velocity_list), length(wi_list)) = solution;

for velocity = 1:length(velocity_list)   %Number of Velocities
    Ux = velocity_list(velocity);
for force = 1:length(load_list) 
    F = load_list(force);
for wi = 1:length(wi_list)
    Wi = wi_list(wi);
for bi = 1:length(beta_list)
    beta = beta_list(bi);

    
    loadflag = true;
    dispflag = true;
    
    ehl_solution(force, velocity, wi).domain = ehl_channel(force, velocity);
    ehl_solution(force, velocity, wi) = ehl_solution(force, velocity, wi).set_solution_parameters(Wi, beta, Ux, F);
    
    disp("Solving for U = " + string(Ux) + " Wi = " + string(Wi) + " beta = " + string(beta) + " load = " + string(F))
    
    if velocity == 0.0001 && F >= 100
        ehl_solution(force, velocity, wi).relaxD = 1e-5;
        ehl_solution(force, velocity, wi).relaxDmax = 1e-3;
        ehl_solution(force, velocity, wi).relaxDmin = 1e-4;
    end
    
    ehl_solution(force, velocity, wi) = ehl_solution(force, velocity, wi).initilizer();
    ehl_solution(force, velocity, wi) = ehl_solution(force, velocity, wi).set_hertizan_pressure(F);
    
    diplacement_counter = 0;
    load_counter = 0;
    
    tic
    while loadflag 
        load_counter = load_counter + 1;
        while dispflag
            diplacement_counter = diplacement_counter + 1;
            
%             [ehl_solution(force, velocity, wi)] = ehl_solution(force, velocity, wi).Reynolds();
            [ehl_solution(force, velocity, wi)] = ehl_solution(force, velocity, wi).LIN("cavitation", "on", 0);
%             [ehl_solution(force, l, wi), pressure] = ehl_solution(force, l, wi).VR(domain, h, Wi, beta, Ux);
            
            [ehl_solution(force, velocity, wi), dispflag] = ehl_solution(force, velocity, wi).displacementsolver();
            
            if mod(diplacement_counter, 100) == 0
                disp("Error: " + string(ehl_solution(force, velocity, wi).errD) + " at itteration " + string(diplacement_counter) + " Load balance at: " + string(load_counter))
            end
            
%             if diplacement_counter == 1
%                 break
%             end
            
        end
        
        [ehl_solution(force, velocity, wi), loadflag] = ehl_solution(force, velocity, wi).loadsolver(F);
        dispflag = true;
        diplacement_counter = 0;
        
%         if diplacement_counter == 1
%             break
%         end
        disp("-------------     v   -----------")
        disp("Error: " + string(ehl_solution(force, velocity, wi).errL) + " Load balance at: " + string(load_counter) + " U = " + string(Ux) + " Wi = " + string(Wi) + " L = " + string(F))
        disp("--------------------------------------")

    end
    
    ehl_solution(force, velocity, wi) = ehl_solution(force, velocity, wi).FieldProperties(501,120);

    disp("the EHL is solved")
    toc
    disp("           ")
    disp("           ")
    
end
end
end
end


%Save the data
ehl_savedata(load_list, beta_list, wi_list, velocity_list, ehl_solution, "Wi", savenote)


end