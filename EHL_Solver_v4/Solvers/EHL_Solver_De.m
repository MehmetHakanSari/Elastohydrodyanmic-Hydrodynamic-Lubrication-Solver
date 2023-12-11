function EHL_Solver_De(velocity_list, de_list, beta_list, load_list, h_cent, savenote)

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
        ehl_channel(i, j) = ehl_channel(i, j).set_epsilon(h_cent(j,  i), "None");       
%         ehl_channel(i, j) = ehl_channel(i, j).set_epsilon("min", "None");      
%         ehl_channel(i, j) = ehl_channel(i, j).set_epsilon("curvature", "None");   
    end
end
         
ehl_solution(length(load_list), length(velocity_list), length(de_list)) = solution;

for velocity = 1:length(velocity_list)   %Number of Velocities
    Ux = velocity_list(velocity);
for force = 1:length(load_list) 
    F = load_list(force);
for de = 1:length(de_list)
    De = de_list(de);
for bi = 1:length(beta_list)
    beta = beta_list(bi);

    
    loadflag = true;
    dispflag = true;
    
    ehl_solution(force, velocity, de).domain = ehl_channel(force, velocity);
    ehl_solution(force, velocity, de) = ehl_solution(force, velocity, de).set_solution_parameters(De, beta, Ux, F);
    
    disp("Solving for U = " + string(Ux) + " De = " + string(De) + " beta = " + string(beta) + " load = " + string(F))
    
    if Ux < 0.0007623 && F >= 100
        ehl_solution(force, velocity, de).relaxD = 1e-5;
        ehl_solution(force, velocity, de).relaxDmax = 1e-3;
        ehl_solution(force, velocity, de).relaxDmin = 1e-4;
    end    
    if Ux > 0.001 && F <= 6
        ehl_solution(force, velocity, de).relaxL = 0.1;
        ehl_solution(force, velocity, de).relaxD = 0.3;
    end
    
    
    ehl_solution(force, velocity, de) = ehl_solution(force, velocity, de).initilizer();
    ehl_solution(force, velocity, de) = ehl_solution(force, velocity, de).set_hertizan_pressure(F);
    
    diplacement_counter = 0;
    load_counter = 0;
    
    tic
    while loadflag 
        load_counter = load_counter + 1;
        while dispflag
            diplacement_counter = diplacement_counter + 1;
            
%             [ehl_solution(force, velocity, de)] = ehl_solution(force, velocity, de).Reynolds();
            [ehl_solution(force, velocity, de)] = ehl_solution(force, velocity, de).LIN("cavitation", "on", De);
%             [ehl_solution(force, l, de), pressure] = ehl_solution(force, l, wi).VR(domain, h, Wi, beta, Ux);
            
            [ehl_solution(force, velocity, de), dispflag] = ehl_solution(force, velocity, de).displacementsolver();
            
            if mod(diplacement_counter, 100) == 0
                disp("Error: " + string(ehl_solution(force, velocity, de).errD) + " at itteration " + string(diplacement_counter) + " Load balance at: " + string(load_counter))
            end
            
%             if diplacement_counter == 1
%                 break
%             end
            
        end
        
        [ehl_solution(force, velocity, de), loadflag] = ehl_solution(force, velocity, de).loadsolver(F);
        dispflag = true;
        diplacement_counter = 0;
        
%         if diplacement_counter == 1
%             break
%         end
        disp("-------------     v   -----------")
        disp("Error: " + string(ehl_solution(force, velocity, de).errL) + " Load balance at: " + string(load_counter) + " U = " + string(Ux) + " De = " + string(De) + " L = " + string(F))
        disp("--------------------------------------")

    end
    
    ehl_solution(force, velocity, de) = ehl_solution(force, velocity, de).FieldProperties(501,120);

    disp("the EHL is solved")
    toc
    disp("           ")
    disp("           ")
    
end
end
end
end


%Save the data
ehl_savedata(load_list, beta_list, de_list, velocity_list, ehl_solution, "De", savenote)


end