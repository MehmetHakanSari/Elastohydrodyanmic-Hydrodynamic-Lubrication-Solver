function ehl_savedata(load_list, beta_list, lineraztion_list, velocity_list, ehl_solution, relaxation_param, messeage)
%load list: vector of dobules
%beta list: vector of dobules
%wi list: vector of dobules
%velocity list: vector of dobules
%message: string of stated message 

if ~isfolder("Results")
       mkdir("Results")
       if ~isfolder("Results\EHL_Results")
           mkdir("Results\EHL_Results")
       end
end

file_name_string = "F";
for i = 1:length(load_list)
    file_name_string = file_name_string + "_" + string(load_list(i));
end
file_name_string = file_name_string + "_beta";
for i = 1:length(beta_list)
    file_name_string = file_name_string + "_" + string(beta_list(i));
end

if relaxation_param == "Wi"
    file_name_string = file_name_string + "_Wi";
    for i = 1:length(lineraztion_list)
        file_name_string = file_name_string + "_" + string(lineraztion_list(i));
    end
elseif relaxation_param == "De"
    file_name_string = file_name_string + "_De";
    for i = 1:length(lineraztion_list)
        file_name_string = file_name_string + "_" + string(lineraztion_list(i));
    end 
end

file_name_string = file_name_string + "_U";
for i = 1:length(velocity_list)
    file_name_string = file_name_string + "_" + string(velocity_list(i));
end

if i == 1
file_name_string = file_name_string + "_bH" + string(ehl_solution(1,1,1).domain.domain_coeff);
end

save("Results/EHL_Results/" + file_name_string + "_" + messeage + ".mat", "ehl_solution")


end