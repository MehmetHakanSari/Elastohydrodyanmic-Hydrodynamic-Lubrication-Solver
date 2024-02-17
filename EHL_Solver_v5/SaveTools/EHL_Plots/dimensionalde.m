function dimensionalde(solution_matrix, savepath)

[k, l, m] = size(solution_matrix); %force, velocity, wi

load_list = [solution_matrix(:,1,1).applied_load];  
beta_list = [solution_matrix(1,1,1).viscocity_ratio]; 
de_list = [solution_matrix(1,1,:).deborah_Number]; 
velocity_list = [solution_matrix(1,:,1).velocity];

file_name_string = "decompD_F";

for i = 1:length(load_list)
    file_name_string = file_name_string + "-" + string(load_list(i));
end
file_name_string = file_name_string + "_beta";
for i = 1:length(beta_list)
    file_name_string = file_name_string + "-" + string(beta_list(i));
end
file_name_string = file_name_string + "_De";
for i = 1:length(de_list)
    file_name_string = file_name_string + "-" + string(de_list(i));
end
file_name_string = file_name_string + "_U";
for i = 1:length(velocity_list)
    file_name_string = file_name_string + "-" + string(velocity_list(i));
end

if ~isfolder(savepath + "file_name_string")
       mkdir(savepath + "file_name_string")
end

savepath_in = savepath + "\" + file_name_string  + "\";

for j = 1:l
    
    if ~isfolder(savepath_in + string(solution_matrix(1, j ,1).velocity))
       mkdir(savepath_in + string(solution_matrix(1, j, 1).velocity))
    end
    
    savepath_2in = savepath_in + "\" + string(solution_matrix(1,j,1).velocity) + "\";
    
    for i = 1:k
        EHLCompareWi_Dimensional(solution_matrix(k, j, :), savepath_2in)
    end
   
    
    legend_list = [];
    
    %Legend and close h_min and friction
    for p = 1:length(solution_matrix(:,1,1))
        legend_list = [legend_list, "\bf{L = " + string(solution_matrix(1,1,p).applied_load) + "}"];
    end
%    legend_list = ["\bf{Wi = 0}", "\bf{Wi = 10}", "\bf{Wi = 20}", "\bf{Wi = 40}"];
   figure(3);
   title("\bf{U = }" + string(solution_matrix(1,1,j).velocity),"Interpreter","latex")
   legend(legend_list,"box","off","Interpreter","latex", "location","best")
   figure(4);
   legend(legend_list,"box","off","Interpreter","latex", "location","best")
   
   minh_plot = figure(3);
   friction_plot = figure(4);
   
   if ~isfolder(savepath_in + "hmin")
       mkdir(savepath_in + "hmin")
   end

   saveas(minh_plot, savepath_in + "hmin\minh" +  "beta" + string(solution_matrix(1,j,1).viscocity_ratio) + "U" + string(solution_matrix(1,j,1).velocity) + ".svg")
   saveas(minh_plot, savepath_in + "hmin\minh" +  "beta" + string(solution_matrix(1,j,1).viscocity_ratio) + "U" + string(solution_matrix(1,j,1).velocity) + ".fig")
   
   if ~isfolder(savepath_in + "friction")
       mkdir(savepath_in + "friction")
   end

   saveas(friction_plot, savepath_in + "friction\" +  "beta" + string(solution_matrix(1,j,1).viscocity_ratio) +  "U" + string(solution_matrix(1,j,1).velocity) + ".svg")
   saveas(friction_plot, savepath_in + "friction\" +  "beta" + string(solution_matrix(1,j,1).viscocity_ratio) +  "U" + string(solution_matrix(1,j,1).velocity) + ".fig")
   
   close all
   
end
end