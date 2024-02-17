function dimensionalload(solution_matrix, savepath)

[k, l, m] = size(solution_matrix); %force, velocity, wi

load_list = [solution_matrix(:,1,1).applied_load];  
beta_list = [solution_matrix(1,1,1).viscocity_ratio]; 

if min([solution_matrix(1,1,:).wiessenberg_Number]) < 8
    wi_list = [solution_matrix(1,1,:).deborah_Number]; 
else
    wi_list = [solution_matrix(1,1,:).wiessenberg_Number]; 
end

velocity_list = [solution_matrix(1,:,1).velocity];

file_name_string = "loadcompD_F";

for i = 1:length(load_list)
    file_name_string = file_name_string + "-" + string(load_list(i));
end
file_name_string = file_name_string + "_beta";
for i = 1:length(beta_list)
    file_name_string = file_name_string + "-" + string(beta_list(i));
end
if min([solution_matrix(1,1,:).wiessenberg_Number]) < 8
    file_name_string = file_name_string + "_De";
else
    file_name_string = file_name_string + "_Wi";
end
for i = 1:length(wi_list)
    file_name_string = file_name_string + "-" + string(wi_list(i));
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
    
    if ~isfolder(savepath_in + string(solution_matrix(1,j,1).velocity))
       mkdir(savepath_in + string(solution_matrix(1,j,1).velocity))
    end
    
    savepath_2in = savepath_in + "\" + string(solution_matrix(1,j,1).velocity) + "\";
    
    for i = 1:m
        EHLCompareL_Dimensional(solution_matrix(:,j,i), savepath_2in)
    end
   
   legend_list = [];
   %Legend and close h_min and friction
   for p = 1:length(solution_matrix(1,1,:))
       if min([solution_matrix(1,1,:).wiessenberg_Number]) < 8
            legend_list = [legend_list, "\bf{De = " + string(solution_matrix(1,1,p).deborah_Number) + "}"];
       else
            legend_list = [legend_list, "\bf{Wi = " + string(solution_matrix(1,1,p).wiessenberg_Number) + "}"];
       end
   end
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
   saveas(minh_plot, savepath_in + "hmin\fminh" +  "beta" + string(solution_matrix(1,j,1).viscocity_ratio) + "U" + string(solution_matrix(1,j,1).velocity) + ".fig")

      
   if ~isfolder(savepath_in + "friction")
       mkdir(savepath_in + "friction")
   end
   
   saveas(friction_plot, savepath_in + "friction\" +  "beta" + string(solution_matrix(1,j,1).viscocity_ratio) +  "U" + string(solution_matrix(1,j,1).velocity) + ".svg")
   saveas(friction_plot, savepath_in + "friction\" +  "beta" + string(solution_matrix(1,j,1).viscocity_ratio) +  "U" + string(solution_matrix(1,j,1).velocity) + ".fig")
   
   close all
   
end
end