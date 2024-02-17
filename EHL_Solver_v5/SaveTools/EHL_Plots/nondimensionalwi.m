function nondimensionalwi(solution_matrix, savepath)

[k, l, m] = size(solution_matrix); %force, velocity, wi

load_list = [solution_matrix(:,1,1).applied_load];  
beta_list = [solution_matrix(1,1,1).viscocity_ratio]; 
if min([solution_matrix(1,1,:).wiessenberg_Number]) < 8
    wi_list = [solution_matrix(1,1,:).deborah_Number]; 
else
    wi_list = [solution_matrix(1,1,:).wiessenberg_Number]; 
end
velocity_list = [solution_matrix(1,:,1).velocity];

if min([solution_matrix(1,1,:).wiessenberg_Number]) < 8
    file_name_string = "decompND_F";
else
    file_name_string = "wicompND_F";
end


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
    
    for i = 1:k
        if min([solution_matrix(1,1,:).wiessenberg_Number]) < 8
            EHLCompareDe_NoNDimensional(solution_matrix(i, j, :), savepath_2in)
        else
            EHLCompareWi_NoNDimensional(solution_matrix(i, j, :), savepath_2in)
        end
    end
   
   legend_list = [];
   %Legend and close h_min and friction
    for p = 1:length(solution_matrix(:,1,1))
        legend_list = [legend_list, "\bf{W = " + string(solution_matrix(p,1,1).applied_load  / solution_matrix(p,1,1).domain.Rx / solution_matrix(p,1,1).domain.E / solution_matrix(p,1,1).domain.L) + "}"];
    end
   figure(3);
   title("\bf{U = }" + string(solution_matrix(1,j,1).velocity * solution_matrix(1,j,1).domain.mu / solution_matrix(1,j,1).domain.Rx / solution_matrix(1,j,1).domain.E),"Interpreter","latex")
   legend(legend_list,"box","off","Interpreter","latex", "location","best")
   figure(4);
   legend(legend_list,"box","off","Interpreter","latex", "location","best")
   figure(9);
   legend(legend_list,"box","off","Interpreter","latex", "location","best")
   figure(10);
   legend(legend_list,"box","off","Interpreter","latex", "location","best")
   figure(11);
   legend(legend_list,"box","off","Interpreter","latex", "location","best")
   figure(12);
   legend(legend_list,"box","off","Interpreter","latex", "location","best")
   
   minh_plot = figure(3);
   friction_plot = figure(4);
   pNNpT_plot = figure(9);
   f_plot = figure(10);
   f_b = figure(11);
   f_b_v = figure(12);
   f_b_n = figure(13);
   
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
   saveas(f_plot, savepath_in + "friction\" +  "f_beta" + string(solution_matrix(1,j,1).viscocity_ratio) +  "U" + string(solution_matrix(1,j,1).velocity) + ".svg")
   saveas(f_plot, savepath_in + "friction\" +  "f_beta" + string(solution_matrix(1,j,1).viscocity_ratio) +  "U" + string(solution_matrix(1,j,1).velocity) + ".fig")
   
   saveas(f_b, savepath_in + "friction\" +  "f_b" + string(solution_matrix(1,j,1).viscocity_ratio) +  "U" + string(solution_matrix(1,j,1).velocity) + ".fig")
   saveas(f_b, savepath_in + "friction\" +  "f_b" + string(solution_matrix(1,j,1).viscocity_ratio) +  "U" + string(solution_matrix(1,j,1).velocity) + ".svg")
   saveas(f_b_v, savepath_in + "friction\" +  "f_b_v" + string(solution_matrix(1,j,1).viscocity_ratio) +  "U" + string(solution_matrix(1,j,1).velocity) + ".fig")
   saveas(f_b_v, savepath_in + "friction\" +  "f_b_v" + string(solution_matrix(1,j,1).viscocity_ratio) +  "U" + string(solution_matrix(1,j,1).velocity) + ".svg")
   saveas(f_b_n, savepath_in + "friction\" +  "f_b_n" + string(solution_matrix(1,j,1).viscocity_ratio) +  "U" + string(solution_matrix(1,j,1).velocity) + ".fig")
   saveas(f_b_n, savepath_in + "friction\" +  "f_b_n" + string(solution_matrix(1,j,1).viscocity_ratio) +  "U" + string(solution_matrix(1,j,1).velocity) + ".svg")

   
   
   if ~isfolder(savepath_in + "pNNpT")
       mkdir(savepath_in + "pNNpT")
   end
   
   saveas(pNNpT_plot, savepath_in + "pNNpT\" +  "beta" + string(solution_matrix(1,j,1).viscocity_ratio) +  "U" + string(solution_matrix(1,j,1).velocity) + ".svg")
   saveas(pNNpT_plot, savepath_in + "pNNpT\" +  "beta" + string(solution_matrix(1,j,1).viscocity_ratio) +  "U" + string(solution_matrix(1,j,1).velocity) + ".fig")
   
   close all
   

end
end