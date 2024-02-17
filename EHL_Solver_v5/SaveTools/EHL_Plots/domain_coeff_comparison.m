function domain_coeff_comparison(solutions_different_coeff)
%solution_different_coeff  includes same parameters but different bH solution
%class. The solution_matrix includes different bH solutions. 
% size of the matrix is {k x l x m x n} where n is the number of different
% bH solutions. 

if ~isfolder("Results/Figures")
       mkdir("Results/Figures")
end

savepath = "Results\Figures\";

if length(solutions_different_coeff) == 2
    solution_matrix  = cat(4, solutions_different_coeff{1} , solutions_different_coeff{2});
end

if length(solutions_different_coeff) == 3
    solution_matrix  = cat(4, solutions_different_coeff{1} , solutions_different_coeff{2}, solutions_different_coeff{3});
end

[k, l, m, n] = size(solution_matrix); 

if ~isfolder(savepath + "EHL_domain_coeff_comparison")
       mkdir(savepath + "EHL_domain_coeff_comparison")
end

savepath_in = savepath + "EHL_domain_coeff_comparison\";


for j = 1:l
    
    if ~isfolder(savepath_in + string(solution_matrix(1,1,j,1).velocity))
       mkdir(savepath_in + string(solution_matrix(1,j,1,1).velocity))
    end
    
    savepath_2in = savepath_in +  string(solution_matrix(1,j,1,1).velocity) + "\";
    
    legend_list = [];
    applied_load_list = [solution_matrix(:,j,1,1).applied_load];
    
    for i = 1:length(applied_load_list)
    	legend_list = [legend_list, "\bf{L} " + string(applied_load_list(i))];
    end
    
    main_legend_list = [];
    coeff_list = "";
    
     wi_index = randi(4);
    
    for i = 1:n
        
        wi_number = domain_coeff_comparison_plots(solution_matrix(:, j, :, i),  i, savepath_2in, wi_index);
        
        new_legend_list =  ["\bf{"  + string(solution_matrix(1, j, 1, i).domain.domain_coeff) +  "bH}  " + legend_list(1), legend_list(2:end)]; 
        
        coeff_list = coeff_list + "_" + string(solution_matrix(1, j, 1, i).domain.domain_coeff);
        
        main_legend_list = [main_legend_list, new_legend_list];
    end
   
    pressure_plot = figure(1);
    gap_plot = figure(2);
   
   figure(1); legend(main_legend_list,"box","off","Interpreter","latex",'Orientation','horizontal', 'NumColumns',5); 
   figure(2); legend(main_legend_list,"box","off","Interpreter","latex",'Orientation','horizontal', 'NumColumns',5); 
   
   

   saveas(pressure_plot, savepath_2in + "svg\pressure\" + "Wi" + string(wi_number) +   "_beta" + string(solution_matrix(1, j, 1, 1).viscocity_ratio) + string(coeff_list) + ".svg")
   saveas(gap_plot, savepath_2in + "svg\gap\" + "Wi" + string(wi_number) +   "_beta" + string(solution_matrix(1, j, 1, 1).viscocity_ratio) + string(coeff_list) + ".svg")
   
   saveas(pressure_plot, savepath_2in + "fig\pressure\" + "Wi" + string(wi_number) +   "_beta" + string(solution_matrix(1, j, 1, 1).viscocity_ratio) + string(coeff_list) + ".fig")
   saveas(gap_plot, savepath_2in + "fig\gap\" + "Wi" + string(wi_number) +   "_beta" + string(solution_matrix(1, j, 1, 1).viscocity_ratio) + string(coeff_list) + ".fig")
   
   close(pressure_plot)
   close(gap_plot)

end
end