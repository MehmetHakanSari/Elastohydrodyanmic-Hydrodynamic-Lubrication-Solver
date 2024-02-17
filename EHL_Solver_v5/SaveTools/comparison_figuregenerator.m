function comparison_figuregenerator(solution_matrix,  demand)

if ~isfolder("Results/Figures")
       mkdir("Results/Figures")
end

savepath = "Results\Figures\";

for type_index = 1:length(demand)
    %plot different loads for same velocity and Wi
   type = demand(type_index);
   if type == "load_D"
      dimensionalload(solution_matrix, savepath)        
   end
   
   if type == "load_ND"
      nondimensionalload(solution_matrix, savepath)   
   end
   
   if type == "velocity_D"
       dimensionalvelocity(solution_matrix, savepath)       
   end
   
   if type == "velocity_ND"
       nondimensionalvelocity(solution_matrix, savepath)  
   end
   
   if type == "Wiessenberg_D"
       dimensionalwi(solution_matrix, savepath)        
   end
   
   if type == "Wiessenberg_ND"
       nondimensionalwi(solution_matrix, savepath)     
   end
   
   if type == "beta_ND"
       nondimensionalbeta(solution_matrix, savepath)     
   end
   
   if type == "StibeckCurve"
       hersey_number_analyses(solution_matrix, savepath)     
   end
   
end
end