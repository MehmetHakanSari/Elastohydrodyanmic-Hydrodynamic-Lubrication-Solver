function valueplot_figuregenerator(solution,  demand)
%solution is one solution class family. 
%demand is demanded figures

%demand options:
% "all" - "pressure" - "theta" - "massflow"


if ~isfolder("Results/Figures")
       mkdir("Results/Figures")
end

savepath = "Results\Figures\";

l_demand = length(demand);

for type_index = 1:l_demand
    %plot different loads for same velocity and Wi
   type = demand(type_index);
   if type == "all"
       pressureplot(solution, savepath);
       close all 
       heightplot(solution, savepath)
       close all 
       stressfieldplot(solution, savepath)
       close all
       velocityfieldplot(solution, savepath)
       close all
   end
   
   if type == "pressure_D"
       pressureplot(solution, savepath);
       if (l_demand > 1)
           close all
       end
   end
   
   if type == "pressure_ND"
       if (l_demand > 1)
           close all
       end
       
   end
   
   if type == "theta"
       if (l_demand > 1)
           close all
       end
   end
   
   if type == "massflow"
       if (l_demand > 1)
           close all
       end
   end
   
    if type == "height"
        heightplot(solution, savepath)
       if (l_demand > 1)
           close all
       end
   end
   
   if type == "stressfield"
       stressfieldplot(solution, savepath)
       if (l_demand > 1)
           close all
       end
   end
   
    if type == "velocityfield"
       velocityfieldplot(solution, savepath)
       if (l_demand > 1)
           close all
       end
   end
       
%plot different wi for same velocity and load


%plot different loads for same velocity and Wi

end
end