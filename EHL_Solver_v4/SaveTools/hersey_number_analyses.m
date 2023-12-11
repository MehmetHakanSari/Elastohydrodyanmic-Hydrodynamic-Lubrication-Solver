function hersey_number_analyses(data, savepath)

[k, l, m] = size(data); %force, velocity, wi

if ~isfolder(savepath + "HerseyNumberAnalyses")
       mkdir(savepath + "HerseyNumberAnalyses")
end

savepath_in = savepath + "HerseyNumberAnalyses"  + "\";

ms = 16; malpha = 0.9;
de_index = [1,2,4];
de_marker = ["v", "o", "square"];

de_color = ["#0c2c84","#225ea8","#1d91c0"];

for de = 1:length(de_index)
    
    for j = 1:k

       f_coeff = [];
       he_num = [];

       for i = 1:l
       
       f_coeff = [f_coeff data(j, i, de).friction{4}];
       he_num = [he_num data(j, i, de).hersey_Number];
            
       end   
   
       figure(1); hold on, box on;
       he_num
       lh = plot(he_num, f_coeff,de_marker(de),"Color",de_color(de),"MarkerSize", ms, "MarkerFaceColor", de_color(de)); lh.Color(4)=malpha;
       
    end

end

figure(1); box on;
xlabel("\boldmath{$He$}", "Interpreter","latex")
ylabel("\boldmath{friction}  ","Interpreter","latex")
title(" ","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
set(gca,'yscale','log','xscale','log')
set(gca, "Linewidth", 1.1)
   
stibeckcurve = figure(1); 
set(stibeckcurve,'Position',[10 90 850 720])
   
saveas(stibeckcurve, savepath_in + "stibeckcurve" +  "beta" + string(data(1,1,j).viscocity_ratio) + ".svg")

end