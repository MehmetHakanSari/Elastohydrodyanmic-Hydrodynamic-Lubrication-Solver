function EHLCompareDe_Dimensional(solution, savepath)
% solution contains the solution class EHL values for different loads but
% same velocity and wi. The length of the solution is equals the length of
% the results that will be plotted. 

savepath = string(savepath);
legend_list = [];
hmin_list = [];

f_t = [];
f_b = [];
f_v = [];
f_coeff = [];


de_list = [solution(1,1,:).deborah_Number];

newcolors = [0 0.4470 0.7410; 0 0.4470 0.7410;  ...
    0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; ...
    0.9290 0.6940 0.1250; 0.9290 0.6940 0.1250; ...
    0.4940 0.1840 0.5560; 0.4940 0.1840 0.5560; ...
    0.4660 0.6740 0.1880; 0.4660 0.6740 0.1880];

for i = 1:length(solution(1,1,:))
    
    figure(1); hold on
    plot([solution(1,1,i).domain.x] ./ [solution(1,1,i).domain.bH], [solution(1,1,i).pressure], "LineWidth",2)
    figure(2); hold on
    plot([solution(1,1,i).domain.x] ./ [solution(1,1,i).domain.bH], [solution(1,1,i).h], "LineWidth",2)
    figure(5); hold on
    plot([solution(1,1,i).domain.x] ./ [solution(1,1,i).domain.bH], [solution(1,1,i).pressure0_old], "LineWidth",1.7)
    figure(6); hold on
    plot([solution(1,1,i).domain.x] ./ [solution(1,1,i).domain.bH], [solution(1,1,i).pressure1_old * solution(1,1,i).domain.epsilon * solution(1,1,i).deborah_Number], "LineWidth",1.7)
    figure(7); hold on
    plot([solution(1,1,i).domain.x] ./ [solution(1,1,i).domain.bH], [OneDcentraldiff(solution(1,1,i).h, solution(1,1,i).domain.dx)],"--","LineWidth",1.7)
    figure(8); hold on
    plot([solution(1,1,i).domain.x] ./ [solution(1,1,i).domain.bH], [OneDcentraldiff(OneDcentraldiff(solution(1,1,i).h, solution(1,1,i).domain.dx),solution(1,1,i).domain.dx)], "--","LineWidth",1.7)
   
    
    h = solution(1,1,i).h;
    hmin_list(i) = min(h);
   
    legend_list = [legend_list, "\bf{De} " + string(de_list(i))];
    
    f_t = [f_t solution(1,1,i).friction{1}];
    f_b = [f_b solution(1,1,i).friction{2}];
    f_v = [f_v solution(1,1,i).friction{3}];
    f_coeff = [f_coeff solution(1,1,i).friction{4}];
    
%     friction_top_list = [friction_top_list solution(1,1,i).friction{1}];
%     friction_bottom_list = [friction_bottom_list solution(1,1,i).friction{2}];
%     friction_sum_list = [friction_sum_list solution(1,1,i).friction{1}+solution(1,1,i).friction{2}];
%     friction_visco_list = [friction_visco_list solution(1,1,i).friction{3}];
end

figure(3); hold on
plot(de_list,  hmin_list/(hmin_list(1)),".-","MarkerSize", 24)

figure(4);
subplot(2,2,1); hold on
set(gca, 'ColorOrder', newcolors);
plot(de_list, f_t, "^--","LineWidth",1.6,"MarkerSize", 24)
plot(de_list, f_b,"v--","LineWidth",1.6,"MarkerSize", 24)
subplot(2,2,2); hold on
plot(de_list, f_coeff,".--","LineWidth",1.6,"MarkerSize", 24)
subplot(2,2,3); hold on
plot(de_list, f_t + f_b, ".--","LineWidth",1.6,"MarkerSize", 24)
subplot(2,2,4); hold on
plot(de_list, f_v, ".--","LineWidth",1.6,"MarkerSize", 24)

figure(1); 
xlabel("\boldmath{$x/b_{H}$}","Interpreter","latex")
ylabel("\bf{p (Pa)}","Interpreter","latex")
title("\bf{U = }" + string(solution(1,1,1).velocity) + ", \bf{L = }" + string(solution(1,1,1).applied_load),"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend(legend_list,"box","off","Interpreter","latex")

figure(5); 
xlabel("\boldmath{$x/b_{H}$}","Interpreter","latex")
ylabel("\bf{p (Pa)}","Interpreter","latex")
title("\bf{$p^N$ \ U = }" + string(solution(1,1,1).velocity) + ", \bf{L = }" + string(solution(1,1,1).applied_load),"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend(legend_list,"box","off","Interpreter","latex")

figure(6); 
xlabel("\boldmath{$x/b_{H}$}","Interpreter","latex")
ylabel("\bf{p (Pa)}","Interpreter","latex")
title("\bf{$p^{NN}$ \ U = }" + string(solution(1,1,1).velocity) + ", \bf{L = }" + string(solution(1,1,1).applied_load),"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend(legend_list,"box","off","Interpreter","latex")

figure(7); 
xlabel("\boldmath{$x$}","Interpreter","latex")
ylabel("\boldmath{}", "Interpreter","latex")
title("\bf{$\frac{\partial h}{\partial x}$ \ U = }" + string(solution(1,1,1).velocity) + ", \bf{L = }" + string(solution(1,1,1).applied_load),"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend(legend_list,"box","off","Interpreter","latex")
xlim([-3 3])

figure(8); 
xlabel("\boldmath{$x$}","Interpreter","latex")
ylabel("boldmath{$\frac{1}{m}$}", "Interpreter","latex")
ylabel("\bf{p (Pa)}","Interpreter","latex")
title("\bf{$\frac{\partial^2 h}{\partial x^2}$ \ U = }" + string(solution(1,1,1).velocity) + ", \bf{L = }" + string(solution(1,1,1).applied_load),"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend(legend_list,"box","off","Interpreter","latex")
xlim([-3 3])


figure(2); 
xlabel("\boldmath{$x/b_{H}$}","Interpreter","latex")
xlim([-3 3])
ylabel("\bf{h (m)}","Interpreter","latex")
title("\bf{U} = " + string(solution(1,1,1).velocity) + ", \bf{L} = " + string(solution(1,1,1).applied_load) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend(legend_list,"box","off","Interpreter","latex")

figure(3); 
xlabel("\bf{Deborah}","Interpreter","latex")
ylabel("\boldmath{$h_{min} / h_{min, De = " + string(min(de_list)) +  "}$}","Interpreter","latex")
title("","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)

% legend(legend_list,"box","off")

figure(4); 
subplot(2,2,1); hold on
% xlabel("\bf{Wi}","Interpreter","latex")
% ylabel("\boldmath{$f_t / f_{t,L=" + string(force_list(1)) + "}$}","Interpreter","latex")
ylabel("\boldmath{$N/m$} ","Interpreter","latex")
title("\boldmath{$f_t \ and \ f_b$ },   \bf{U} = " + string(solution(1,1).velocity),"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)

subplot(2,2,2); hold on
% xlabel("\bf{L}","Interpreter","latex")
% ylabel("\boldmath{$f_b / f_{b,L=" + string(force_list(1)) + "}$}","Interpreter","latex")
ylabel("\boldmath{$\mu$}","Interpreter","latex")
title(" ","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)

subplot(2,2,3); hold on
xlabel("\bf{De}", "Interpreter","latex")
ylabel(" \boldmath{$N/m$}", "Interpreter","latex")
title("\boldmath{$f_t + f_b$} ","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)

subplot(2,2,4); hold on
xlabel("\bf{De}","Interpreter","latex")
ylabel("\boldmath{$N/m$}","Interpreter","latex")
title("\boldmath{$f_{viscoelastic}$} ","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)


pressure_plot = figure(1);
gap_plot = figure(2);
minh_plot = figure(3);
friction_plot = figure(4);
pressure0_plot = figure(5);
pressure1_plot = figure(6);
dgap_plot = figure(7);
d2gap_plot = figure(8);

set(pressure_plot,'Position',[10 10 980 720])
set(gap_plot,'Position',[10 10 980 720])
set(minh_plot,'Position',[35 35 1080 800])
set(friction_plot,'Position',[835 235 1600 1400])
set(pressure0_plot,'Position',[10 10 980 720])
set(pressure1_plot,'Position',[10 10 980 720])
set(dgap_plot,'Position',[10 10 980 720])
set(d2gap_plot,'Position',[10 10 980 720])

% savepath = 'E:\Bilkent Dökümanları\Masterımsı\Haftalık\Dec 15 22\';

if ~isfolder(savepath + "/gap")
    mkdir(savepath + "/gap")
end
if ~isfolder(savepath + "/pressure")
    mkdir(savepath + "/pressure")
end

saveas(pressure_plot, savepath + "pressure/" + "EHLPressureDe_" +  "beta" + string(solution(1,1).viscocity_ratio) + "L" + string(solution(1,1,1).applied_load) + "U" + string(solution(1,1).velocity) + ".svg")
saveas(gap_plot, savepath + "gap/" + "EHLGapDe_" +  "beta" + string(solution(1,1).viscocity_ratio) + "L" + string(solution(1,1,1).applied_load) + "U" + string(solution(1,1).velocity) + ".svg")
saveas(dgap_plot, savepath + "gap/" + "dGapDe_" +  "beta" + string(solution(1,1).viscocity_ratio) + "L" + string(solution(1,1,1).applied_load) + "U" + string(solution(1,1).velocity) + ".svg")
saveas(d2gap_plot, savepath + "gap/" + "d2GapDe_" +  "beta" + string(solution(1,1).viscocity_ratio) + "L" + string(solution(1,1,1).applied_load) + "U" + string(solution(1,1).velocity) + ".svg")
saveas(pressure0_plot, savepath + "pressure/" + "EHLPressureNDe_" +  "beta" + string(solution(1,1).viscocity_ratio) + "L" + string(solution(1,1,1).applied_load) + "U" + string(solution(1,1).velocity) + ".svg")
saveas(pressure1_plot, savepath + "pressure/" + "EHLPressureNNDe_" +  "beta" + string(solution(1,1).viscocity_ratio) + "L" + string(solution(1,1,1).applied_load) + "U" + string(solution(1,1).velocity) + ".svg")

saveas(pressure_plot, savepath + "pressure/" + "EHLPressureDe_" +  "beta" + string(solution(1,1).viscocity_ratio) + "L" + string(solution(1,1,1).applied_load) + "U" + string(solution(1,1).velocity) + ".fig")
saveas(gap_plot, savepath + "gap/" + "EHLGapDe_" +  "beta" + string(solution(1,1).viscocity_ratio) + "L" + string(solution(1,1,1).applied_load) + "U" + string(solution(1,1).velocity) + ".fig")
saveas(dgap_plot, savepath + "gap/" + "dGapDe_" +  "beta" + string(solution(1,1).viscocity_ratio) + "L" + string(solution(1,1,1).applied_load) + "U" + string(solution(1,1).velocity) + ".fig")
saveas(d2gap_plot, savepath + "gap/" + "d2GapDe_" +  "beta" + string(solution(1,1).viscocity_ratio) + "L" + string(solution(1,1,1).applied_load) + "U" + string(solution(1,1).velocity) + ".fig")
saveas(pressure0_plot, savepath + "pressure/" + "EHLPressureNDe_" +  "beta" + string(solution(1,1).viscocity_ratio) + "L" + string(solution(1,1,1).applied_load) + "U" + string(solution(1,1).velocity) + ".fig")
saveas(pressure1_plot, savepath + "pressure/" + "EHLPressureNNDe_" +  "beta" + string(solution(1,1).viscocity_ratio) + "L" + string(solution(1,1,1).applied_load) + "U" + string(solution(1,1).velocity) + ".fig")


close(pressure_plot)
close(gap_plot)
close(dgap_plot)
close(d2gap_plot)
close(pressure0_plot)
close(pressure1_plot)


% U_nondimensional =
% 
%    1.8750e-11
%                   for U = 0.0001, mu = 0.1, E1 = 1e7, v = 0.5, Rx = 0.02

% force_nondmensional =  F / Rx / 20 / bH / E 
% 
% ans =
% 
%     0.0136

end