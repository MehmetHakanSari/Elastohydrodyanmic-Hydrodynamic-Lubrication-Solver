function EHLCompareL_Dimensional(solution, savepath)
% solution contains the solution class EHL values for different loads but
% same velocity and wi. The length of the solution is equals the length of
% the results that will be plotted. 


savepath = string(savepath);
legend_list = [];
hmin_list = [];
applied_load_list = [solution(:,1).applied_load];

f_t = [];
f_b = [];
f_v = [];
f_coeff = [];

newcolors = [0 0.4470 0.7410; 0 0.4470 0.7410;  ...
    0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; ...
    0.9290 0.6940 0.1250; 0.9290 0.6940 0.1250; ...
    0.4940 0.1840 0.5560; 0.4940 0.1840 0.5560; ...
    0.4660 0.6740 0.1880; 0.4660 0.6740 0.1880];


for i = 1:length(solution(:,1))
    figure(1); hold on
    plot([solution(i,1).domain.x] ./ [solution(i,1).domain.bH], [solution(i,1).pressure], "LineWidth",2)
    figure(2); hold on
    plot([solution(i,1).domain.x] ./ [solution(i,1).domain.bH], [solution(i,1).h], "LineWidth",2)
    
    h = solution(i,1).h;
    hmin_list(i) = min(h);
   
    legend_list = [legend_list, "\bf{L} " + string(applied_load_list(i))];
    
    f_t = [f_t solution(i,1).friction{1}];
    f_b = [f_b solution(i,1).friction{2}];
    f_v = [f_v solution(i,1).friction{3}];
    f_coeff = [f_coeff solution(i,1).friction{4}];
end

figure(3); hold on
plot(applied_load_list,  hmin_list/(hmin_list(1)),".-","MarkerSize", 24)

figure(4);
subplot(2,2,1); hold on
set(gca, 'ColorOrder', newcolors);
plot(applied_load_list, f_t, "^--","LineWidth",1.6,"MarkerSize", 24)
plot(applied_load_list, f_b,"v--","LineWidth",1.6,"MarkerSize", 24)
subplot(2,2,2); hold on
plot(applied_load_list, f_coeff,".--","LineWidth",1.6,"MarkerSize", 24)
subplot(2,2,3); hold on
plot(applied_load_list, f_t + f_b, ".--","LineWidth",1.6,"MarkerSize", 24)
subplot(2,2,4); hold on
plot(applied_load_list, f_v, ".--","LineWidth",1.6,"MarkerSize", 24)


figure(1); 
xlabel("\boldmath{$x/b_{H}$}", "Interpreter","latex")
ylabel("\bf{p (Pa)}", "Interpreter","latex")
% title("\bf{Wi} = " + string(solution(1,1).wiessenberg_Number) + ", \bf{U} = " + string(solution(1,1).velocity) ,"Interpreter","latex")
title("\bf{De} = " + string(solution(1,1).deborah_Number) + ", \bf{U} = " + string(solution(1,1).velocity) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend(legend_list,"box","off","Interpreter","latex")

figure(2); 
xlabel("\boldmath{$x/b_{H}$}","Interpreter","latex")
xlim([-3 3])
ylabel("\bf{h (m)}","Interpreter","latex")
title("\bf{De} = " + string(solution(1,1).deborah_Number) + ", \bf{U} = " + string(solution(1,1).velocity) ,"Interpreter","latex")
% title("\bf{Wi} = " + string(solution(1,1).wiessenberg_Number) + ", \bf{U} = " + string(solution(1,1).velocity) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend(legend_list,"box","off","Interpreter","latex")

figure(3); 
xlabel("\bf{load}","Interpreter","latex")
ylabel("\boldmath{$h_{min} / h_{min, L = "+ string(min(applied_load_list)) +"}$}","Interpreter","latex")
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
ylabel("\boldmath{$\mu$} ","Interpreter","latex")
title("","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)

subplot(2,2,3); hold on
xlabel("\bf{L}", "Interpreter","latex")
ylabel(" \boldmath{$N/m$}", "Interpreter","latex")
title("\boldmath{$f_t + f_b$} ","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)

subplot(2,2,4); hold on
xlabel("\bf{L}","Interpreter","latex")
ylabel("\boldmath{$N/m$}","Interpreter","latex")
title("\boldmath{$f_{viscoelastic}$} ","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)


pressure_plot = figure(1);
gap_plot = figure(2);
minh_plot = figure(3);
friction_plot = figure(4);
set(pressure_plot,'Position',[10 10 980 720])
set(gap_plot,'Position',[10 10 980 720])
set(minh_plot,'Position',[35 35 1080 800])
set(friction_plot,'Position',[835 235 1600 1400])

% savepath = 'E:\Bilkent Dökümanları\Masterımsı\Haftalık\Dec 15 22\';

if ~isfolder(savepath + "/gap")
    mkdir(savepath + "/gap")
end
if ~isfolder(savepath + "/pressure")
    mkdir(savepath + "/pressure")
end

saveas(pressure_plot, savepath + "pressure/" + "EHLPressure_" +  "beta" + string(solution(1,1).viscocity_ratio) + "Wi" + string(solution(1,1).wiessenberg_Number) + "U" + string(solution(1,1).velocity) + ".svg")
saveas(gap_plot, savepath + "gap/" + "EHLGap_" +  "beta" + string(solution(1,1).viscocity_ratio) + "Wi" + string(solution(1,1).wiessenberg_Number) + "U" + string(solution(1,1).velocity) + ".svg")

saveas(pressure_plot, savepath + "pressure/" + "EHLPressure_" +  "beta" + string(solution(1,1).viscocity_ratio) + "Wi" + string(solution(1,1).wiessenberg_Number) + "U" + string(solution(1,1).velocity) + ".fig")
saveas(gap_plot, savepath + "gap/" + "EHLGap_" +  "beta" + string(solution(1,1).viscocity_ratio) + "Wi" + string(solution(1,1).wiessenberg_Number) + "U" + string(solution(1,1).velocity) + ".fig")

% saveas(minh_plot, savepath + "EHLminhU_" +  "beta" + string(beta) + "Wi" + string(Wi) + "L" + string(force) + ".svg")

close(pressure_plot)
close(gap_plot)


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