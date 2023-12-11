function EHLComparebeta_NoNDimensional(solution, savepath)
% solution contains the solution class EHL values for different loads but
% same velocity and wi. The length of the solution is equals the length of
% the results that will be plotted. 


savepath = string(savepath);
legend_list = [];
hmin_list = [];

applied_beta_list = [solution(1,:).viscocity_ratio];

% size(solution)

f_t = [];
f_b = [];
f_v = [];
f_coeff = [];

newcolors = [0 0.4470 0.7410; 0 0.4470 0.7410;  ...
    0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; ...
    0.9290 0.6940 0.1250; 0.9290 0.6940 0.1250; ...
    0.4940 0.1840 0.5560; 0.4940 0.1840 0.5560; ...
    0.4660 0.6740 0.1880; 0.4660 0.6740 0.1880];

renkler = [43,140,190 ; 215,48,31; 35,139,69; 136,65,157; 204,76,2] / 256;

for i = 1:length(solution(1,:))
    figure(1); hold on; set(gca, 'ColorOrder', renkler);
    plot([solution(1,i).domain.x] ./ [solution(1,i).domain.bH], [solution(1,i).pressure ./ max(real(solution(1,i).hertzian_pressure))], "LineWidth",2)
    figure(2); hold on; set(gca, 'ColorOrder', renkler);
    plot([solution(1,i).domain.x] ./ [solution(1,i).domain.bH], [solution(1,i).h ./ solution(1,i).domain.href], "LineWidth",2)
    
    h = solution(1,i).h;
    hmin_list(i) = min(h);
   
    legend_list = [legend_list, "\bf{$\beta$} =" + string(applied_beta_list(i))];
    
    
    f_t = [f_t solution(1,i).friction{1} * solution(1,i).domain.href / (solution(1,i).velocity * solution(1,i).domain.mu * solution(1,i).domain.L)];
    f_b = [f_b solution(1,i).friction{2} * solution(1,i).domain.href / (solution(1,i).velocity * solution(1,i).domain.mu * solution(1,i).domain.L)];
    f_v = [f_v solution(1,i).friction{3} * solution(1,i).domain.href / (solution(1,i).velocity * solution(1,i).domain.mu * solution(1,i).domain.L)];
    f_coeff = [f_coeff solution(1,i).friction{4}];
end

figure(3); hold on; box on;
set(gca, 'ColorOrder', renkler);
plot(applied_beta_list,  hmin_list/(hmin_list(1)),".-","LineWidth",1.9,"MarkerSize", 24)
set(gca,'linewidth',1.1)

figure(4);
subplot(2,2,1); hold on; box on;
set(gca, 'ColorOrder', newcolors);
plot(applied_beta_list, f_t, "^--","LineWidth",1.9,"MarkerSize", 24)
plot(applied_beta_list, f_b,"v--","LineWidth",1.9,"MarkerSize", 24)
subplot(2,2,2); hold on
set(gca, 'ColorOrder', renkler);
plot(applied_beta_list, f_coeff,".--","LineWidth",1.9,"MarkerSize", 24)
subplot(2,2,3); hold on
set(gca, 'ColorOrder', renkler);
plot(applied_beta_list, f_t + f_b, ".--","LineWidth",1.9,"MarkerSize", 24)
subplot(2,2,4); hold on
set(gca, 'ColorOrder', renkler);
plot(applied_beta_list, f_v, ".--","LineWidth",1.9,"MarkerSize", 24)

figure(5); hold on; box on;
set(gca, 'ColorOrder', renkler);
plot(applied_beta_list, f_coeff,".-","LineWidth", 2, "MarkerSize", 24)
set(gca,'linewidth',1.1)


figure(1); box on;
xlabel("\boldmath{$x/b_{H}$}", "Interpreter","latex")
ylabel("\boldmath{$p / p_{H}$}", "Interpreter","latex")
title("\bf{Wi} = " + string(solution(1,1).wiessenberg_Number) + ...
    ", \bf{U} = " + string(solution(1,1).velocity * solution(1,1).domain.mu / solution(1,1).domain.Rx / solution(1,1).domain.E), "Interpreter","latex")
% title("\bf{De} = " + string(solution(1,1).deborah_Number) + ", \bf{U} = " + string(solution(1,1).velocity) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend(legend_list,"box","off","Interpreter","latex")
% if velocity == 0.1 
%     xlim([-20 20])
% elseif velocity == 1
%     xlim([-20 20])
% end
set(gca,'linewidth',1.1)

figure(2); box on;
xlabel("\boldmath{$x/b_{H}$}","Interpreter","latex")
xlim([-3 3])
ylabel("\boldmath{$h / h_{ref}$}","Interpreter","latex")
title("\bf{Wi} = " + string(solution(1,1).wiessenberg_Number) + ...
    ", \bf{\Gamma} = " + string(solution(1,1).velocity * solution(1,1).domain.mu / solution(1,1).domain.Rx / solution(1,1).domain.E), "Interpreter","latex")
% title("\bf{De} = " + string(solution(1,1).deborah_Number) + ", \bf{U} = " + string(solution(1,1).velocity) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend(legend_list,"box","off","Interpreter","latex")
set(gca,'linewidth',1.1)

figure(3); box on;
xlabel("\boldmath{$\beta$}","Interpreter","latex")
ylabel("\boldmath{$h_{min} / h_{min, " +  string(applied_beta_list(1))  +"}$}","Interpreter","latex")
title("","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
set(gca,'linewidth',1.1)

figure(5); box on;
xlabel("\bf{$\beta$}", "Interpreter","latex")
ylabel("\bf{$f$}  ","Interpreter","latex")
title(" ","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
set(gca,'linewidth',1.1)

figure(4); 
subplot(2,2,1); hold on; box on;
% xlabel("\bf{Wi}","Interpreter","latex")
% ylabel("\boldmath{$f_t / f_{t,L=" + string(force_list(1)) + "}$}","Interpreter","latex")
ylabel("\boldmath{$friction$} ","Interpreter","latex")
title("\boldmath{$f_t \ and \ f_b$},   \bf{U} = " + string(solution(1,1).velocity * solution(1,1).domain.mu / solution(1,1).domain.Rx / solution(1,1).domain.E),"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
set(gca,'linewidth',1.1)

subplot(2,2,2); hold on; box on;
% xlabel("\bf{L}","Interpreter","latex")
% ylabel("\boldmath{$f_b / f_{b,L=" + string(force_list(1)) + "}$}","Interpreter","latex")
ylabel("\boldmath{$\mu$}  ","Interpreter","latex")
title(" ","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
set(gca,'linewidth',1.1)

subplot(2,2,3); hold on; box on;
xlabel("\bf{F}", "Interpreter","latex")
ylabel(" \boldmath{$friction$}", "Interpreter","latex")
title("\boldmath{$f_t + f_b$} ","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
set(gca,'linewidth',1.1)

subplot(2,2,4); hold on; box on;
xlabel("\bf{F}","Interpreter","latex")
ylabel("\boldmath{$friction$}","Interpreter","latex")
title("\boldmath{$f_{viscoelastic}$} ","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
set(gca,'linewidth',1.1)

pressure_plot = figure(1);
gap_plot = figure(2);
minh_plot = figure(3);
friction_plot = figure(4);
f_coeff = figure(5);
set(pressure_plot,'Position',[10 10 850 720])
set(gap_plot,'Position',[10 10 850 720])
set(minh_plot,'Position',[35 35 850 720])
set(friction_plot,'Position',[835 235 850 720])
set(f_coeff,'Position',[10 10 850 720])

% savepath = 'E:\Bilkent Dökümanları\Masterımsı\Haftalık\Dec 15 22\';

if ~isfolder(savepath + "/gap")
    mkdir(savepath + "/gap")
end
if ~isfolder(savepath + "/pressure")
    mkdir(savepath + "/pressure")
end

saveas(pressure_plot, savepath + "pressure/" + "EHLPressureNoNWi_" +  "W" + string(solution(1,1).applied_load) + "Wi" + string(solution(1,1).wiessenberg_Number) + "U" + string(solution(1,1).velocity) + ".svg")
saveas(gap_plot, savepath + "gap/" + "EHLGapNoNWi_" +  "W" + string(solution(1,1).applied_load) + "Wi" + string(solution(1,1).wiessenberg_Number) + "U" + string(solution(1,1).velocity) + ".svg")

saveas(pressure_plot, savepath + "pressure/" + "EHLPressureNoNWi_" +  "W" + string(solution(1,1).applied_load) + "Wi" + string(solution(1,1).wiessenberg_Number) + "U" + string(solution(1,1).velocity) + ".fig")
saveas(gap_plot, savepath + "gap/" + "EHLGapNoNWi_" +  "W" + string(solution(1,1).applied_load) + "Wi" + string(solution(1,1).wiessenberg_Number) + "U" + string(solution(1,1).velocity) + ".fig")

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