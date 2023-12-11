function EHLCompareDe_NoNDimensional(solution, savepath)
% solution contains the solution class EHL values for different loads but
% same velocity and wi. The length of the solution is equals the length of
% the results that will be plotted. 

savepath = string(savepath);
legend_list = [];
hmin_list = [];

pN_pNN_list = [];

f_t = [];
f_b = [];
f_v_b = [];
f_coeff = [];

de_list = [solution(1,1,:).deborah_Number];

newcolors = [0 0.4470 0.7410; 0 0.4470 0.7410;  ...
    0.8500 0.3250 0.0980; 0.8500 0.3250 0.0980; ...
    0.9290 0.6940 0.1250; 0.9290 0.6940 0.1250; ...
    0.4940 0.1840 0.5560; 0.4940 0.1840 0.5560; ...
    0.4660 0.6740 0.1880; 0.4660 0.6740 0.1880];

renkler = [43,140,190 ; 215,48,31; 35,139,69; 136,65,157; 204,76,2] / 256;

renkler_double = [43,140,190 ; 43,140,190 ; 215,48,31; 215,48,31; 35,139,69;  35,139,69; 136,65,157;   136,65,157 ] / 256;

for i = 1:length(solution(1,1,:))
%     if i ~= 2
    figure(1); hold on; set(gca, 'ColorOrder', renkler);
    plot([solution(1,1,i).domain.x] ./ [solution(1,1,i).domain.bH], [solution(1,1,i).pressure] ./ [max(real(solution(1,1,i).hertzian_pressure))], "LineWidth",2)
    figure(2); hold on ; set(gca, 'ColorOrder', renkler);
    plot([solution(1,1,i).domain.x] ./ [solution(1,1,i).domain.bH], [solution(1,1,i).h ./ solution(1,1,i).domain.href], "LineWidth",2)
    figure(5); hold on ; set(gca, 'ColorOrder', renkler);
    plot([solution(1,1,i).domain.x] ./ [solution(1,1,i).domain.bH], [solution(1,1,i).pressure0_old ./ max(real(solution(1,1,i).hertzian_pressure))], "LineWidth",1.7)
    figure(6); hold on ; set(gca, 'ColorOrder', renkler);
    plot([solution(1,1,i).domain.x] ./ [solution(1,1,i).domain.bH], [solution(1,1,i).pressure1_old * solution(1,1,i).deborah_Number ./ max(real(solution(1,1,i).hertzian_pressure))], "LineWidth",1.7)
    figure(7); hold on ; set(gca, 'ColorOrder', renkler);
    plot([solution(1,1,i).domain.x] ./ [solution(1,1,i).domain.bH], [OneDcentraldiff(solution(1,1,i).h, solution(1,1,i).domain.dx)],"--","LineWidth",1.7)
    figure(8); hold on ; set(gca, 'ColorOrder', renkler);
    plot([solution(1,1,i).domain.x] ./ [solution(1,1,i).domain.bH], [OneDcentraldiff(OneDcentraldiff(solution(1,1,i).h, solution(1,1,i).domain.dx),solution(1,1,i).domain.dx)], "--","LineWidth",1.7)
   
    pN_pNN_list = [pN_pNN_list  max(solution(1,1,i).pressure1_old) * solution(1,1,i).deborah_Number / max(solution(1,1,i).pressure)];
    
    h = solution(1,1,i).h;
    hmin_list(i) = min(h);
   
    legend_list = [legend_list, "\bf{De = } " + string(de_list(i))];
      
    f_t = [f_t solution(1,1,i).friction{1} * solution(1,1,i).domain.href / (solution(1,1,i).velocity * solution(1,1,i).domain.mu * solution(1,1,i).domain.L)];
    f_b = [f_b solution(1,1,i).friction{2} * solution(1,1,i).domain.href / (solution(1,1,i).velocity * solution(1,1,i).domain.mu * solution(1,1,i).domain.L)];
    f_v_b = [f_v_b solution(1,1,i).friction{5} * solution(1,1,i).domain.href / (solution(1,1,i).velocity * solution(1,1,i).domain.mu * solution(1,1,i).domain.L)];
    f_coeff = [f_coeff solution(1,1,i).friction{4}];
%     end
end
pN_pNN_list
% de_list = [de_list(1) de_list(3) de_list(4) de_list(5)];
de_list = [de_list(1) de_list(2) de_list(3) de_list(4)];
% hmin_list = [hmin_list(1) hmin_list(3) hmin_list(4) hmin_list(5)];
hmin_list = [hmin_list(1) hmin_list(2) hmin_list(3) hmin_list(4)];


figure(3); hold on; box on;
set(gca, 'ColorOrder', renkler);
plot(de_list,  hmin_list/(hmin_list(1)),".-","LineWidth",1.9,"MarkerSize", 24)
set(gca,'linewidth',1.1)

figure(9); hold on; box on;
set(gca, 'ColorOrder', renkler);
plot(de_list, pN_pNN_list,".-","LineWidth",1.9, "MarkerSize", 24)
set(gca,'linewidth',1.1)

figure(10); hold on; box on;
set(gca, 'ColorOrder', renkler);
plot(de_list, f_coeff,".-","LineWidth",1.9, "MarkerSize", 24)
set(gca,'linewidth',1.1)

% He = 
% figure(61); hold on; box on;
% plot()

figure(4);
subplot(2,2,1); hold on; box on; 
set(gca, 'ColorOrder', newcolors);
plot(de_list, f_t, "^--","LineWidth",1.9,"MarkerSize", 24)
plot(de_list, f_b,"v--","LineWidth",1.9,"MarkerSize", 24)
subplot(2,2,2); hold on
set(gca, 'ColorOrder', renkler);
plot(de_list, f_coeff,".--","LineWidth",1.9,"MarkerSize", 24)
subplot(2,2,3); hold on
set(gca, 'ColorOrder', renkler);
plot(de_list, f_b, ".--","LineWidth",1.9,"MarkerSize", 24)
subplot(2,2,4); hold on
set(gca, 'ColorOrder', renkler);
plot(de_list, f_v_b, ".--","LineWidth",1.9,"MarkerSize", 24)


figure(11); hold on; box on;
set(gca, 'ColorOrder', renkler_double);
plot(de_list, f_b, ".-","LineWidth",1.8,"MarkerSize", 22)
plot(de_list, f_b - f_v_b, ".--","LineWidth",1.8,"MarkerSize", 20)

figure(12); hold on; box on;
set(gca, 'ColorOrder', renkler);
plot(de_list, f_v_b, ".-","LineWidth",1.9,"MarkerSize", 24)

figure(13); hold on; box on;
set(gca, 'ColorOrder', renkler);
plot(de_list, f_b - f_v_b, ".-","LineWidth",1.9,"MarkerSize", 24)


figure(1); box on; 
xlabel("\boldmath{$x/b_{H}$}","Interpreter","latex")
ylabel("\boldmath{$p / p_{H}$}","Interpreter","latex")
title("\bf{U = }" + string(solution(1,1,1).velocity * solution(1,1,1).domain.mu / solution(1,1,1).domain.Rx / solution(1,1,1).domain.E)  + ...
    ", \bf{L = }" + string(solution(1,1,1).applied_load / solution(1,1,1).domain.Rx / solution(1,1,1).domain.E / solution(1,1,1).domain.L),"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend(legend_list,"box","off","Interpreter","latex")
set(gca,'linewidth',1.1)

figure(5); box on;
xlabel("\boldmath{$x/b_{H}$}","Interpreter","latex")
ylabel("","Interpreter","latex")
title("\bf{$p^N$ \ U = }" + string(solution(1,1,1).velocity * solution(1,1,1).domain.mu / solution(1,1,1).domain.Rx / solution(1,1,1).domain.E) + ...
    ", \bf{L = }" + string(solution(1,1,1).applied_load / solution(1,1,1).domain.Rx / solution(1,1,1).domain.E / solution(1,1,1).domain.L),"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend(legend_list,"box","off","Interpreter","latex")
set(gca,'linewidth',1.1)

figure(6); box on;
xlabel("\boldmath{$x/b_{H}$}","Interpreter","latex")
ylabel("\boldmath{$p / p_{H}$}","Interpreter","latex")
title("\bf{$p^{NN}$ \ U = }" + string(solution(1,1,1).velocity * solution(1,1,1).domain.mu / solution(1,1,1).domain.Rx / solution(1,1,1).domain.E) + ...
    ", \bf{L = }" + string(solution(1,1,1).applied_load / solution(1,1,1).domain.Rx / solution(1,1,1).domain.E / solution(1,1,1).domain.L),"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend(legend_list,"box","off","Interpreter","latex")
set(gca,'linewidth',1.1)

figure(7); box on; 
xlabel("\boldmath{$x/b_{H}$}","Interpreter","latex")
ylabel("\boldmath{}", "Interpreter","latex")
title("\bf{$\frac{\partial h}{\partial x}$ \ U = }" + string(solution(1,1,1).velocity * solution(1,1,1).domain.mu / solution(1,1,1).domain.Rx / solution(1,1,1).domain.E) + ...
    ", \bf{L = }" + string(solution(1,1,1).applied_load / solution(1,1,1).domain.Rx / solution(1,1,1).domain.E / solution(1,1,1).domain.L),"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend(legend_list,"box","off","Interpreter","latex")
xlim([-3 3])
set(gca,'linewidth',1.1)

figure(8); box on;
xlabel("\boldmath{$x/b_{H}$}","Interpreter","latex")
ylabel("boldmath{$\frac{1}{m}$}", "Interpreter","latex")
ylabel("\bf{p (Pa)}","Interpreter","latex")
title("\bf{$\frac{\partial^2 h}{\partial x^2}$ \ U = }" + string(solution(1,1,1).velocity * solution(1,1,1).domain.mu / solution(1,1,1).domain.Rx / solution(1,1,1).domain.E) + ...
    ", \bf{L = }" + string(solution(1,1,1).applied_load / solution(1,1,1).domain.Rx / solution(1,1,1).domain.E / solution(1,1,1).domain.L),"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend(legend_list,"box","off","Interpreter","latex")
xlim([-3 3])
set(gca,'linewidth',1.1)


figure(2); box on;
xlabel("\boldmath{$x/b_{H}$}","Interpreter","latex")
xlim([-3 3])
ylabel("\boldmath{$h / h_{ref}$}","Interpreter","latex")
title("\bf{U} = " + string(solution(1,1,1).velocity * solution(1,1,1).domain.mu / solution(1,1,1).domain.Rx / solution(1,1,1).domain.E) + ...
    ", \bf{L} = " + string(solution(1,1,1).applied_load / solution(1,1,1).domain.Rx / solution(1,1,1).domain.E / solution(1,1,1).domain.L) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend(legend_list,"box","off","Interpreter","latex")
set(gca,'linewidth',1.1)

figure(3); box on;
xlabel("\bf{De}","Interpreter","latex")
ylabel("\boldmath{$h_{min} / h_{min," +  string(de_list(1))  +"}$}","Interpreter","latex")
title("","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
set(gca,'linewidth',1.1)

figure(9); box on;
xlabel("\bf{De}","Interpreter","latex")
ylabel("\boldmath{$ max(p^{NN}) / max(p^{T}) $}","Interpreter","latex")
title("","Interpreter","latex")
set(gca,'fontsize',16)
set(gca,'TickLabelInterpreter','latex')
set(gca,'linewidth',1.1)

figure(10); box on;
xlabel("\bf{De}","Interpreter","latex")
ylabel("\boldmath{$f$}","Interpreter","latex")
title("","Interpreter","latex")
set(gca,'fontsize',16)
set(gca,'TickLabelInterpreter','latex')
set(gca,'linewidth',1.1)

% legend(legend_list,"box","off")

figure(4); 
subplot(2,2,1); hold on; box on;
ylabel("\boldmath{$friction$} ","Interpreter","latex")
title("\boldmath{$f_t$ },   \bf{U} = " + string(solution(1,1,1).velocity * solution(1,1,1).domain.mu / solution(1,1,1).domain.Rx / solution(1,1,1).domain.E),"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
set(gca,'linewidth',1.1)
subplot(2,2,2); hold on; box on;
ylabel("\boldmath{$\mu$} ","Interpreter","latex")
title("","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
set(gca,'linewidth',1.1)
subplot(2,2,3); hold on; box on;
xlabel("\bf{De}", "Interpreter","latex")
ylabel(" \boldmath{$friction$}", "Interpreter","latex")
title("\boldmath{$f_t + f_b$} ","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
set(gca,'linewidth',1.1)
subplot(2,2,4); hold on; box on;
xlabel("\bf{De}","Interpreter","latex")
ylabel("\boldmath{$friction$}","Interpreter","latex")
title("\boldmath{$f_{viscoelastic}$} ","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
set(gca,'linewidth',1.1)

figure(11);
xlabel("\bf{De}", "Interpreter","latex")
ylabel(" \boldmath{$friction$}", "Interpreter","latex")
title("\boldmath{$f_b$} ","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
set(gca,'linewidth',1.1)

figure(12);
xlabel("\bf{De}","Interpreter","latex")
ylabel("\boldmath{$friction$}","Interpreter","latex")
title("\boldmath{$f_{NN}$} ","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
set(gca,'linewidth',1.1)

figure(13);
xlabel("\bf{De}","Interpreter","latex")
ylabel("\boldmath{$friction$}","Interpreter","latex")
title("\boldmath{$f_{N}$} ","Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
set(gca,'linewidth',1.1)


pressure_plot = figure(1);
gap_plot = figure(2);
minh_plot = figure(3);
friction_plot = figure(4);
pressure0_plot = figure(5);
pressure1_plot = figure(6);
dgap_plot = figure(7);
d2gap_plot = figure(8);
pNNpT_plot = figure(9);
f_plot = figure(10);
f_b = figure(11);
f_b_v = figure(12);
f_b_n = figure(13);

set(pressure_plot,'Position',[10 10 980 720])
set(gap_plot,'Position',[10 10 980 720])
set(minh_plot,'Position',[35 35 1080 800])
set(friction_plot,'Position',[835 235 1600 1400])
set(pressure0_plot,'Position',[10 10 980 720])
set(pressure1_plot,'Position',[10 10 980 720])
set(dgap_plot,'Position',[10 10 980 720])
set(d2gap_plot,'Position',[10 10 980 720])
set(pNNpT_plot,'Position',[35 35 1080 800])
set(f_plot,'Position',[10 90 850 720])
% savepath = 'E:\Bilkent Dökümanları\Masterımsı\Haftalık\Dec 15 22\';
set(f_b, 'Position', [10 90 850 720])
set(f_b_v, 'Position', [10 90 850 720])
set(f_b_n, 'Position', [10 90 850 720])

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

% if (pN_pNN_list(1)) == Inf
%     disp([solution.applied_load solution.velocity])
%     pN_pNN_list(1) = 100;
% end
%     
% pN_pNN_list    

% U_nondimensional =
% 
%    1.8750e-11
%                   for U = 0.0001, mu = 0.1, E1 = 1e7,   v = 0.5, Rx = 0.02

% U = mu  * u / E / R

% force_nondmensional =  F / Rx / 20 / bH / E 
% 
% ans =
% 
%     0.0136

% W = W / E / R / L

end