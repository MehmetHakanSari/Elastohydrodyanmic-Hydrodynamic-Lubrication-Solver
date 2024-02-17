function heightplot(solution, savepath)

file_name = "L";

file_name = file_name + string(solution.applied_load);
file_name = file_name + "_Wi" + string(solution.wiessenberg_Number);
file_name = file_name + "_beta" + string(solution.viscocity_ratio);
file_name = file_name + "_U" + string(solution.velocity);
file_name = file_name + "_bH" + string(solution.domain.domain_coeff);


if ~isfolder("Results/Figures/" + file_name)
       mkdir("Results/Figures/" + file_name);
end

if ~isfolder("Results/Figures/" + file_name + "/height")
       mkdir("Results/Figures/" + file_name + "/height");
end

savepath = savepath + file_name + "\height\";
 
figure(1);
plot(solution.domain.x, solution.h, "--", "LineWidth", 1.7);
figure(2); subplot(1, 2, 1); hold on
plot(solution.domain.x, OneDcentraldiff(solution.h, solution.domain.dx), "--", "LineWidth", 1.7);
figure(2); subplot(1, 2, 2); hold on
plot(solution.domain.x, OneDcentraldiff(OneDcentraldiff(solution.h, solution.domain.dx), solution.domain.dx), "--", "LineWidth", 1.7);
figure(3); hold on
plot(solution.domain.x, OneDcentraldiff(solution.h, solution.domain.dx), "--", "LineWidth", 1.7);
figure(4);  hold on
plot(solution.domain.x, OneDcentraldiff(OneDcentraldiff(solution.h, solution.domain.dx), solution.domain.dx), "--", "LineWidth", 1.7);


figure(1); 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\bf{h (m)}", "Interpreter","latex")
title("\bf{h, \ \ \ Wi} = " + string(solution(1,1).wiessenberg_Number) + ", \bf{U} = " + string(solution(1,1).velocity) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])


figure(2); subplot(1, 2, 1) 
xlabel("\boldmath{$x$}","Interpreter","latex")
ylabel("\boldmath{}", "Interpreter","latex")
title("\boldmath{$\frac{\partial h}{\partial x}$, Wi} = " + string(solution(1,1).wiessenberg_Number) + ", \bf{U} = " + string(solution(1,1).velocity) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
% legend("\bf{p_N}", "\bf{p_NN}", "box","off","Interpreter","latex")

figure(2);  subplot(1, 2, 2)
xlabel("\boldmath{$x$}","Interpreter","latex")
ylabel("boldmath{$\frac{1}{m}$}", "Interpreter","latex")
title("\boldmath{$\frac{\partial^2 h}{\partial x^2}$, Wi} = " + string(solution(1,1).wiessenberg_Number) + ", \bf{U} = " + string(solution(1,1).velocity) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
xlim([min(solution.domain.x) max(solution.domain.x)])
% legend("\bf{p_N}", "\bf{p_NN}", "box","off","Interpreter","latex")
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])

figure(3); 
xlabel("\boldmath{$x$}","Interpreter","latex")
ylabel("\boldmath{}", "Interpreter","latex")
title("\boldmath{$\frac{\partial h}{\partial x}$, Wi} = " + string(solution(1,1).wiessenberg_Number) + ", \bf{U} = " + string(solution(1,1).velocity) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
% legend("\bf{p_N}", "\bf{p_NN}", "box","off","Interpreter","latex")

figure(4);  
xlabel("\boldmath{$x$}","Interpreter","latex")
ylabel("boldmath{$\frac{1}{m}$}", "Interpreter","latex")
title("\boldmath{$\frac{\partial^2 h}{\partial x^2}$, Wi} = " + string(solution(1,1).wiessenberg_Number) + ", \bf{U} = " + string(solution(1,1).velocity) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
xlim([min(solution.domain.x) max(solution.domain.x)])
% legend("\bf{p_N}", "\bf{p_NN}", "box","off","Interpreter","latex")
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])

h_plot = figure(1);
dhd2h_plot = figure(2);
dh_plot = figure(3);
d2h_plot = figure(4); 

set(h_plot,'Position',[280 200 980 720])
set(dhd2h_plot,'Position',[80 50 1680 870])
set(h_plot,'Position',[280 200 980 720])
set(d2h_plot,'Position',[80 50 1680 870])

saveas(h_plot, savepath + "h.svg")
saveas(dh_plot, savepath + "h_derivatives.svg")
saveas(h_plot, savepath + "h_firstderivative.svg")
saveas(dh_plot, savepath + "h_secondderivative.svg")

saveas(h_plot, savepath + "h.fig")
saveas(dhd2h_plot, savepath + "h_derivatives.fig")
saveas(dh_plot, savepath + "h_firstderivative.fig")
saveas(d2h_plot, savepath + "h_secondderivative.fig")


end