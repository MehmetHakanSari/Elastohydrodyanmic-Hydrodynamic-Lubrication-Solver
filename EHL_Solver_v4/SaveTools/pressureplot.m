function pressureplot(solution, savepath)

file_name = "L";

file_name = file_name + string(solution.applied_load);
file_name = file_name + "_Wi" + string(solution.wiessenberg_Number);
file_name = file_name + "_beta" + string(solution.viscocity_ratio);
file_name = file_name + "_U" + string(solution.velocity);
file_name = file_name + "_bH" + string(solution.domain.domain_coeff);


if ~isfolder("Results/Figures/" + file_name)
       mkdir("Results/Figures/" + file_name);
end

if ~isfolder("Results/Figures/" + file_name + "/pressure")
       mkdir("Results/Figures/" + file_name + "/pressure");
end

savepath = savepath + file_name + "\pressure\";

solution.domain.epsilon * solution.wiessenberg_Number;

figure(1);
plot(solution.domain.x, solution.pressure, "-", "LineWidth", 1.7);
figure(2); subplot(1, 2, 1); hold on
plot(solution.domain.x, solution.pressure0_old, "--", "LineWidth", 1.7);
plot(solution.domain.x, solution.pressure1_old  * solution.domain.epsilon * solution.wiessenberg_Number, ":", "LineWidth", 1.9);
figure(2); subplot(1, 2, 2); hold on
plot(solution.domain.x, solution.pressure, "-", "LineWidth", 1.7);
plot(solution.domain.x, solution.pressure0_old, "--", "LineWidth", 1.7);
plot(solution.domain.x, solution.pressure1_old  * solution.domain.epsilon * solution.wiessenberg_Number, ":", "LineWidth", 1.9);
figure(3); hold on
plot(solution.domain.x, solution.pressure0_old, "--", "LineWidth", 1.7);
figure(4); hold on
plot(solution.domain.x, solution.pressure1_old  * solution.domain.epsilon * solution.wiessenberg_Number, ":", "LineWidth", 1.9);


figure(1); 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\bf{p (Pa)}", "Interpreter","latex")
title("\bf{Wi} = " + string(solution(1,1).wiessenberg_Number) + ", \bf{U} = " + string(solution(1,1).velocity) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
xlim([min(solution.domain.x) max(solution.domain.x)])
legend("\boldmath{$p_{T}$}", "box","off","Interpreter","latex")
% legend(legend_list,"box","off","Interpreter","latex")
% if velocity == 0.1 
%     xlim([-20 20])
% elseif velocity == 1
%     xlim([-20 20])
% end

figure(2); 
subplot(1, 2, 1)
xlabel("\boldmath{$x$}","Interpreter","latex")
ylabel("\bf{p (Pa)}","Interpreter","latex")
title("\bf{Wi} = " + string(solution(1,1).wiessenberg_Number) + ", \bf{U} = " + string(solution(1,1).velocity) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend("\boldmath{$p_{N}$}", "\boldmath{$p_{NN}$}", "box","off","Interpreter","latex")
xlim([min(solution.domain.x) max(solution.domain.x)])

figure(2); 
subplot(1, 2, 2)
xlabel("\boldmath{$x$}","Interpreter","latex")
ylabel("\bf{p (Pa)}","Interpreter","latex")
title("\bf{Wi} = " + string(solution(1,1).wiessenberg_Number) + ", \bf{U} = " + string(solution(1,1).velocity) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
legend("\boldmath{$p_{T}$}", "\boldmath{$p_{N}$}", "\boldmath{$p_{NN}$}", "box","off","Interpreter","latex")
xlim([min(solution.domain.x) max(solution.domain.x)])

figure(3); 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\bf{p (Pa)}", "Interpreter","latex")
title("\bf{Wi} = " + string(solution(1,1).wiessenberg_Number) + ", \bf{U} = " + string(solution(1,1).velocity) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
xlim([min(solution.domain.x) max(solution.domain.x)])
legend("\boldmath{$p_{N}$}", "box","off","Interpreter","latex")

figure(4); 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\bf{p (Pa)}", "Interpreter","latex")
title("\bf{Wi} = " + string(solution(1,1).wiessenberg_Number) + ", \bf{U} = " + string(solution(1,1).velocity) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
xlim([min(solution.domain.x) max(solution.domain.x)])
legend("\boldmath{$p_{NN}$}", "box","off","Interpreter","latex")

pT_plot = figure(1);
pAll_plot = figure(2);
pN_plot = figure(3);
pNN_plot = figure(4);

set(pT_plot,'Position',[280 200 980 720])
set(pAll_plot,'Position',[80 50 1680 870])
set(pN_plot,'Position',[280 200 980 720])
set(pNN_plot,'Position',[80 50 1680 870])

saveas(pT_plot, savepath + "pT.svg")
saveas(pAll_plot, savepath + "pNN.svg")
saveas(pN_plot, savepath + "pT.svg")
saveas(pNN_plot, savepath + "pNN.svg")

saveas(pT_plot, savepath + "pT.fig")
saveas(pAll_plot, savepath + "pNN.fig")
saveas(pN_plot, savepath + "pT.fig")
saveas(pNN_plot, savepath + "pNN.fig")


end