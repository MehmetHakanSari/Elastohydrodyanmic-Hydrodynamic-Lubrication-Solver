function stressfieldplot(solution, savepath)

file_name = "L";

file_name = file_name + string(solution.applied_load);
file_name = file_name + "_Wi" + string(solution.wiessenberg_Number);
file_name = file_name + "_beta" + string(solution.viscocity_ratio);
file_name = file_name + "_U" + string(solution.velocity);
file_name = file_name + "_bH" + string(solution.domain.domain_coeff);


if ~isfolder("Results/Figures/" + file_name)
       mkdir("Results/Figures/" + file_name);
end

if ~isfolder("Results/Figures/" + file_name + "/stressfield")
       mkdir("Results/Figures/" + file_name + "/stressfield");
end

savepath = savepath + file_name + "\stressfield\";


[J, N] = size(solution.stress_field{1,1});
h = solution.h;

%GRIDING
X = zeros(J, N);
Y = zeros(J, N);

for i = 1:J
    X(i,:) = solution.domain.x;
end
dy = solution.h / (J - 1);
% dy_nondimensional = solution.h / (J - 1) / solution.domain.href;
for i = 1:J
    Y(i,:) = solution.h - dy .* (i-1);
%     Y(i,:) = solution.h / solution.domain.h_ref - dy_nondimensional .* (i-1);
end

dydx = -Y ./ (h / solution.domain.href^2) .* OneDcentraldiff(h / solution.domain.href,  solution.domain.dx / solution.domain.L);

% solution.domain.epsilon * solution.wiessenberg_Number * (1 - solution.viscosity_ratio)
figure(1); subplot(2,2,1)
surf(X, Y, solution.stress_field{1,1}, "LineStyle", "None");
figure(1); hold on; subplot(2,2,3)
surf(X, Y, solution.stress_field{1,2}, "LineStyle", "None");
% % figure(1); hold on; subplot(3,2,3)
surf(X, Y, solution.stress_field{1,3}, "LineStyle", "None");

figure(1);  subplot(2,2,2)
surf(X, Y, solution.stress_field{2,1} * solution.domain.epsilon * solution.wiessenberg_Number, "LineStyle", "None");
figure(1); hold on;  subplot(2,2,4)
surf(X, Y, solution.stress_field{2,2} * solution.domain.epsilon * solution.wiessenberg_Number, "LineStyle", "None");
% figure(1); hold on;  subplot(3,2,6)
% surf(X, Y, solution.stress_field{2,3} * solution.domain.epsilon * solution.wiessenberg_Number, "LineStyle", "None");

gradTDxx = TwoDcentraldiff(solution.stress_field{2,1} , solution.domain.dx / solution.domain.L, solution.h / (J - 1) / solution.domain.href);
gradTDxy = TwoDcentraldiff(solution.stress_field{2,2} , solution.domain.dx / solution.domain.L, solution.h / (J - 1) / solution.domain.href);
gradTxx = TwoDcentraldiff(solution.stress_field{1,1} , solution.domain.dx / solution.domain.L, solution.h / (J - 1) / solution.domain.href);
gradTxy = TwoDcentraldiff(solution.stress_field{1,2} , solution.domain.dx / solution.domain.L, solution.h / (J - 1) / solution.domain.href);
gradu = TwoDcentraldiff(solution.velocity_field{1,1} , solution.domain.dx / solution.domain.L, solution.h / (J - 1) / solution.domain.href);
graduD = TwoDcentraldiff(solution.velocity_field{2,1} , solution.domain.dx / solution.domain.L, solution.h / (J - 1) / solution.domain.href);
grad2u = TwoDcentraldiff(gradu{2} , solution.domain.dx / solution.domain.L, solution.h / (J - 1) / solution.domain.href);
grad2uD = TwoDcentraldiff(graduD{2} , solution.domain.dx / solution.domain.L, solution.h / (J - 1) / solution.domain.href);

DTDxxDx = gradTDxx{1} + gradTDxx{2} .* dydx;
DTDxyDx = gradTDxy{1} + gradTDxy{2} .* dydx;
DTxxDx = gradTxx{1} + gradTxx{2} .* dydx;

% figure(7); hold on
% surf(X, Y, DTDxxDx * solution.domain.epsilon * solution.wiessenberg_Number, "LineStyle", "None");
% figure(8); hold on
% surf(X, Y, gradTDxy{2} * solution.domain.epsilon * solution.wiessenberg_Number, "LineStyle", "None");

% figure(9); hold on
% surf(X, Y, (DTDxxDx + gradTDxy{2}) * solution.domain.epsilon * solution.wiessenberg_Number, "LineStyle", "None");

dpdx = (solution.viscocity_ratio * grad2u{2} + DTxxDx + gradTxy{2});
dpDdx = (solution.viscocity_ratio * grad2uD{2} + DTDxxDx + gradTDxy{2}) * solution.deborah_Number;
figure(10); hold on;
surf(X, Y,  dpdx + dpDdx,"LineStyle", "None");
% 
% figure(61);
% contour(X,Y, dpdx + dpDdx)

% 
% figure(11); hold on;
% plot(solution.domain.x, solution.domain.dx * cumtrapz(dy .* trapz(dpdx + dpDdx)), "LineWidth" , 2)

hm = trapz(solution.h.^-2) / trapz(solution.h.^-3); 
F1 = 6 * OneDcentraldiff(solution.h .* solution.thetaD_old, 1);
F2 = 2 * (1 - solution.viscocity_ratio) * OneDcentraldiff(solution.theta_old .* OneDcentraldiff(solution.h, 1) .* (1 - 3 * hm ./ solution.h) .* ...
    (2 - 3 * hm ./ solution.h), 1);

figure(11), hold on;
% plot(solution.domain.x, F1, "Linewidth", 2)
plot(solution.domain.x / solution.domain.bH, F2, "Linewidth", 2)


%

uT = solution.velocity_field{1,1} + solution.velocity_field{2,1} * solution.deborah_Number;
vT = solution.velocity_field{1,2} + solution.velocity_field{2,2} * solution.deborah_Number;
DTTxxDx = DTxxDx + DTDxxDx * solution.deborah_Number;
dTTxxdy = gradTxx{2} + gradTDxx{2} * solution.deborah_Number;
convective_term_xx_D = solution.velocity_field{2,1} .* DTDxxDx + solution.velocity_field{2,2} .* gradTDxx{2}; 
convective_term_xy_D = solution.velocity_field{2,1} .* DTDxyDx + solution.velocity_field{2,2} .* gradTDxy{2}; 
convective_term_T = uT .* DTTxxDx + vT .* dTTxxdy; 
row_len = length(convective_term_T(:,1));
colors = ["#0c2c84", "#990000"];
% colors = ["#990000", "#0c2c84"];
ls = ":";
fs = 24;
figure(12); 
subplot(3,1,1); hold on; box on;
txxline1 = plot(solution.domain.x / solution.domain.bH, convective_term_xx_D(round(row_len/3),:)/ max(abs(convective_term_xx_D),[],"all"),ls,"LineWidth", 1.8, "Color",colors(1));
xlim([-3 3])
set(gca,'linewidth',1.2)
ylabel("\boldmath{$y=2h/3$}","Interpreter", "latex")
subplot(3,1,3); hold on; box on;
txxline2 = plot(solution.domain.x / solution.domain.bH, convective_term_xx_D(end,:)/max(abs(convective_term_xx_D),[],"all"),ls,"LineWidth", 1.8,"Color",colors(1));
xlim([-3 3])
set(gca,'linewidth',1.2)
ylabel("\boldmath{$y=0$}","Interpreter", "latex")
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
subplot(3,1,2); hold on; box on;
txxline3 = plot(solution.domain.x / solution.domain.bH, convective_term_xx_D(round(2*row_len/3),:) / max(abs(convective_term_xx_D),[],"all"),ls,"LineWidth", 1.8,"Color",colors(1));  
xlim([-3 3])
set(gca,'linewidth',1.2)
ylabel("\boldmath{$y=h/3$}","Interpreter", "latex")
% legend("W = 800N, De = 0.05", "W = 800N, De = 0", "W = 100N, De = 0", "W = 100N, De = 0.05", "box", "off", "Interpreter", "latex")
% legend("W = 800N, De = 0.05", "W = 800N, De = 0", "W = 100N, De = 0", "W = 100N, De = 0.05", "box", "off", "Interpreter", "latex")
hh = figure(12);
set(findall([hh], '-property', 'Interpreter'), 'Interpreter', 'latex')
set(findall([hh], '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex', 'FontSize', fs-6)
set(findobj([hh], 'Type', 'legend'), 'FontSize', fs-4)
% set(findall([txxline1 txxline2 txxline3], '-property', 'Interpreter'), 'Interpreter', 'latex')
% set(findall([txxline1 txxline2 txxline3], '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex', 'FontSize', fs-6)
% set(findobj([txxline1 txxline2 txxline3], 'Type', 'legend'), 'FontSize', fs-4)

% figure(13); 
% subplot(3,1,1); hold on; box on;
% plot(solution.domain.x / solution.domain.bH, convective_term_xy_D(round(row_len/3),:)/ max(convective_term_xy_D,[],"all"),"-","LineWidth", 1.8, "Color",colors(1))
% xlim([-5 5])
% subplot(3,1,2); hold on; box on;
% plot(solution.domain.x / solution.domain.bH, convective_term_xy_D(round(2*row_len/3),:) / max(convective_term_xy_D,[],"all"),"-","LineWidth", 1.8,"Color",colors(1))
% xlim([-5 5])
% subplot(3,1,3); hold on; box on;
% plot(solution.domain.x / solution.domain.bH, convective_term_xy_D(end,:)/ max(convective_term_xy_D,[],"all"),"-","LineWidth", 1.8, "Color",colors(1))
% xlim([-5 5])

figure(1);  subplot(2,2,1)
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\boldmath{$y$}", "Interpreter","latex")
title("\boldmath{$T_{xx}^{N}$}", "Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
colorbar
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
% xlim([-6 * (solution.domain.L) / (solution.domain.domain_coeff), -1 * (solution.domain.L) / (solution.domain.domain_coeff)])
ylim([0 h(251) * 2.5])
% set(gca, "Yscale", "log")

figure(1);  subplot(2,2,3)
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\boldmath{$y$}", "Interpreter","latex")
title("\boldmath{$T_{xy}^{N}$}", "Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
colorbar
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
% xlim([-6 * (solution.domain.L) / (solution.domain.domain_coeff), -1 * (solution.domain.L) / (solution.domain.domain_coeff)])
ylim([0 h(251) * 2.5])
% set(gca, "Yscale", "log")

% figure(1);  subplot(3,2,3);
% xlabel("\boldmath{$x$}", "Interpreter","latex")
% ylabel("\boldmath{$y$}", "Interpreter","latex")
% title("\boldmath{$T_{yy}^{N}$}", "Interpreter","latex")
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'fontsize',24)
% colorbar
% xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
% % xlim([-6 * (solution.domain.L) / (solution.domain.domain_coeff), -1 * (solution.domain.L) / (solution.domain.domain_coeff)])
% ylim([0 h(251) * 2.5])
% % set(gca, "Yscale", "log")

figure(1);  subplot(2,2,2) 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\boldmath{$y$}", "Interpreter","latex")
title("\boldmath{$T_{xx}^{NN}$}", "Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
colorbar
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
% xlim([-6 * (solution.domain.L) / (solution.domain.domain_coeff), -1 * (solution.domain.L) / (solution.domain.domain_coeff)])
ylim([0 h(251) * 2.5])
% set(gca, "Yscale", "log")

figure(1);  subplot(2,2,4) 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\boldmath{$y$}", "Interpreter","latex")
title("\boldmath{$T_{xy}^{NN}$}", "Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
colorbar
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
% xlim([-6 * (solution.domain.L) / (solution.domain.domain_coeff), -1 * (solution.domain.L) / (solution.domain.domain_coeff)])
ylim([0 h(251) * 2.5])
% set(gca, "Yscale", "log")
% 
% figure(1);  subplot(3,2,6)
% xlabel("\boldmath{$x$}", "Interpreter","latex")
% ylabel("\boldmath{$y$}", "Interpreter","latex")
% title("\boldmath{$T_{yy}^{NN}$}", "Interpreter","latex")
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'fontsize',20)
% colorbar
% xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
% % xlim([-6 * (solution.domain.L) / (solution.domain.domain_coeff), -1 * (solution.domain.L) / (solution.domain.domain_coeff)])
% ylim([0 h(251) * 2.5])
% % set(gca, "Yscale", "log")

figure(7); 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\boldmath{$y$}", "Interpreter","latex")
title("\boldmath{$dT_{xx}dx^{NN}$}", "Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',20)
colorbar
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
% xlim([-6 * (solution.domain.L) / (solution.domain.domain_coeff), -1 * (solution.domain.L) / (solution.domain.domain_coeff)])
ylim([0 h(251) * 2.5])


figure(8); 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\boldmath{$y$}", "Interpreter","latex")
title("\boldmath{$dT_{xy}dy^{NN}$}", "Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
colorbar
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
ylim([0 h(251) * 2.5])

figure(9); 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\boldmath{$y$}", "Interpreter","latex")
title("\boldmath{$dT_{xx}dx^{NN} + dT_{xy}dy^{NN}$}", "Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
colorbar
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
ylim([0 h(251) * 2.5])

figure(10); 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\boldmath{$y$}", "Interpreter","latex")
title("\boldmath{$pT$}", "Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
colorbar
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
ylim([0 h(251) * 2.5])

Txx_plot = figure(1);
% Txy_plot = figure(2);
% Tyy_plot = figure(3);
% TDxx_plot = figure(4);
% TDxy_plot = figure(5);
% TDyy_plot = figure(6);
DTDxxDx_plot = figure(7);
dTDxydy_plot = figure(8);

% set(Txx_plot,'Position',[160 90 980 720])
% set(Txy_plot,'Position',[160 90 980 720])
% set(Tyy_plot,'Position',[160 90 980 720])
% set(TDxx_plot,'Position',[160 90 980 720])
% set(TDxy_plot,'Position',[160 90 980 720])
% set(TDyy_plot,'Position',[160 90 980 720])
% set(DTDxxDx_plot,'Position',[160 90 980 720])
% set(dTDxydy_plot,'Position',[160 90 980 720])

% saveas(Txx_plot, savepath + "T_xx.fig")
% saveas(Txy_plot, savepath + "T_xy.fig")
% saveas(Tyy_plot, savepath + "T_yy.fig")
% saveas(TDxx_plot, savepath + "TD_xx.fig")
% saveas(TDxy_plot, savepath + "TD_xy.fig")
% saveas(TDyy_plot, savepath + "TD_yy.fig")
% saveas(DTDxxDx_plot, savepath + "DTD_xxDx.fig")
% saveas(dTDxydy_plot, savepath + "dTD_xydy.fig")

close(figure(10));close(figure(9)); close(figure(1)), close(figure(8)), close(figure(7)), close(figure(11))

end