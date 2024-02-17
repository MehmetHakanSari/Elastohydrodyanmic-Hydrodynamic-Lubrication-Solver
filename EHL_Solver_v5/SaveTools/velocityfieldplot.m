function velocityfieldplot(solution, savepath)

file_name = "L";

file_name = file_name + string(solution.applied_load);
file_name = file_name + "_Wi" + string(solution.wiessenberg_Number);
file_name = file_name + "_beta" + string(solution.viscocity_ratio);
file_name = file_name + "_U" + string(solution.velocity);
file_name = file_name + "_bH" + string(solution.domain.domain_coeff);


if ~isfolder("Results/Figures/" + file_name)
       mkdir("Results/Figures/" + file_name);
end

if ~isfolder("Results/Figures/" + file_name + "/velocityfield")
       mkdir("Results/Figures/" + file_name + "/velocityfield");
end

savepath = savepath + file_name + "\velocityfield\";

[N, J] = size(solution.velocity_field{1,1});
h = solution.h;

%GRIDING
X = zeros(J, N);
Y = zeros(J, N);

for i = 1:J
    X(i,:) = solution.domain.x;
end
dy = solution.h / (J - 1);
for i = 1:J
    Y(i,:) = solution.h - dy .* (i-1);
end



% solution.domain.epsilon * solution.wiessenberg_Number * (1 - solution.viscosity_ratio)

figure(1);
surf(X, Y, solution.velocity_field{1,1}, "LineStyle", "None");
figure(2); hold on
surf(X, Y, solution.velocity_field{1,2}, "LineStyle", "None");

figure(3);
surf(X, Y, solution.velocity_field{2,1} * solution.domain.epsilon * solution.wiessenberg_Number, "LineStyle", "None");
figure(4); hold on
surf(X, Y, solution.velocity_field{2,2} * solution.domain.epsilon * solution.wiessenberg_Number, "LineStyle", "None");

figure(5);
surf(X, Y, solution.velocity_field{1,1} + solution.velocity_field{2,1} * solution.domain.epsilon * solution.wiessenberg_Number, "LineStyle", "None");
figure(6); hold on
surf(X, Y, solution.velocity_field{1, 2} + solution.velocity_field{2,2} * solution.domain.epsilon * solution.wiessenberg_Number, "LineStyle", "None");

gradu = TwoDcentraldiff(solution.velocity_field{1,1}, solution.domain.dx / solution.domain.L, solution.h / (J - 1) / solution.domain.href);
graduD = TwoDcentraldiff(solution.velocity_field{2,1}, solution.domain.dx / solution.domain.L, solution.h / (J - 1) / solution.domain.href);

figure(7);
surf(X, Y, gradu{2}, "LineStyle", "None");
figure(8);
surf(X, Y, graduD{2} * solution.domain.epsilon * solution.wiessenberg_Number, "LineStyle", "None");


figure(1); 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\boldmath{$y$}", "Interpreter","latex")
title("\boldmath{$u^{N}$}", "Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
colorbar
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
ylim([0 h(251) * 2.5])


figure(2); 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\boldmath{$y$}", "Interpreter","latex")
title("\boldmath{$v^{N}$}", "Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
colorbar
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
ylim([0 h(251) * 2.5])


figure(3); 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\boldmath{$y$}", "Interpreter","latex")
title("\boldmath{$u^{NN}$}", "Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
colorbar
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
ylim([0 h(251) * 2.5])

figure(4); 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\boldmath{$y$}", "Interpreter","latex")
title("\boldmath{$v^{NN}$}", "Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
colorbar
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
ylim([0 h(251) * 2.5])


figure(5); 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\boldmath{$y$}", "Interpreter","latex")
title("\boldmath{$u^{T}$}", "Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
colorbar
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
ylim([0 h(251) * 2.5])

figure(6); 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\boldmath{$y$}", "Interpreter","latex")
title("\boldmath{$v^{T}$}", "Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
colorbar
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
ylim([0 h(251) * 2.5])

figure(7); 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\boldmath{$y$}", "Interpreter","latex")
title("\boldmath{$dudy^{N}$}", "Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
colorbar
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
ylim([0 h(251) * 2.5])

figure(8); 
xlabel("\boldmath{$x$}", "Interpreter","latex")
ylabel("\boldmath{$y$}", "Interpreter","latex")
title("\boldmath{$dudy^{NN}$}", "Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)
colorbar
xlim([-3 * (solution.domain.L) / (solution.domain.domain_coeff), 3 * (solution.domain.L) / (solution.domain.domain_coeff)])
ylim([0 h(251) * 2.5])



u_plot = figure(1);
v_plot = figure(2);
uD_plot = figure(3);
vD_plot = figure(4);
uT_plot = figure(5);
vT_plot = figure(6);
dudy_plot = figure(7);
duDdy_plot = figure(8);


set(u_plot,'Position',[160 90 980 720])
set(v_plot,'Position',[160 90 980 720])
set(uD_plot,'Position',[160 90 980 720])
set(vD_plot,'Position',[160 90 980 720])
set(uT_plot,'Position',[160 90 980 720])
set(vT_plot,'Position',[160 90 980 720])
set(dudy_plot,'Position',[160 90 980 720])
set(duDdy_plot,'Position',[160 90 980 720])

saveas(u_plot, savepath + "u.fig")
saveas(v_plot, savepath + "v.fig")
saveas(uD_plot, savepath + "uD.fig")
saveas(vD_plot, savepath + "vD.fig")
saveas(uT_plot, savepath + "uT.fig")
saveas(vT_plot, savepath + "vT.fig")
saveas(dudy_plot, savepath + "dudy.fig")
saveas(duDdy_plot, savepath + "duDdy.fig")

% gradu = TwoDcentraldiff(solution.velocity_field{1,1}, solution.domain.dx, dy);
% dudy = gradu{2};
% u  = solution.velocity_field{1,1};

% figure(5);
% plot(dudy(:,225), linspace(0,1,N), "-", "LineWidth", 2)
% figure(7);
% plot(u(:,220), linspace(0,1,N), "-", "LineWidth", 2)
% figure(8);
% plot(u(:,230), linspace(0,1,N), "-", "LineWidth", 2)
% figure(9);
% plot(u(:,240), linspace(0,1,N), "-", "LineWidth", 2)
% figure(10);
% plot(u(:,250), linspace(0,1,N), "-", "LineWidth", 2)
% figure(11);
% plot(u(:,260), linspace(0,1,N), "-", "LineWidth", 2)



end