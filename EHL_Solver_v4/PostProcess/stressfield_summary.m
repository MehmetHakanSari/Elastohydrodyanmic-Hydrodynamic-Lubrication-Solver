function stressfield_summary(solution)

% file_name = "L";
% 
% file_name = file_name + string(solution.applied_load);
% file_name = file_name + "_Wi" + string(solution.wiessenberg_Number);
% file_name = file_name + "_beta" + string(solution.viscocity_ratio);
% file_name = file_name + "_U" + string(solution.velocity);
% file_name = file_name + "_bH" + string(solution.domain.domain_coeff);
% 
% 
% if ~isfolder("Results/Figures/" + file_name)
%        mkdir("Results/Figures/" + file_name);
% end
% 
% if ~isfolder("Results/Figures/" + file_name + "/stressfield")
%        mkdir("Results/Figures/" + file_name + "/stressfield");
% end
% 
% savepath = savepath + file_name + "\stressfield\";



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

dydx = -Y ./ h .* OneDcentraldiff(h / solution.domain.href,  solution.domain.dx / solution.domain.L);

gradTDxx = TwoDcentraldiff(solution.stress_field{2,1} , solution.domain.dx / solution.domain.L, solution.h / (J - 1) / solution.domain.href);
gradTDxy = TwoDcentraldiff(solution.stress_field{2,2} , solution.domain.dx / solution.domain.L, solution.h / (J - 1) / solution.domain.href);
gradTxx = TwoDcentraldiff(solution.stress_field{1,1} , solution.domain.dx / solution.domain.L, solution.h / (J - 1) / solution.domain.href);
gradTxy = TwoDcentraldiff(solution.stress_field{1,2} , solution.domain.dx / solution.domain.L, solution.h / (J - 1) / solution.domain.href);
graduD = TwoDcentraldiff(solution.velocity_field{2,1} , solution.domain.dx / solution.domain.L, solution.h / (J - 1) / solution.domain.href);
gradvD = TwoDcentraldiff(solution.velocity_field{2,2} , solution.domain.dx / solution.domain.L, solution.h / (J - 1) / solution.domain.href);

DuDDx = graduD{1} + graduD{2} .* dydx;
DvDDx = gradvD{1} + gradvD{2} .* dydx;

DTDxxDx = gradTDxx{1} + gradTDxx{2} .* dydx;
DTDxyDx = gradTDxy{1} + gradTDxy{2} .* dydx;
DTxxDx = gradTxx{1} + gradTxx{2} .* dydx;
DTxyDx = gradTxy{1} + gradTxy{2} .* dydx;

% First order xx direction element plotter
% stress_xx_plotter(solution.velocity_field{2,1}, solution.velocity_field{2,2}, ...
%                   DuDDx, graduD{2}, ...
%                   solution.stress_field{2,1}, solution.stress_field{2,2}, ...
%                   DTDxxDx, gradTDxx{2}, dy / solution.domain.href, solution.domain.x / solution.domain.bH)
              
%First order xy direction element plotter
% stress_firstorder_xy_plotter(solution.velocity_field{2,1}, solution.velocity_field{2,2}, ...
%                   DvDDx, graduD{2}, ...
%                   solution.stress_field{2,1}, solution.stress_field{2,3}, ...
%                   DDTxyDx, gradTDxy{2}, dy / solution.domain.href, solution.domain.x / solution.domain.bH, ...
%                   solution.viscocity_ratio, solution.deborah_Number)
              
stress_zeroorder_xy_plotter(solution.velocity_field{1,1}, solution.velocity_field{1,2}, ...
                  solution.velocity_field{4,1}, graduD{2}, ...
                  solution.stress_field{1,3}, ...
                  DTxyDx, gradTxy{2}, dy / solution.domain.href, solution.domain.x / solution.domain.bH, ...
                  solution.viscocity_ratio, solution.stress_field{2,2})

              

% solution.stress_field{1,2} + solution.deborah_Number * solution.stress_field{2,2}         
              
              
figure(31); hold on; box on;
subplot(4,1,4); hold on; box on; set(gca, "Fontsize", 14, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex') 
plot(solution.domain.x / solution.domain.bH, solution.h / solution.domain.href,"--" ,"LineWidth", 1.9)
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
ylabel("\bf{h/href}","Interpreter", "latex")
xlim([-3 3])
legend("\bf{De=0.01}", "\bf{De=0.05}", "box", "off", "Interpreter","latex")

figure(32); hold on; box on;
subplot(4,1,4); hold on; box on; set(gca, "Fontsize", 14, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex') 
plot(solution.domain.x / solution.domain.bH, solution.h / solution.domain.href,"--" ,"LineWidth", 1.9)
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
ylabel("\bf{h/href}","Interpreter", "latex")
xlim([-3 3])
legend("\bf{De=0.01}", "\bf{De=0.05}", "box", "off", "Interpreter","latex")
figure(13); hold on; box on;
subplot(4,1,4); hold on; box on; set(gca, "Fontsize", 14, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex') 
plot(solution.domain.x / solution.domain.bH, solution.h / solution.domain.href,"--" ,"LineWidth", 1.9)
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
ylabel("\bf{h/href}","Interpreter", "latex")
xlim([-3 3])
legend("\bf{De=0.01}", "\bf{De=0.05}", "box", "off", "Interpreter","latex")

% subplot(1,2,2); hold on;
% plot(solution.domain.x / solution.domain.bH, solution.pressure, "LineWidth", 1.9)
% xlabel("\bf{x/bH}","Interpreter", "latex")
% ylabel("\bf{p/pH}","Interpreter", "latex")
% set(gca, "Linewidth", 1.1)
% set(gca, "Fontsize", 20)
% set(gca, 'TickLabelInterpreter', 'latex')
% xlim([-3 3])

% dy .* trapz(solution.stress_field{2,1})  
% 
% hm = trapz(solution.h.^-2) / trapz(solution.h.^-3); 
% 
% uT = solution.velocity_field{1,1} + solution.velocity_field{2,1} * solution.deborah_Number;
% vT = solution.velocity_field{1,2} + solution.velocity_field{2,2} * solution.deborah_Number;
% DTTxxDx = DTxxDx + DTDxxDx * solution.deborah_Number;
% dTTxxdy = gradTxx{2} + gradTDxx{2} * solution.deborah_Number;
% convective_term_xx_D = solution.velocity_field{2,1} .* DTDxxDx + solution.velocity_field{2,2} .* gradTDxx{2}; 
% convective_term_xy_D = solution.velocity_field{2,1} .* DTDxyDx + solution.velocity_field{2,2} .* gradTDxy{2}; 
% convective_term_T = uT .* DTTxxDx + vT .* dTTxxdy; 
% row_len = length(convective_term_T(:,1));
% colors = ["#0c2c84", "#990000"];
% % colors = ["#990000", "#0c2c84"];
% ls = ":";
% fs = 24;
% figure(12); 
% subplot(3,1,1); hold on; box on;
% txxline1 = plot(solution.domain.x / solution.domain.bH, convective_term_xx_D(round(row_len/3),:)/ max(abs(convective_term_xx_D),[],"all"),ls,"LineWidth", 1.8, "Color",colors(1));
% xlim([-3 3])
% set(gca,'linewidth',1.2)
% ylabel("\boldmath{$y=2h/3$}","Interpreter", "latex")
% subplot(3,1,3); hold on; box on;
% txxline2 = plot(solution.domain.x / solution.domain.bH, convective_term_xx_D(end,:)/max(abs(convective_term_xx_D),[],"all"),ls,"LineWidth", 1.8,"Color",colors(1));
% xlim([-3 3])
% set(gca,'linewidth',1.2)
% ylabel("\boldmath{$y=0$}","Interpreter", "latex")
% xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
% subplot(3,1,2); hold on; box on;
% txxline3 = plot(solution.domain.x / solution.domain.bH, convective_term_xx_D(round(2*row_len/3),:) / max(abs(convective_term_xx_D),[],"all"),ls,"LineWidth", 1.8,"Color",colors(1));  
% xlim([-3 3])
% set(gca,'linewidth',1.2)
% ylabel("\boldmath{$y=h/3$}","Interpreter", "latex")
% % legend("W = 800N, De = 0.05", "W = 800N, De = 0", "W = 100N, De = 0", "W = 100N, De = 0.05", "box", "off", "Interpreter", "latex")
% % legend("W = 800N, De = 0.05", "W = 800N, De = 0", "W = 100N, De = 0", "W = 100N, De = 0.05", "box", "off", "Interpreter", "latex")
% hh = figure(12);
% set(findall([hh], '-property', 'Interpreter'), 'Interpreter', 'latex')
% set(findall([hh], '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex', 'FontSize', fs-6)
% set(findobj([hh], 'Type', 'legend'), 'FontSize', fs-4)
% % set(findall([txxline1 txxline2 txxline3], '-property', 'Interpreter'), 'Interpreter', 'latex')
% % set(findall([txxline1 txxline2 txxline3], '-property', 'TickLabelInterpreter'), 'TickLabelInterpreter', 'latex', 'FontSize', fs-6)
% % set(findobj([txxline1 txxline2 txxline3], 'Type', 'legend'), 'FontSize', fs-4)
% 
% 
% 
% Txx_plot = figure(1);
% % Txy_plot = figure(2);
% % Tyy_plot = figure(3);
% % TDxx_plot = figure(4);
% % TDxy_plot = figure(5);
% % TDyy_plot = figure(6);
% DTDxxDx_plot = figure(7);
% dTDxydy_plot = figure(8);
% 
% % set(Txx_plot,'Position',[160 90 980 720])
% % set(Txy_plot,'Position',[160 90 980 720])
% % set(Tyy_plot,'Position',[160 90 980 720])
% % set(TDxx_plot,'Position',[160 90 980 720])
% % set(TDxy_plot,'Position',[160 90 980 720])
% % set(TDyy_plot,'Position',[160 90 980 720])
% % set(DTDxxDx_plot,'Position',[160 90 980 720])
% % set(dTDxydy_plot,'Position',[160 90 980 720])
% 
% % saveas(Txx_plot, savepath + "T_xx.fig")
% % saveas(Txy_plot, savepath + "T_xy.fig")
% % saveas(Tyy_plot, savepath + "T_yy.fig")
% % saveas(TDxx_plot, savepath + "TD_xx.fig")
% % saveas(TDxy_plot, savepath + "TD_xy.fig")
% % saveas(TDyy_plot, savepath + "TD_yy.fig")
% % saveas(DTDxxDx_plot, savepath + "DTD_xxDx.fig")
% % saveas(dTDxydy_plot, savepath + "dTD_xydy.fig")
% 
% close(figure(10));close(figure(9)); close(figure(1)), close(figure(8)), close(figure(7)), close(figure(11))

end