function stress_zeroorder_xy_plotter(u, v, dudy, duDdy, Tyy, dTxydx, dTxydy, dy, x, viscos, TDxy)
%dy is non dimensional y spacing along x. 
%x is non dimensional "x" w.r.t to Hertiazn semi width

fl_i = round(2 * length(u(:,1)) / 3);
sl_i = round(1.5 * length(u(:,1)) / 3);
tl_i = round(1 * length(u(:,1)) / 3);

u_dTxydx = dy .* trapz(u .* dTxydx);
v_dTxydy = dy .* trapz(v .* dTxydy);

Tyy_dudy = dy .* trapz(Tyy .* dudy);

eta_dudy = dy .* trapz((1-viscos) * duDdy);


colors = ["#0c2c84", "#990000", "#005824", "#6e016b", "#016450"];
lstyle = "-";
figure(2); hold on; box on;
plot(x, -u_dTxydx, lstyle ,"LineWidth", 1.9, "color", colors(1));
plot(x, -v_dTxydy, lstyle ,"LineWidth", 1.9, "color", colors(2));
plot(x, Tyy_dudy, lstyle ,"LineWidth", 1.9, "color", colors(3));
plot(x, eta_dudy, lstyle ,"LineWidth", 1.9, "color", colors(4));
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
set(gca, "Fontsize", 20, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex')


legend("\boldmath{$u \cdot \frac{\partial \tau_{xy}}{\partial x}$}",...
        "\boldmath{$v \cdot \frac{\partial \tau_{xy}}{\partial y}$}", ...
        "\boldmath{$\tau_{yy} \cdot \frac{\partial u}{\partial y}$}", ...
        "\boldmath{$\eta_{p} \cdot \frac{\partial u^D}{\partial y}$}", ...
        "box", "off", "Interpreter", "latex")

SF = max(abs(TDxy), [], "all"); %scaling factor    


u_dTxydx = u .* dTxydx / SF;
v_dTxydy = v .* dTxydy / SF;

Tyy_dudy = Tyy .* dudy / SF;

eta_dudy = (1-viscos) * duDdy / SF;

Txy_sum = Tyy_dudy + eta_dudy - u_dTxydx - v_dTxydy;


% figure(61); hold on;
% plot(x, trapz(eta_dudy), lstyle ,"LineWidth", 1.9, "color", colors(1));
    
figure(31); hold on; box on;
% tl = tiledlayout(3,1,'Tilespacing','tight', 'Padding', 'tight');
lstyle = "--";
% t1 = nexttile(3); hold on;
subplot(4,1,1); hold on; box on; set(gca, "Fontsize", 14, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex') 
plot(x, -(u_dTxydx(fl_i,:)), lstyle ,"LineWidth", 1.9, "color", colors(1));
plot(x, -(v_dTxydy(fl_i,:)), lstyle ,"LineWidth", 1.9, "color", colors(2));
plot(x, Tyy_dudy(fl_i,:), lstyle ,"LineWidth", 1.9, "color", colors(3));
% plot(x, Txx.*dvdx, lstyle ,"LineWidth", 1.9, "color", colors(4));
plot(x, eta_dudy(fl_i,:), lstyle ,"LineWidth", 1.9, "color", colors(4));
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
ylabel("\boldmath{$2h/3$}","Interpreter", "latex")
xlim([-3 3])
% t2 = nexttile(2); hold on;
subplot(4,1,2); hold on; box on; set(gca, "Fontsize", 14, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex') 
plot(x, -(u_dTxydx(sl_i,:)), lstyle ,"LineWidth", 1.9, "color", colors(1));
plot(x, -(v_dTxydy(sl_i,:)), lstyle ,"LineWidth", 1.9, "color", colors(2));
plot(x, Tyy_dudy(sl_i,:), lstyle ,"LineWidth", 1.9, "color", colors(3));
% plot(x, Txx.*dvdx, lstyle ,"LineWidth", 1.9, "color", colors(4));
plot(x, eta_dudy(sl_i,:), lstyle ,"LineWidth", 1.9, "color", colors(4));
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
ylabel("\boldmath{$h/2$}","Interpreter", "latex")
xlim([-3 3])
% t3 = nexttile(1); hold on;
subplot(4,1,3); hold on; box on; set(gca, "Fontsize", 14, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex') 
plot(x, -(u_dTxydx(tl_i,:)), lstyle ,"LineWidth", 1.9, "color", colors(1));
plot(x, -(v_dTxydy(tl_i,:)), lstyle ,"LineWidth", 1.9, "color", colors(2));
plot(x, Tyy_dudy(tl_i,:), lstyle ,"LineWidth", 1.9, "color", colors(3));
% plot(x, Txx.*dvdx, lstyle ,"LineWidth", 1.9, "color", colors(4));
plot(x, eta_dudy(tl_i,:), lstyle ,"LineWidth", 1.9, "color", colors(4));
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
ylabel("\boldmath{$h/3$}","Interpreter", "latex")
xlim([-3 3])


set(gca, "Fontsize", 14, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex')
legend("\boldmath{$u \cdot \frac{\partial \tau_{xy}}{\partial x}$}",...
        "\boldmath{$v \cdot \frac{\partial \tau_{xy}}{\partial y}$}", ...
        "\boldmath{$\tau_{yy} \cdot \frac{\partial u}{\partial y}$}", ...
        "\boldmath{$(1-\beta) \cdot \frac{\partial u^D}{\partial y}$}", ...
        "box", "off", "Interpreter", "latex")
    
    

figure(32); hold on; box on;
% tl = tiledlayout(3,1,'Tilespacing','tight', 'Padding', 'tight');
lstyle = "--";
% t1 = nexttile(3); hold on;
subplot(4,1,1); hold on; box on; set(gca, "Fontsize", 14, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex') 
plot(x, -(u_dTxydx(1,:)), lstyle ,"LineWidth", 1.9, "color", colors(1));
plot(x, -(v_dTxydy(1,:)), lstyle ,"LineWidth", 1.9, "color", colors(2));
plot(x, Tyy_dudy(1,:), lstyle ,"LineWidth", 1.9, "color", colors(3));
% plot(x, Txx.*dvdx, lstyle ,"LineWidth", 1.9, "color", colors(4));
plot(x, eta_dudy(1,:), lstyle ,"LineWidth", 1.9, "color", colors(4));
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
ylabel("\boldmath{$y = h$}","Interpreter", "latex")
xlim([-3 3])
subplot(4,1,2); hold on; box on; set(gca, "Fontsize", 14, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex') 
plot(x, -(u_dTxydx(sl_i,:)), lstyle ,"LineWidth", 1.9, "color", colors(1));
plot(x, -(v_dTxydy(sl_i,:)), lstyle ,"LineWidth", 1.9, "color", colors(2));
plot(x, Tyy_dudy(sl_i,:), lstyle ,"LineWidth", 1.9, "color", colors(3));
% plot(x, Txx.*dvdx, lstyle ,"LineWidth", 1.9, "color", colors(4));
plot(x, eta_dudy(sl_i,:), lstyle ,"LineWidth", 1.9, "color", colors(4));
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
ylabel("\boldmath{$h/2$}","Interpreter", "latex")
xlim([-3 3])
% t3 = nexttile(1); hold on;
subplot(4,1,3); hold on; box on; set(gca, "Fontsize", 14, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex') 
plot(x, -(u_dTxydx(end,:)), lstyle ,"LineWidth", 1.9, "color", colors(1));
plot(x, -(v_dTxydy(end,:)), lstyle ,"LineWidth", 1.9, "color", colors(2));
plot(x, Tyy_dudy(end,:), lstyle ,"LineWidth", 1.9, "color", colors(3));
% plot(x, Txx.*dvdx, lstyle ,"LineWidth", 1.9, "color", colors(4));
plot(x, eta_dudy(end,:), lstyle ,"LineWidth", 1.9, "color", colors(4));
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
ylabel("\boldmath{$y = 0$}","Interpreter", "latex")
xlim([-3 3])


set(gca, "Fontsize", 14, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex')
legend("\boldmath{$u \cdot \frac{\partial \tau_{xy}}{\partial x}$}",...
        "\boldmath{$v \cdot \frac{\partial \tau_{xy}}{\partial y}$}", ...
        "\boldmath{$\tau_{yy} \cdot \frac{\partial u}{\partial y}$}", ...
        "\boldmath{$(1-\beta) \cdot \frac{\partial u^D}{\partial y}$}", ...
        "box", "off", "Interpreter", "latex")
    
    
    


eta_du0dy = (1-viscos) * dudy / SF;


% figure(41); hold on; box on;
% tl = tiledlayout(3,1,'Tilespacing','tight', 'Padding', 'tight');
% lstyle = "--";
% t1 = nexttile(3); hold on;
% plot(x, eta_du0dy(fl_i,:), lstyle ,"LineWidth", 1.9, "color", colors(4));
% xlim([-4 1.1])
% t2 = nexttile(2); hold on;
% plot(x, eta_du0dy(sl_i,:), lstyle ,"LineWidth", 1.9, "color", colors(4));
% xlim([-4 1.1])
% t3 = nexttile(1); hold on;
% plot(x, eta_du0dy(tl_i,:), lstyle ,"LineWidth", 1.9, "color", colors(4));
% xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
% xlim([-4 1.1])
% 
% 
% set(gca, "Fontsize", 14, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex')
% legend("\boldmath{$(1-\beta) \cdot \frac{\partial u}{\partial y}$}", ...
%         "box", "off", "Interpreter", "latex")


% lstyle = "--";
figure(13); 
subplot(4,1,1); hold on; box on; title("\boldmath{$\tau_{xy}^{NN}$}", "Interpreter","latex"); set(gca, "Fontsize", 14, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex') 
plot(x, Txy_sum(1,:), lstyle ,"LineWidth", 1.9)
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
ylabel("\boldmath{$y = h$}","Interpreter", "latex")
xlim([-3 3])
set(gca, "Fontsize", 14)
subplot(4,1,2); hold on; box on;set(gca, "Fontsize", 14, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex') 
plot(x, Txy_sum(sl_i,:), lstyle ,"LineWidth", 1.9)
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
ylabel("\boldmath{$y = h/2$}","Interpreter", "latex")
xlim([-3 3])
set(gca, "Fontsize", 14)
subplot(4,1,3); hold on; box on;set(gca, "Fontsize", 14, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex') 
plot(x, Txy_sum(end,:), lstyle ,"LineWidth", 1.9)
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
ylabel("\boldmath{$y = 0$}","Interpreter", "latex")
xlim([-3 3])
set(gca, "Fontsize", 14)


% figure(4)
% surf(dTxydy, "linestyle", "none")
end