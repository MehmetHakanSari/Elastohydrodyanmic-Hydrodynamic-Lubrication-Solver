function stress_firstorder_xy_plotter(u, v, dvdx, dudy, Txx, Tyy, dTxydx, dTxydy, dy, x, viscos, De)
%dy is non dimensional y spacing along x. 
%x is non dimensional "x" w.r.t to Hertiazn semi width

fl_i = round(2 * length(u(1,:)) / 3);
sl_i = round(1.5 * length(u(1,:)) / 3);
tl_i = round(1 * length(u(1,:)) / 3);

u_dTxydx = dy .* trapz(u .* dTxydx);
v_dTxydy = dy .* trapz(v .* dTxydy);

Tyy_dudy = dy .* trapz(Tyy .* dudy);
Txx_dvdx = dy .* trapz(Txx .* dvdx);

eta_dudy = dy .* trapz((1-viscos) * dudy);

colors = ["#0c2c84", "#990000", "#005824", "#6e016b", "#016450"];
lstyle = ":";
figure(2); hold on; box on;
plot(x, -u_dTxydx, lstyle ,"LineWidth", 1.9, "color", colors(1));
plot(x, -v_dTxydy, lstyle ,"LineWidth", 1.9, "color", colors(2));
plot(x, Tyy_dudy, lstyle ,"LineWidth", 1.9, "color", colors(3));
plot(x, Txx_dvdx, lstyle ,"LineWidth", 1.9, "color", colors(4));
plot(x, eta_dudy, lstyle ,"LineWidth", 1.9, "color", colors(5));
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
set(gca, "Fontsize", 20, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex')



legend("\boldmath{$u \cdot \frac{\partial \tau_{xy}}{\partial x}$}",...
        "\boldmath{$v \cdot \frac{\partial \tau_{xy}}{\partial y}$}", ...
        "\boldmath{$\tau_{yy} \cdot \frac{\partial u}{\partial y}$}", ...
        "\boldmath{$\tau_{xx} \cdot \frac{\partial v}{\partial x}$}", ...
        "\boldmath{$\frac{\eta_{p}}{De} \cdot \frac{\partial u}{\partial y}$}", ...
        "box", "off", "Interpreter", "latex")

figure(3); hold on; box on;
subplot(1,3,1)
plot(x, -(u.*dTxydx(fs_i,:)), lstyle ,"LineWidth", 1.9, "color", colors(1));
plot(x, -(v.*dTxydy(fs_i,:)), lstyle ,"LineWidth", 1.9, "color", colors(2));
plot(x, Tyy.*dudy(fs_i,:), lstyle ,"LineWidth", 1.9, "color", colors(3));
% plot(x, Txx.*dvdx, lstyle ,"LineWidth", 1.9, "color", colors(4));
plot(x, eta.*dudy(fs_i,:), lstyle ,"LineWidth", 1.9, "color", colors(5));

end