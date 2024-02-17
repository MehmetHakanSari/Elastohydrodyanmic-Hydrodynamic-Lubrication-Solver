function stress_xx_plotter(u, v, dudx, dudy, Txx, Txy, dTxxdx, dTxxdy, dy, x)
%dy is non dimensional y spacing along x. 
%x is non dimensional "x" w.r.t to Hertiazn semi width


u_dTxxdx = dy .* trapz(u .* dTxxdx);
v_dTxxdy = dy .* trapz(v .* dTxxdy);

Txy_dudy = dy .* trapz(2 * Txy .* dudy);
Txx_dudx = dy .* trapz(2 * Txx .* dudx);

colors = ["#0c2c84", "#990000", "#005824", "#6e016b"];
lstyle = ":";
figure(1); hold on; box on;
plot(x, -u_dTxxdx, lstyle ,"LineWidth", 1.9, "color", colors(1));
plot(x, -v_dTxxdy, lstyle ,"LineWidth", 1.9, "color", colors(2));
plot(x, Txy_dudy, lstyle ,"LineWidth", 1.9, "color", colors(3));
plot(x, Txx_dudx, lstyle ,"LineWidth", 1.9, "color", colors(4));
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
set(gca, "Fontsize", 20, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex')

legend("\boldmath{$u \cdot \frac{\partial \tau_{xx}}{\partial x}$}",...
        "\boldmath{$v \cdot \frac{\partial \tau_{xx}}{\partial y}$}", ...
        "\boldmath{$2 \tau_{xy} \cdot \frac{\partial u}{\partial y}$}", ...
        "\boldmath{$2 \tau_{xx} \cdot \frac{\partial u}{\partial x}$}", ...
        "box", "off", "Interpreter", "latex")

    
figure(3); hold on, box on;
plot(x, dy .* trapz(Txx),":", "color", colors(1), "LineWidth", 1.8)
plot(x, dy .* trapz(Txy),":","color", colors(2), "LineWidth", 1.8)
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
set(gca, "Fontsize", 20, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex')
legend("\boldmath{$\tau_{xx}$}","\boldmath{$\tau_{xy}$}", "box", "off", "Interpreter", "latex")

figure(4); hold on; box on;
plot(x, (-Txy_dudy - Txx_dudx + u_dTxxdx + v_dTxxdy) * 0.005, "-", x, dy .* trapz(Txx), ":", "LineWidth", 1.8)
xlabel("\boldmath{$x/b_H$}","Interpreter", "latex")
set(gca, "Fontsize", 20, "Linewidth", 1.2, 'TickLabelInterpreter', 'latex')
legend("\boldmath{$Sum$}","\boldmath{$\tau_{xx}$}", "box", "off", "Interpreter", "latex")


end