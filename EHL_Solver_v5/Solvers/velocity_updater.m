function velocityfield_dic = velocity_updater(h, mesh, fluid, p,DTxxDx, T_xy, dudy)


dpdx = OneDcentraldiff(p, mesh.dx, "CD2");

beta = fluid.beta;
T_xy = T_xy;
[m, n] = size(T_xy);
dy = mesh.dy;
for i = 1:m
    dydy(i,:) = mesh.dy;
end

TD = mesh.Y ./ h; % Transformed Grid 

I11 = TD .* (dy .* trapz(dydy .*cumtrapz(flip(DTxxDx)))) - dydy .* cumtrapz(dydy .* cumtrapz(flip(DTxxDx)));
I12 = TD .* (dy .* trapz(flip(T_xy))) -  dydy .* cumtrapz((flip(T_xy)));


u_k = 1 / 2 * dpdx .* (mesh.Y.^2 - mesh.Y .* h) + (1 - mesh.Y ./ h);
u = 1 / 2 * dpdx .* (mesh.Y.^2  -  mesh.Y .* h) + (1 - mesh.Y ./ h) ...
    + I11 + I12; 

close all
figure(4); hold on
plot(u_k(16, :) , "-.")
plot(u(16, :), "--")
plot(I11(16, :) + I12(16, :), "k:")

gradu = CD4(u, mesh.dx, mesh.dy, "CD4");
gradu = TwoDcentraldiff(u, mesh.dx, mesh.dy);
dudx = gradu{1};
dudy = gradu{2};
DuDx = dudx + dudy .* mesh.Jacobian("dydx");
v = dydy .* cumtrapz(DuDx);
% gradv = CD4(v, mesh.dx, mesh.dy, "CD4");
gradv = TwoDcentraldiff(v, mesh.dx, mesh.dy);
DvDx = gradv{1} + gradv{2} .* mesh.Jacobian("dydx");
dvdy = gradv{2};

step_flag = "off";
if step_flag == "on"
    disp("Step Approximation")
end

uD = 0;
vD = 0;

velocityfield_dic = containers.Map({char('u'), char('v'), char('uD'), char('vD'), char('DuDx'), char('DvDx'), char('dudy'), char('dvdy')}, ...
                            {u, v, uD, vD, DuDx, DvDx, dudy, dvdy});

end