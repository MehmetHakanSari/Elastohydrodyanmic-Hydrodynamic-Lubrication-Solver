function [velocityfield_dic, stressfield_dic, force] = StressVelocity_LIN(solution, fluid, mesh)

De = fluid.De;
p =  solution.pressure_LIN("p0");
pD = solution.pressure_LIN("p1");
h = solution.h;

hm = trapz(mesh.x, h.^-2) / trapz(mesh.x, h.^-3);    

dpdx = OneDcentraldiff(p, mesh.dx, "CD6");
d2pdx2 =OneDcentraldiff(dpdx, mesh.dx, "CD6");

dpDdx = OneDcentraldiff(pD, mesh.dx, "CD6");
d2pDdx2 = OneDcentraldiff(dpDdx, mesh.dx, "CD6");
dhdx = OneDcentraldiff(h, mesh.dx);
Y = mesh.Y; J = mesh.N_y; N = mesh.N_x;  

u = 1/ 2 * dpdx .* (Y.^2 - Y .* h) + (1 - Y ./ h);
dudx = 1/2 * d2pdx2 .* (Y.^2 - Y .* h) - 1/2 * dpdx .* Y .* dhdx + Y .* h.^(-2) .* dhdx;
dydy = reshape(repelem(mesh.dy, J), J, N);
v = dydy .* cumtrapz(dudx);

dudy = dpdx .* (Y - h / 2) - 1 ./ h;                                       
dvdy =  - dudx;
dydx = mesh.Jacobian("dydx");
gradv = CD4(v, mesh.dx, mesh.dy, "CD6");
DvDx = (gradv{1} + gradv{2} .* dydx);

gradu = CD4(u, mesh.dx, mesh.dy, "CD6");
DuDx = (gradu{1} + gradu{2} .* dydx);
% 
% figure(1);
% surf(DuDx + gradv{2}, "LineStyle", "none")
% figure(2);
% surf(gradu{2} - dudy, "LineStyle", "none")

T_xx = zeros(J, N);
T_xy = (1 - fluid.beta) .* dudy;
T_yy = 2 * (1 - fluid.beta) .* dvdy;

gradT_xy = CD4(T_xy, mesh.dx, mesh.dy);
gradT_yy = CD4(T_yy, mesh.dx, mesh.dy);

DTxyDx = gradT_xy{1} + gradT_xy{2}.* mesh.Jacobian("dydx");
DTyyDx = gradT_yy{1} + gradT_yy{2}.* mesh.Jacobian("dydx");
dTxydy = gradT_xy{2};
dTyydy = gradT_yy{2};
hm = trapz(h.^-2) / trapz(h.^-3);  
uD = h.^2 / 2 .* dpDdx .* (Y.^2 ./ h.^2 - Y ./ h);
uD = uD +  (1 - fluid.beta) ./ h .* dhdx .* ( 1 - 3 * hm./h ) .* ( 2 - 3 * hm./h ) .* (Y.^2 ./ h.^2 - Y ./ h); 

uD_anal_Step = (1/2 * (1 - fluid.beta) * dhdx ./h - (1 - fluid.beta) / 8 * h.^3 .* dhdx .* (dpdx).^2 + 1/2 * h.^2 .* dpDdx ) .* (Y.^2 ./ h.^2 - Y ./ h);
vD_anal_Step =  (2 * Y.^3 - 3 * Y.^2 .* h) / 12 .* (d2pDdx2) - dpDdx .* dhdx .* Y.^2 / 4;
% uD = (1/2 * (1 - fluid.beta) * dhdx ./h - (1 - fluid.beta) / 8 * h.^3 .* dhdx .* (dpdx).^2 + 1/2 * h.^2 .* dpDdx ) .* (Y.^2 ./ h.^2 - Y ./ h);
% surf(uD - uD_anal_Step, "LineStyle", "none")

if solution.step_flag == "on"
    %Constant U approximation until step
    [M, step_locataion_index] = max(abs(dhdx));
    u(:, 1:step_locataion_index) = u(:, 1:step_locataion_index) + De * uD_anal_Step(:, 10);
    u(:, step_locataion_index+1:end) = u(:, step_locataion_index+1:end) + De * uD_anal_Step(:, end-10);
    v(:, 1:step_locataion_index) = v(:, 1:step_locataion_index) + De * vD_anal_Step(:, 10);
    v(:, step_locataion_index+10:end) = v(:, step_locataion_index+10:end) + De * vD_anal_Step(:, step_locataion_index+10:end);
%     u = u + De * uD_anal_Step;
%     v = v + De * vD_anal_Step;

    disp("First order velocity correction")
    uD = uD_anal_Step;
    vD = vD_anal_Step;
end

graduD = CD4(uD, mesh.dx, mesh.dy);
duDdx = graduD{1};
duDdy = graduD{2};
DuDDx = duDdx + duDdy .* mesh.Jacobian("dydx");
vD = dydy .* cumtrapz(DuDDx);
gradvD = CD4(vD, mesh.dx, mesh.dy);

TD_xx = 2 * T_xy .* dudy;
TD_xy = -(u .* DTxyDx+ v .* dTxydy - T_yy .* dudy) + (1 - fluid.beta) * graduD{2};
TD_yy = -(u .* DTyyDx + v .* dTyydy - 2 * T_xy .* DvDx - 2 * T_yy .* dvdy) + 2 *fluid.beta * gradvD{2};

q = mesh.dy .* trapz(u);
qD = mesh.dy .* trapz(uD);

Wax = trapz(mesh.x, h .* (dpdx + De * dpDdx));

velocityfield_dic = containers.Map({char('u'), char('v'), char('uD'), char('vD'), char('DuDx'), char('DvDx'), char('dudy'), char('dvdy'), char("q"), char("qD")}, ...
                            {u, v, uD, vD, dudx, DvDx, dudy, dvdy, q, qD});
stressfield_dic = containers.Map({char('tau_xx'), char('tau_xy'), char('tau_yy'), char('tauD_xx'), char('tauD_xy'), char('tauD_yy')}, ... 
                                   {T_xx, T_xy, T_yy, TD_xx, TD_xy, TD_yy});
force = containers.Map({char('Wax')},{Wax});

% figure(1);
% surf(De * uD, "linestyle", "none")
% figure(2);
% surf(De * uD_anal_Step, "linestyle", "none")
% figure(3);
% surf(v, "linestyle", "none")
% figure(4);
% surf(De * vD, "linestyle", "none")

end