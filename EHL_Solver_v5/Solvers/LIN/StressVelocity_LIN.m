function [velocityfield_dic, stressfield_dic, force] = StressVelocity_LIN(solution, fluid, mesh)

De = fluid.De;
p =  solution.pressure0;
pD = solution.pressure1;
h = solution.h;
dhdx = OneDcentraldiff(h, mesh.dx);
hm = trapz(mesh.x, h.^-2) / trapz(mesh.x, h.^-3);    

dpdx = OneDcentraldiff(p, mesh.dx, "CD6");
d2pdx2 =OneDcentraldiff(dpdx, mesh.dx, "CD6");
% d2pdx2 = OneDSecDcentraldiff(p, mesh.dx, "CD2");
% d2pdx2_an = -6 * dhdx .*( -3 * hm./h.^4 + 2 * 1./h.^3);
% dpdx_an = 6 *   (h - hm) ./ h.^3;   

dpDdx = OneDcentraldiff(pD, mesh.dx);
dhdx = OneDcentraldiff(h, mesh.dx);
Y = mesh.Y; J = mesh.N_y; N = mesh.N_x;  

u = 1/ 2 * dpdx .* (Y.^2 - Y .* h) + (1 - Y ./ h);
dudx = 1/2 * d2pdx2 .* (Y.^2 - Y .* h) - 1/2 * dpdx .* Y .* dhdx + Y .* h.^(-2) .* dhdx;
dydy = reshape(repelem(mesh.dy, J), J, N);
v = dydy .* cumtrapz(dudx);

dudy = dpdx .* (Y - h / 2) - 1 ./ h;                                       
dvdy =  - dudx;
dydx = mesh.Jacobian("dydx");
gradv = CD4(v, mesh.dx, mesh.dy);
% DvDx = (gradv{1} + gradv{2} .* mesh.Jacobian("dydx"));
% v = dhdx .* (2 - 3 * hm./h) .* (Y.^3 ./ h.^3 - Y.^2 ./ h.^2);
% dvdy = dhdx .* (2 - 3 * hm./h) .* (3 * Y.^2 ./ h.^3 - 2 * Y ./ h.^2);  
% V = v ./ h - Y .* h.^(-2) .* dhdx .* u;
% gradV = CD4(V, dx, dy);
% gradv = TwoDcentraldiff(v, mesh.dx, mesh.dy);
DvDx = (gradv{1} + gradv{2} .* dydx);
% DvDx = (gradv{1} + dvdy .* mesh.Jacobian("dydx"));
% gradu =  CD4(u, dx, dy);
gradu = TwoDcentraldiff(u, mesh.dx, mesh.dy);
DuDx = gradu{1} + gradu{2} .* dydx;
% surf(X, Y, dudx - DuDx, "linestyle", "none")  


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

% uD_anal_Step = (Y.^2  - Y .* h) ./ h.^3 .* (1/2 * (1 - fluid.beta) * dhdx - (1 - fluid.beta) / 8 * h.^4 .* dhdx .* (dpdx).^2 + 1/2 * h.^3 .* dpDdx );

% surf(uD - uD_anal_Step, "LineStyle", "none")

MASS = mesh.dy .* trapz(uD);
MASS_anal = mesh.dy .* trapz(uD_anal_Step);
ANAL =  (1 - fluid.beta) / 576 * (dpdx(end-1).^2 * (h(end)^2 - h(1)^2) / 2 + 2 * (1/h(end)^2 - 1/h(1)^2) ) .* ( mesh.dx * trapz( h.^(-3) ) )^(-1);

pD_anal =  (1 - fluid.beta) / 48  * ( dpdx(end-1)^2  * (h.^2 - h(1).^2) / 2  +  2 * (1./ h.^2 - 1./ h(1)^2) ) - 12 * ANAL * cumtrapz(mesh.dx ,1./ h.^3) ;

% figure(1)
% plot(dpDdx)
% hold on
dpDdx_ANAL =  (1 - fluid.beta) / 48 * dhdx .* ( h .* dpdx.^2 - 4 ./ h.^3) - 12 ./ h.^3 .* ANAL;
% plot(dpDdx_ANAL)
% plot(OneDcentraldiff(pD, dx))
% figure(2); hold on
% plot(dx * cumtrapz(dpDdx) , "LineWidth" , 1.8)
% plot(dx * cumtrapz(48 * dpDdx_ANAL), "--", "LineWidth" , 1.8)
% plot(pD, "-", "LineWidth" , 1.8)
% plot(48 * pD_anal, ":", "LineWidth", 1.8)
% figure(3); hold on
% plot( (1 - solution.viscocity_ratio) / 48 * dhdx .* ( h .* dpdx.^2 - 4 ./ h.^3) , "--")
% plot(- 12 ./ h.^3 .* MASS_anal, "--")


graduD = CD4(uD, mesh.dx, mesh.dy);
duDdx = graduD{1};
duDdy = graduD{2};
DuDDx = duDdx + duDdy .* mesh.Jacobian("dydx");
vD = dydy .* cumtrapz(DuDDx);
gradvD = CD4(vD, mesh.dx, mesh.dy);

TD_xx = 2 * T_xy .* dudy;
TD_xy = -(u .* DTxyDx+ v .* dTxydy - T_yy .* dudy) + (1 - fluid.beta) * graduD{2};
TD_yy = -(u .* DTyyDx + v .* dTyydy - 2 * T_xy .* DvDx - 2 * T_yy .* dvdy) + 2 *fluid.beta * gradvD{2};

Wax = trapz(mesh.x, h .* (dpdx + De * dpDdx));
% velocityfield_cell = containers.Map({"u","v","uD","vD","dudx","dvdx","dudy","dvdy"}, ...

velocityfield_dic = containers.Map({char('u'), char('v'), char('uD'), char('vD'), char('DuDx'), char('DvDx'), char('dudy'), char('dvdy')}, ...
                            {u, v, uD, vD, dudx, DvDx, dudy, dvdy});
stressfield_dic = containers.Map({char('tau_xx'), char('tau_xy'), char('tau_yy'), char('tauD_xx'), char('tauD_xy'), char('tauD_yy')}, ... 
                                   {T_xx, T_xy, T_yy, TD_xx, TD_xy, TD_yy});
force = {Wax};

end