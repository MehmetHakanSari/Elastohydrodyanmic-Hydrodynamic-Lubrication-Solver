function [friction_cell, velocityfield_cell, stressfield_cell, force] = StressVelocity_LIN(solution, N_x, N_y)
% p: zero order pressure, 1D array
% pD: first order presusre: 1D array
% h: height  , 1D array
% viscos: mu_solvent / (mu_solvet + mu_polymer)
% De: deborah number
% Ux: velocity of the channel
% mu: viscosity of the solvent
% h0: reference height to perform thin film approximation

N = N_x;
J = N_y;

De = solution.domain.epsilon * solution.wiessenberg_Number;  

p =  solution.pressure0_old;
pD = solution.pressure1_old;
h = solution.h;

if length(solution.h) ~= N
    h = interp1(linspace(0,1,length(h)), h, linspace(0,1,N));
    p = interp1(linspace(0,1,length(p)), p, linspace(0,1,N));
    pD = interp1(linspace(0,1,length(pD)), pD, linspace(0,1,N));
end

%non-dimensionalize
p = p * solution.domain.href^2  / (solution.domain.mu * solution.velocity * solution.domain.L); 
pD = pD * solution.domain.href^2  / (solution.domain.mu * solution.velocity * solution.domain.L); 
h = h / solution.domain.href;

x_nondimensional = linspace(0,1,N);
dx = 1 / (N-1);

%GRIDING
X = zeros(J, N);
Y = zeros(J, N);

for i = 1 : J
    X(i,:) = x_nondimensional;
end
dy = h / (J - 1);
for i = 1 : J
    Y(i,:) = h - dy .* (i-1);
end


dhdx = OneDcentraldiff(h, dx);
grady = TwoDcentraldiff(Y, dx, dy);
gradx = TwoDcentraldiff(X, dx ,dy);
dydx = -Y ./ h.^2 .* dhdx;          %Coordiante transform y' = y/h where the rest of the code uses y' as y.
dydx = -Y ./ h .* dhdx;             %Modified by multiplying "h". where it come from [d(.)/dy' = h*d(.)/dy]

hm = trapz(x_nondimensional, h.^-2) / trapz(x_nondimensional, h.^-3); 
dpdx = OneDcentraldiff(p, dx);
dpDdx = OneDcentraldiff(pD, dx);
d2pdx2 = OneDcentraldiff(dpdx, dx);

% dpdx = 6 *   (h - hm) ./ h.^3;                        %pressure derivative in l order
% d2pdx2 = -6 * dhdx .*( -3 * hm./h.^4 + 2 * 1./h.^3);
% gradu = TwoDcentraldiff(u, dx, dy);
% DuDx = (gradu{1} + gradu{2} .* dydx);
% v_AN_NUM = dydy .* cumtrapz(dudx);
% gradv_AN_NUM = TwoDcentraldiff(v_AN_NUM, dx, dy);
% v_AN = dhdx .* (2 - 3 * hm./h) .* (Y.^3 ./ h.^3 - Y.^2 ./ h.^2);
% dvdy_AN = dhdx .* (2 - 3 * hm./h) .* (3 * Y.^2 ./ h.^3 - 2 * Y ./ h.^2); 


u = 1/ 2 * dpdx .* (Y.^2 - Y .* h) + (1 - Y ./ h);
dudy = dpdx .* (Y - h / 2) - 1 ./ h;  
dudx = 1/2 * d2pdx2 .* (Y.^2 - Y .* h) - 1/2 * dpdx .* Y .* dhdx + Y .* h.^(-2) .* dhdx;

dvdy =  - dudx;
dydy = reshape(repelem(dy, J), J, N);
v = dydy .* cumtrapz(dudx);
gradv = TwoDcentraldiff(v, dx, dy);
DvDx = (gradv{1} + gradv{2} .* dydx);

% figure(1)
% surf(X, Y, v_AN, "LineStyle", "None")
% colorbar
% figure(2)
% surf(X, Y, v_AN_NUM, "LineStyle", "None")
% colorbar
% figure(3)
% surf(X, Y, dvdy - gradv{2}, "LineStyle", "None")
% colorbar
% figure(4)
% surf(X, Y, dudx + gradv_AN_NUM{2}, "LineStyle", "None")
% colorbar
% figure(5)
% surf(X, Y, dudx + gradv{2}, "LineStyle", "None")
% colorbar

% plot(dpdx); hold on
% plot(dpdx_AN)
% surf(x_MAT, y_MAT, u, "LineStyle", "None")

T_xx = zeros(J,N);
T_xy = (1 - solution.viscocity_ratio) .* dudy;
T_yy = 2 * (1 - solution.viscocity_ratio) .* dvdy;

gradT_xy = TwoDcentraldiff(T_xy, dx, dy);
gradT_yy = TwoDcentraldiff(T_yy, dx, dy);

DTxyDx = gradT_xy{1} + gradT_xy{2} .* dydx;
DTyyDx = gradT_yy{1} + gradT_yy{2} .* dydx;
dTxydy = gradT_xy{2};
dTyydy = gradT_yy{2};

uD = h.^2 / 2 .* dpDdx .* (Y.^2 ./ h.^2 - Y ./ h);
uD = uD + 1 ./ h .* dhdx .* (1 - 3 * hm./h) .* (2 - 3 * hm./h) .* (Y.^2 ./ h.^2 - Y ./ h); 

graduD = TwoDcentraldiff(uD, dx, dy);
duDdx = graduD{1};
duDdy = graduD{2};
DuDDx = duDdx + duDdy .* dydx;

%By continuity duDdx + dvDdy = 0
%Taking derivative then integrating over y would give me vD.
vD = dydy .* cumtrapz(DuDDx);
gradvD = TwoDcentraldiff(vD, dx, dy);

TD_xx = 2 * T_xy .* dudy;
TD_xy = -(u .* DTxyDx+ v .* dTxydy - T_yy .* dudy) + (1 - solution.viscocity_ratio) * graduD{2};
TD_yy = -(u .* DTyyDx + v .* dTyydy - 2 * T_xy .* DvDx - 2 * T_yy .* dvdy) + 2 * solution.viscocity_ratio * gradvD{2};

T_xy_full = dudy;

% friction_top = (trapz(x,T_xy_full(1,:)) + De * trapz(x,T_xy_full(1,:))) * Ux * mu / h0;
% friction_bottom = (trapz(x,T_xy(end,:)) + De * trapz(x,TD_xy(end,:))) * Ux * mu / h0;

x_phys = solution.domain.x;

if length(x_phys) ~= N
    x_phys = interp1(linspace(0,1,length(x_phys)), x_phys, linspace(0,1,N));
end

%Dimensional friction
friction_top = - (trapz(x_phys,T_xy_full(1,:)) + De * trapz(x_phys,TD_xy(1,:))) * solution.velocity * solution.domain.mu / solution.domain.href;
friction_bottom = (trapz(x_phys,T_xy_full(end,:)) + De * trapz(x_phys,TD_xy(end,:))) * solution.velocity * solution.domain.mu / solution.domain.href;
% friction_viscoelastic = (- De * trapz(x_phys,TD_xy(1,:)) + De * trapz(x_phys,TD_xy(end,:))) * solution.velocity * solution.domain.mu / solution.domain.href;
friction_viscoelastic_bottom = De * trapz(x_phys,TD_xy(end,:)) * solution.velocity * solution.domain.mu / solution.domain.href;
friction_viscoelastic_top = - De * trapz(x_phys,TD_xy(1,:)) * solution.velocity * solution.domain.mu / solution.domain.href;


if isempty(solution.applied_load)
    generated_load = (x_phys(2) - x_phys(1)) * trapz(p);
    friction_coeff = - (friction_bottom) / generated_load;
else
    friction_coeff = - (friction_bottom) / solution.applied_load;
end

Wax = trapz(x_phys, h .* (dpdx + De * dpDdx));

% figure(1)
% surf(x_MAT,y_MAT,T_xy,"linestyle","none")
% xlim([0.4 0.6])
% ylim([0 0.1])
% colorbar
% 
% figure(2)
% surf(x_MAT,y_MAT,De * TD_xy,"linestyle","none")
% xlim([0.4 0.6])
% ylim([0 0.1])
% colorbar

% friction_top = (trapz(x_phys_narrowed,T_xy_full(1,226:276)) + De * trapz(x_phys_narrowed,T_xy_full(1,226:276))) * Ux * mu / h0;
% friction_bottom = (trapz(x_phys_narrowed,T_xy(end,226:276)) + De * trapz(x_phys_narrowed,TD_xy(end,226:276))) * Ux * mu / h0;

% disp("Generated Load: " + string(trapz(x*l, p*(U*eta0*l/h0^2) + p0)))
% disp("Friction top: " + string(trapz(x*l, T_xy(1,:)*(U*eta0/h0))))
% disp("Friction bottom: " + string(trapz(x*l, T_xy(end,:)*(U*eta0/h0))))

friction_cell = {friction_top, friction_bottom, friction_viscoelastic_top, friction_coeff, friction_viscoelastic_bottom};
velocityfield_cell = {u, v; uD, vD; dudx, DvDx; dudy, dvdy};
stressfield_cell = {T_xx, T_xy, T_yy; TD_xx, TD_xy, TD_yy};
force = {Wax};

end