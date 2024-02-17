function [pressure_dic] = PressureSolver_v1(h, mesh, fluid, DTxxDx, T_xy)


%-----------------------------------

[m, n] = size(DTxxDx);
dydy = zeros(m, n);
N = n;
dhdx = OneDcentraldiff(h, mesh.dx, "CD2");

for i = 1:m
    dydy(i,:) = mesh.dy;
end

TD = mesh.Y ./ h; % Transformed Grid 

zeta = h.^3 / 12;
W = zeros(1, N);
E = zeros(1, N);
W(2:end) = (zeta(1:end-1) + zeta(2:end)) / (2*mesh.dx^2);
E(1:end-1) = (zeta(1:end-1) + zeta(2:end)) / (2*mesh.dx^2);
Cent = -(E + W);

Q1 = 0.5 *   fluid.beta * dhdx;   %Newtonian stress solvent contrubution
Q4 = 0.5 *  (1- fluid.beta) * dhdx;      %Newtonian stress polymer contrubution
dy = mesh.dy;
I11 = dy .* trapz(TD(:,1) * (dy .* trapz(dydy .*cumtrapz(flip(DTxxDx))))) - dy .* trapz(dydy .* cumtrapz(dydy .* cumtrapz(flip(DTxxDx))));
I12 = dy .* (trapz(TD(:,1) * (dy .* trapz(flip(T_xy)))) -  trapz(dydy .* cumtrapz((flip(T_xy)))));

% plot(I11); hold on
% plot(I12)

Q2 = OneDcentraldiff(I11, mesh.dx);             % xx component of the pressure
Q3 = OneDcentraldiff(I12, mesh.dx);             % xy component of the pressure

p_s = zeros(1, N);
p_xx = zeros(1, N);
p_xy = zeros(1, N);
p_xy_N = zeros(1, N);

p_s(2:end-1) = TDMA(W(2:end-1),Cent(2:end-1),E(2:end-1), Q1(2:end-1));
p_xx(2:end-1) = TDMA(W(2:end-1),Cent(2:end-1),E(2:end-1), Q2(2:end-1));
p_xy(2:end-1) = TDMA(W(2:end-1),Cent(2:end-1),E(2:end-1), Q3(2:end-1));
p_xy_N(2:end-1) = TDMA(W(2:end-1),Cent(2:end-1),E(2:end-1), Q4(2:end-1));

p_xy_NN = p_xy - p_xy_N;         %Non-Newtonian pressure in xy component
p_ve = p_xx + p_xy_NN;           %Total Non-newtonian pressure

%Load
p = p_s + p_xx + p_xy;
pressure_dic = containers.Map({char('p'), char('p_s'), char('p_xx'), char('p_xy'), char('p_xy_N'), char('p_xy_NN'), char('p_ve')},...
    {p, p_s, p_xx, p_xy, p_xy_N, p_xy_NN, p_ve});

end