function [p, p_s, p_xx, p_xy, p_xy_N, p_xy_NN, p_ve] = PressureSolver_v1(coeffs)

DTxxDx = coeffs{1};
x = coeffs{2};
dy = coeffs{3};
y_MAT = coeffs{4};
h = coeffs{5};
dx = coeffs{6};
dhdx = coeffs{7};
T_xy = coeffs{8};
viscos = coeffs{9};
%-----------------------------------

[m, n] = size(DTxxDx);
dydy = zeros(m, n);
N = n;

for i = 1:m
    dydy(i,:) = dy;
end

TD = y_MAT ./ h; % Transformed Grid 

zeta = h.^3 / 12;
W = zeros(1, N);
E = zeros(1, N);
W(2:end) = (zeta(1:end-1) + zeta(2:end)) / (2*dx^2);
E(1:end-1) = (zeta(1:end-1) + zeta(2:end)) / (2*dx^2);
Cent = -(E + W);

Q1 = 0.5 * (1-viscos) * dhdx;   %Newtonian stress solvent contrubution
Q4 = 0.5 * viscos * dhdx;      %Newtonian stress polymer contrubution

I11 = dy .* trapz(TD(:,1) * (dy .* trapz(dydy .*cumtrapz(flip(DTxxDx))))) - dy .* trapz(dydy .* cumtrapz(dydy .* cumtrapz(flip(DTxxDx))));
I12 = dy .* (trapz(TD(:,1) * (dy .* trapz(flip(T_xy)))) -  trapz(dydy .* cumtrapz((flip(T_xy)))));
Q2 = OneDcentraldiff(I11,dx);             % xx component of the pressure
Q3 = OneDcentraldiff(I12,dx);             % xy component of the pressure

p_s = zeros(1,N);
p_xx = zeros(1,N);
p_xy = zeros(1,N);
p_xy_N = zeros(1,N);

p_s(2:end-1) = TDMA(W(2:end-1),Cent(2:end-1),E(2:end-1),Q1(2:end-1));
p_xx(2:end-1) = TDMA(W(2:end-1),Cent(2:end-1),E(2:end-1),Q2(2:end-1));
p_xy(2:end-1) = TDMA(W(2:end-1),Cent(2:end-1),E(2:end-1),Q3(2:end-1));
p_xy_N(2:end-1) = TDMA(W(2:end-1),Cent(2:end-1),E(2:end-1),Q4(2:end-1));

p_xy_NN = p_xy - p_xy_N;         %Non-Newtonian pressure in xy component
p_ve = p_xx + p_xy_NN;           %Total Non-newtonian pressure

%Load
p = p_s + p_xx + p_xy;


end