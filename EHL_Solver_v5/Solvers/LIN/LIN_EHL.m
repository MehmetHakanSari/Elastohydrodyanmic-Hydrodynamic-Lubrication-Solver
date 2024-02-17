function [p_dic_LIN] =  LIN_EHL(h, fluid, mesh)

tic

dhdx = OneDcentraldiff(h, mesh.dx, "CD2");
d2hdx2 = OneDcentraldiff(dhdx, mesh.dx, "CD2");

%Li's direct solution for Deborah Order
h_i = h(1);
h_o = h(end);
hm = trapz(h.^-2) / trapz(h.^-3);  
hd =  (1 - fluid.beta) * ( (3/8 * (hm^2/h_o^4 - hm^2/h_i^4)) + (-1/2 * (hm/h_o^3 - hm/h_i^3)) + (1/6 * (1/h_o^2 - 1/h_i^2)) ) * trapz(mesh.x, h.^-3)^-1;
% pD = (1 - fluid.beta) * ((9/2 * (hm^2./h.^4 - hm^2/h0^4) - 6 * (hm./h.^3 - hm/h0^3) + 2 * (1./h.^2 - 1/h0^2)) - 12 * hd * cumtrapz(mesh.x, h.^-3));
%Li's pressure derivative solution for Deborah Order

% dpDdx =   (-18 * (hm^2./h.^5) + 18 * (hm./h.^4) - 4 * (1./h.^3) ) .* fluid.beta * (1 - viscos)  - 12 * hd./h.^3;
% pDe = cumtrapz(x, dpDdx);

%TDMA Coefficients
zeta = h.^3 / 12;
W = zeros(1, mesh.N_x);
E = zeros(1, mesh.N_x);
W(2:end) = (zeta(1:end-1) + zeta(2:end)) / (2 * mesh.dx^2);
E(1:end-1) = (zeta(1:end-1) + zeta(2:end)) / (2 * mesh.dx^2);
C = -(E + W);

F1 = d2hdx2 .* (-1 / 3 + 1.5 * hm ./ h - 1.5 * (hm ./ h).^2);   %First order pressure. 
F2 = 1 ./ h .* dhdx.^2 .* (-1.5 * hm ./ h + 3 * (hm ./ h).^2);
F = (1 - fluid.beta) * (F1 + F2);
CHECK = 0.5 * dhdx;

p_TDMA = zeros(1, mesh.N_x);
p_TDMA(2: end-1) = TDMA(W(2:end-1),C(2:end-1),E(2:end-1),F(2:end-1));

p_check = zeros(1, mesh.N_x);
p_check(2:end-1) = TDMA(W(2:end-1),C(2:end-1),E(2:end-1),CHECK(2:end-1)); 

dpdx = OneDcentraldiff(p_check, mesh.dx);

q_D_step =  4 * ((1 - fluid.beta) / 48 * (dpdx(5)^2 * (h_o^2 - h_i^2) / 2 + 2 * (1 / h_o^2 - 1 / h_i^2) ) * (mesh.dx * trapz( 1 ./ h.^3) )^(-1));
p_D_step =  ((1 - fluid.beta) / 1 * (dpdx.^2 .* (h.^2 - h_i^2) / 2 + 2 * (1 ./ h.^2 - 1 / h_i^2) ) -  12 * q_D_step * (mesh.dx * cumtrapz( 1 ./ h.^3) ));
% 
% plot(mesh.x,  p_D_step, "LineWidth", 1.8); hold on;
% plot(mesh.x, p_TDMA, "--" ,"LineWidth", 1.8)
disp("-----------------------------------------")
disp("analitical q_D:  " + string(q_D_step))
disp("analitical h_m q_D:  " + string(hd))
disp("-----------------------------------------")
disp("Linearized Reynolds Equation is solved")
disp("-----------------------------------------")
toc

p = p_check; 
pD = p_TDMA; 
pressure = p + pD * fluid.De;

p_dic_LIN = containers.Map({char('p'), char('p1'), char('p0')}, {pressure, pD, p});


end