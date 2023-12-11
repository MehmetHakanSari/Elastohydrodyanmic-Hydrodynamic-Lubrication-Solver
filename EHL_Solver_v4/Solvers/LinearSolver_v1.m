function [u, v, DuDx, dudy, DvDx, dvdy, T_xx, T_xy, T_yy, p_LIN, pD_LIN] = LinearSolver_v1(coeffs)

% I check the linear solver with the previous ones. It seems okay in
% pressurewise. If it is same in pressurewise it is also same in
% velocitywise. However, calculated EHL height is different. This results
% different pressure outcome in error of < 10% 

x = coeffs{1};
h = coeffs{2};
dhdx = coeffs{3};
dx = coeffs{4};
dy = coeffs{5};
% d2hdx2 = coeffs{6};
y_MAT = coeffs{7};
dydx = coeffs{8};
viscos = coeffs{9};

J = length(y_MAT(:,1));
N = length(y_MAT(1,:));

%first order solution to pressure
hm = trapz(x,h.^-2) / trapz(x,h.^-3);               %profile feature
dpdx = 6 *   (h - hm) ./ h.^3;                        %pressure derivative in l order
p_LIN = cumtrapz(x,dpdx); 

h0 = h(1);
h1 = h(end);
hd = ((3/8 * (hm^2/h1^4 - hm^2/h0^4)) + (-1/2 * (hm/h1^3 - hm/h0^3)) + (1/6 * (1/h1^2 - 1/h0^2))) * trapz(x,h.^-3)^-1;
pD_LIN = 9/2 * (hm^2./h.^4 - hm^2/h0^4) - 6 * (hm./h.^3 - hm/h0^3) + 2 * (1./h.^2 - 1/h0^2) - 12 * hd * cumtrapz(x,h.^-3);

%First order velocity Solutions
u = 1/ 2 * dpdx .* (y_MAT.^2 - y_MAT .* h) + (1 - y_MAT./h);
v = dhdx .* (2 - 3 * hm./h) .* (y_MAT.^3 ./ h.^3 - y_MAT.^2 ./ h.^2);

dudy = dpdx .* (y_MAT - h / 2) - 1 ./ h; %This comes from zero order pressure derivative
dvdy = dhdx .* (2 - 3 * hm./h) .* (3 * y_MAT.^2 ./ h.^3 - 2 * y_MAT ./ h.^2);

gradu = TwoDcentraldiff(u, dx, dy);
gradv = TwoDcentraldiff(v, dx, dy);
%Derivative w.r.t new coordinates
DuDx = (gradu{1} + gradu{2}.* dydx); %gradu{2} == dudy (I checked)
DvDx = (gradv{1} + gradv{2} .* dydx);

%First order stress
T_xx = zeros(J,N);
T_xy = viscos .* dudy;
T_yy = 2 * viscos .* dvdy;

end