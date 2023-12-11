function [u, v, DuDx, dudy, DvDx, dvdy, T_xx, T_xy, T_yy,p, pD] = LinearSolver(coeffs)


tic
x = coeffs{1};
h = coeffs{2};
dhdx = coeffs{3};
dx = coeffs{4};
dy = coeffs{5};
d2hdx2 = coeffs{6};
y_MAT = coeffs{7};
dydx = coeffs{8};
viscos = coeffs{9};
path =  "E:\Bilkent Dökümanları\Masterımsı\Kodlar\Data\LIN\";


J = length(y_MAT(:,1));
N = length(y_MAT(1,:));

%first order solution to pressure
hm = trapz(x,h.^-2) / trapz(x,h.^-3);               %profile feature
dpdx = 6 *   (h - hm) ./ h.^3;                        %pressure derivative in l order
d2pdx2 = -6 * dhdx .*( -3 * hm./h.^4 + 2 * 1./h.^3);
p = cumtrapz(x,dpdx);                              %solution to pressure w.r.t x in l order

%Li's direct solution for Deborah Order

h0 = h(1);
h1 = h(end);
hd = ((3/8 * (hm^2/h1^4 - hm^2/h0^4)) + (-1/2 * (hm/h1^3 - hm/h0^3)) + (1/6 * (1/h1^2 - 1/h0^2))) * trapz(x,h.^-3)^-1;
QD = hd;    %Volumetric Flow Rate of the Deborah Order
pD = 9/2 * (hm^2./h.^4 - hm^2/h0^4) - 6 * (hm./h.^3 - hm/h0^3) + 2 * (1./h.^2 - 1/h0^2) - 12 * hd * cumtrapz(x,h.^-3);

%Li's pressure derivative solution for Deborah Order

dpDdx = -18 * (hm^2./h.^5 .* dhdx) + 18 * (hm./h.^4 .* dhdx) - 4 .* (1./h.^3) .* dhdx - 12 * hd./h.^3;
pDe = cumtrapz(x,dpDdx);

%Pressure Derivatives and Pressure Distribution
p_MAT = zeros(J,N);
for i = 1:J
    p_MAT(i,:) = p;          %Pressure Distribution variable name stands for p_MATRİX
end
grad_p = TwoDcentraldiff(p_MAT, dx, dy);
grad2_p = TwoDcentraldiff(grad_p{1}, dx, dy);
error_CDvsAN = (grad2_p{1}(1,:) - d2pdx2) ./ d2pdx2 * 100;


%TDMA Coefficients
zeta = h.^3 / 12;
W = zeros(1,N);
E = zeros(1,N);
W(2:end) = (zeta(1:end-1) + zeta(2:end)) / (2*dx^2);
E(1:end-1) = (zeta(1:end-1) + zeta(2:end)) / (2*dx^2);
C = -(E + W);
F1 = d2hdx2 .* (-1 / 3 + 1.5 * hm ./ h - 1.5 * (hm ./ h).^2);   %First order pressure. 
F2 = 1 ./ h .* dhdx.^2 .* (-1.5 * hm ./ h + 3 * (hm ./ h).^2);
F = F1 + F2;
CHECK = 0.5 * dhdx;

p_1 = zeros(1,N);
p_1(2: end-1) = TDMA(W(2:end-1),C(2:end-1),E(2:end-1),F(2:end-1));

%TDMA results can be different than direct formulation of p for low mesh. I
%prefer TDMA for smaller mesh. 
p_check = zeros(1,N);
p_check(2:end-1) = TDMA(W(2:end-1),C(2:end-1),E(2:end-1),CHECK(2:end-1)); 

%First order velocity Solutions

u = 1/ 2 * dpdx .* (y_MAT.^2 - y_MAT .* h) + (1 - y_MAT./h);
v = dhdx .* (2 - 3 * hm./h) .* (y_MAT.^3 ./ h.^3 - y_MAT.^2 ./ h.^2);
Q = h / 2 - 1 / 12 * dpdx .* h.^3;  %Volumetric Flow Rate Analitically
Q_NUM = trapz(u);                   %Volumetric Flow Rate Numerically

dudy = dpdx .* (y_MAT - h / 2) - 1 ./ h; %This comes from zero order pressure derivative
dvdy = dhdx .* (2 - 3 * hm./h) .* (3 * y_MAT.^2 ./ h.^3 - 2 * y_MAT ./ h.^2);

dudx = 1/2 * d2pdx2 .* (y_MAT.^2 - y_MAT .* h) - 1/2 * dpdx .* y_MAT .* dhdx + y_MAT .* h.^(-2) .* dhdx;

gradu = TwoDcentraldiff(u, dx, dy);
gradv = TwoDcentraldiff(v, dx, dy);
%Derivative w.r.t new coordinates
DuDx = (gradu{1} + gradu{2}.* dydx); %gradu{2} == dudy (I checked)
DvDx = (gradv{1} + gradv{2} .* dydx);

%Validation of DvDx
%     gradDu = TwoDcentraldiff(DuDx, dx ,dy);
%     gradDv = TwoDcentraldiff(DvDx, dx, dy);
%     D2uDx2 = (gradDu{1} - gradDu{2}.* y_MAT ./ h .* dhdx);
%     D2vDxDy = gradDv{2};

%First order stress
T_xx = zeros(J,N);
T_xy = viscos .* dudy;
T_yy = 2 * viscos .* dvdy;

%Derivatives of Stress Components (First Order)
gradT_xy = TwoDcentraldiff(T_xy, dx, dy);
gradT_yy = TwoDcentraldiff(T_yy, dx, dy);

DTxyDx = gradT_xy{1} + gradT_xy{2} .* dydx;
DTyyDx = gradT_yy{1} + gradT_yy{2} .* dydx;
dTxydy = gradT_xy{2};
dTyydy = gradT_yy{2};

%Analitical Solution:
%     dudydx = d2pdx2 .* y_MAT - (1 / 2) * (d2pdx2 .* h + dpdx .* dhdx) + dhdx ./h.^2;
%     grad2u = TwoDcentraldiff(gradu{2}, dx, dy);
%     dT_xy_ANdx = (1 - viscos) * dudydx;
%     dT_xy_CDdx = (1 - viscos) * grad2u{1};

% Deborah Order Velocity Solutions
%Li's solution
uD = h.^2 / 2 .* dpDdx .* (y_MAT.^2 ./ h.^2 - y_MAT ./ h);
uD = uD + 1 ./ h .* dhdx .* (1 - 3 * hm./h) .* (2 - 3 * hm./h) .* (y_MAT.^2 ./ h.^2 - y_MAT ./ h);

graduD = TwoDcentraldiff(uD, dx, dy);
duDdx = graduD{1};
duDdy = graduD{2};
DuDDx = duDdx + duDdy .* dydx;

%By continuity duDdx + dvDdy = 0
%Taking derivative then integrating over y would give me vD.
dydy = reshape(repelem(dy,J),J,N);
vD = dydy .* cumtrapz(DuDDx);
gradvD = TwoDcentraldiff(vD, dx, dy);

%Deborah order stresss
TD_xx = 2 * T_xy .* dudy;
TD_xy = -(u .* DTxyDx+ v .* dTxydy - T_yy .* dudy) + viscos * graduD{2};
TD_yy = -(u .* DTyyDx + v .* dTyydy - 2 * T_xy .* DvDx - 2 * T_yy .* dvdy) + 2 * viscos * gradvD{2};

%     plot(x, h, "-.", "linewidth",1.5, "MarkerSize",1.5 );
%     plot(x,pDe);

LD = -trapz(pD, x);
LN = -trapz(p, x);

ValueCell = {p_check, pD, T_xx, T_yy, T_xy,LD,LN};

save(path + "LIN_variables" + "_" + string(J) + "x" + string(N), "ValueCell")
disp("Linearized Reynolds Equation is solved")
disp("-----------------------------------------")
toc
end