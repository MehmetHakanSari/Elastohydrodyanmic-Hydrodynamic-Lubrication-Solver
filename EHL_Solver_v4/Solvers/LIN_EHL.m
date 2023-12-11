function [pressure, p, pD] =  LIN_EHL(nx, mu, h, pin, poldi, pout, pc, Ux, tollp, L, maxiterR, De, viscos, h_ref)



tic
N = nx;

x = linspace(0, 1, N);
dx = 1 / (N-1);   
h = h / h_ref;
dhdx = OneDcentraldiff(h, dx);
d2hdx2 = OneDcentraldiff(dhdx, dx);

%first order solution to pressure
hm = trapz(x,h.^-2) / trapz(x,h.^-3);               %profile feature
dpdx = 6 *   (h - hm) ./ h.^3;                        %pressure derivative in l order
% d2pdx2 = -6 * dhdx .*( -3 * hm./h.^4 + 2 * 1./h.^3);
p = cumtrapz(x, dpdx);                              %solution to pressure w.r.t x in l    

%Li's direct solution for Deborah Order

h0 = h(1);
h1 = h(end);
hd = ((3/8 * (hm^2/h1^4 - hm^2/h0^4)) + (-1/2 * (hm/h1^3 - hm/h0^3)) + (1/6 * (1/h1^2 - 1/h0^2))) * trapz(x,h.^-3)^-1;
pD = (1 - viscos) * (9/2 * (hm^2./h.^4 - hm^2/h0^4) - 6 * (hm./h.^3 - hm/h0^3) + 2 * (1./h.^2 - 1/h0^2) - 12 * hd * cumtrapz(x,h.^-3));
%Li's pressure derivative solution for Deborah Order

dpDdx = -18 * (hm^2./h.^5 .* dhdx) + 18 * (hm./h.^4 .* dhdx) - 4 .* (1./h.^3) .* dhdx - 12 * hd./h.^3;
pDe = cumtrapz(x,dpDdx);


% %Pressure Derivatives and Pressure Distribution
% p_MAT = zeros(J,N);
% for i = 1:J
%     p_MAT(i,:) = p;          %Pressure Distribution variable name stands for p_MATRİX
% end
% grad_p = TwoDcentraldiff(p_MAT, dx, dy);
% grad2_p = TwoDcentraldiff(grad_p{1}, dx, dy);
% error_CDvsAN = (grad2_p{1}(1,:) - d2pdx2) ./ d2pdx2 * 100;


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

p_TDMA = zeros(1,N);
p_TDMA(2: end-1) = TDMA(W(2:end-1),C(2:end-1),E(2:end-1),F(2:end-1));

% figure(1)
% plot(p_TDMA * De)

%TDMA results can be different than direct formulation of p for low mesh. I
%prefer TDMA for smaller mesh. 
p_check = zeros(1,N);
p_check(2:end-1) = TDMA(W(2:end-1),C(2:end-1),E(2:end-1),CHECK(2:end-1)); 

% figure(2)
% plot(p_check)


disp("Linearized Reynolds Equation is solved")
disp("-----------------------------------------")
toc

% p(0) = pin; p(end) = pout;


% figure(12); hold on
% plot(x, p)
% plot(x,pD)
% 
% e = max(h) - min(h);
% ho = min(h);
% dpdxN = (1 / ho^2) * (e / (2 * ho + 3 * e) );
% dpdxNN_i = (4 * ho + 42 * e) / (12 * ho + 216 * e);
% dpdxNN_o = (8 * ho + 66 * e) / (24 * ho + 216 * e);
% 
% 
% pN_anal_inlet = dpdxN * x(1:round(length(x)/2));
% pN_anal_outlet = 2*max(pN_anal_inlet) + (-dpdxN) * x(round(length(x)/2)+1: end);
% pN_anal = [pN_anal_inlet pN_anal_outlet];
% 
% pNN_anal_inlet = -dpdxNN_i * x(1:round(length(x)/2));
% pNN_anal_outlet = 2*max(abs(pNN_anal_inlet)) + (-dpdxNN_o) * x(round(length(x)/2)+1: end);
% 
% pNN_anal = [pNN_anal_inlet pNN_anal_outlet-min(pNN_anal_outlet)];
% 
% figure(13); hold on
% plot(x, 6*pN_anal)
% plot(x, pNN_anal)
% 
% figure(14);hold on
% plot(x, p, "--")
% plot(x, 6*pN_anal, "-")
% plot(x, pD, "--")
% plot(x, pNN_anal, "-")
% 
% figure(15); hold on
% plot(x, -p ./ pD)
% 
% dpDdx_num = OneDcentraldiff(pD, x(2) - x(1));
% display("dpdx_i / dpdx_o analitical: " + string((min(h) / max(h))^3))
% display("dpdx_i / dpdx_o numerical: " + string(dpDdx_num(15)/dpDdx_num(end-15)))
% display("hm analitical: " + string(ho * (ho + 2 * e) / (ho + 3/2*e)))
% display("hm numerical: " + string(hm))
% display("step term numericcal: " + string(9/2 * (hm^2./h1^4 - hm^2/h0^4) - 6 * (hm./h1^3 - hm/h1^3) + 2 * (1./h1^2 - 1/h0^2)))
% display("step term analitical: " + string((4*h1 + 21 * e) / (12 * h1 + 90*e)))
% display("dpdx_i : " + string(dpDdx_num(15)))
% display("dpdx_o : " + string(dpDdx_num(end-15)))

% figure(1); hold on 
% plot(linspace(0,1,N),De * pD,"-","LineWidth", 2)
p = p * (mu * Ux * L) / (h_ref^2);
pD = pD ;%* (mu * Ux * L) / (h_ref^2);

% figure(3)
% plot(p)
% 
% figure(4)
% plot(pD * De)




pressure = p + pD * De;
end