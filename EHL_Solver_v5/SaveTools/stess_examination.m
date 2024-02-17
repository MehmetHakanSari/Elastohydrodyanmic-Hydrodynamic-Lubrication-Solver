function stess_examination(solution)

savepath = "Results\Figures\";

if ~isfolder("Results/Figures/EHL_stress_examination")
       mkdir("Results/Figures/EHL_stress_examination")
end

savepath = "Results\Figures\EHL_stress_examination\";


stress_cell = solution.stress_field;
velocity_cell = solution.velocity_field;



h = solution.h / solution.domain.href;

T_xx = stress_cell{1, 1}; TD_xx = stress_cell{2, 1};
T_xy = stress_cell{1, 2}; TD_xy = stress_cell{2, 2};
u = velocity_cell{1,1}; uD = velocity_cell{2,1}; 

[N, J] = size(T_xx);

dy = h / (J - 1); % dy spaceing for each location of x.
x = linspace(0, 1, solution.domain.N_x);
dx = x(2) - x(1);
De = solution.deborah_Number;

gradTxx = TwoDcentraldiff(T_xx + De * TD_xx, dx , dy);
gradTxy = TwoDcentraldiff(T_xy + De * TD_xy, dx, dy);
gradu = TwoDcentraldiff(u  + De * uD, dx, dy);
grad2u = TwoDcentraldiff(gradu{2}, dx, dy); 





INT_xx = dy .* trapz(T_xx);
INT_Dxx =  dy .* trapz(De * TD_xx); 

INT_xy = dy .*  trapz(T_xy);
INT_Dxy = dy .*  trapz(De * TD_xy);

INT_all = dy .* trapz(grad2u{2}* solution.viscocity_ratio + gradTxx{1} + gradTxy{2});
INT_u = dy .* trapz(grad2u{2}* solution.viscocity_ratio);
INT_dxx = dy .* trapz(gradTxx{1});
INT_dxy = dy .* trapz(gradTxy{2});

figure(61); hold on; box on;
% surf(grad2u{2}* solution.viscocity_ratio + gradTxx{1} + gradTxy{2}, "linestyle", "none");
plot(x, INT_all)
plot(x, INT_u)
plot(x, INT_dxx)
plot(x, INT_dxy)

lw = 1.9; ms = 20;

figure(1); hold on; box on;
lh = plot(x, INT_xx, "-", "LineWidth", lw);
lh = plot(x, INT_Dxx, "--", "LineWidth", lw);
xlabel("x")
ylabel("T_xx int")
legend("N", "NN","box", "off")

figure(2); hold on; box on;
lh = plot(x, INT_xy, "-", "LineWidth", lw);
lh = plot(x, INT_Dxy, "--", "LineWidth", lw);
xlabel("x")
ylabel("T_xy int")
legend("N", "NN","box", "off")

figure(3); hold on; box on;
lh = plot(x, INT_xx, "-", "LineWidth", lw);
lh = plot(x, INT_Dxx, "--", "LineWidth", lw);
lh = plot(x, INT_xy, "-", "LineWidth", lw);
lh = plot(x, INT_Dxy, "--", "LineWidth", lw);
xlabel("x")
ylabel("T_xx int")
legend("N_xx", "NN_xx", "N_xy", "NN_xy","box", "off")


end