function [K,Dlength,npuntiD,GFcut] = ElasticKernel(E,nu1,L,N,dx,magn)

gf_factor = 10; %2048
N_gf = N * gf_factor;
dx_gf = L / N_gf;
x = linspace(-L/2, L/2, N_gf);

q0 = 2 * pi / L; %fundemental wavelength

displacement_vector = zeros(1, N_gf);

k = [1:(N_gf / 2)-1  flip(-(1:(N_gf/2)))];
% k = [flip(1:(N_gf / 2)-1)  (-(1:(N_gf/2)))];

displacement_vector(2:end) = ...
        ( (-1).^k ) .* (q0 / (2*pi) ) .* ...
        ( 2 * sin( (1/2) .* k * q0 * dx) ./ (k .* q0) ) ...
        .*(-2 * (1 - nu1^2)*(1 ./ abs(k * q0) ));


plot(x, displacement_vector, "k-", "LineWidth", 1.8)
legend("\boldmath{$\chi_{\lambda}(x)$}", "box", "off", "Fontsize", 20, "Interpreter","latex")
xlabel("Magnified Length","Interpreter","latex"); set(gca, "Fontsize", 16); set(gca, 'TickLabelInterpreter', 'latex');
ax2 = axes('XAxisLocation','top', ...
    'Color','none', ...
    'XColor','k');
ax2.YAxis.Visible = 'off';
ax2.XLim = [1 N_gf];
xlabel(ax2, "Node Number","Interpreter","latex")
set(gca, "Fontsize", 16); set(gca, 'TickLabelInterpreter', 'latex');

Q = fft(displacement_vector);
Q = [Q Q(1)] / E;

figure(2)
plot(real(Q), "k-", "LineWidth", 1.8)
legend("\boldmath{$\chi_{\lambda}(q)$}", "box", "off", "Fontsize", 20, "Interpreter","latex") 
xlim([1 N_gf+1])
xlabel("\bf{q}","Interpreter","latex");  set(gca, "Fontsize", 16); set(gca, 'TickLabelInterpreter', 'latex');
ax2 = axes('XAxisLocation','top', ...
    'Color','none', ...
    'XColor','k');
ax2.YAxis.Visible = 'off';
ax2.XLim = [1 N_gf];
xlabel(ax2, "Node Number","Interpreter","latex")
set(gca, "Fontsize", 16); set(gca, 'TickLabelInterpreter', 'latex');



x_green_function = zeros(1, N+1);
green_function = zeros(1, N+1);

x_green_function(1:end-1) = x((0:N-1) * gf_factor + 1); 
x_green_function(end) = x(end); 

green_function(1:end-1) = Q((0:N-1) * gf_factor + 1);
green_function(end) = Q(end); 

figure(3)
plot(x_green_function, real(green_function), ".-", "Linewidth", 1.5, "Markersize", 16)
xlabel("Magnified Length","Interpreter","latex");  



% green_function
Dlength = L / magn;
% Rlength = 1 / magn;
N_d = round(N / magn); %real node number 

GFcut = real( green_function( round((N+1) / 2) : round((N+1) / 2) + N_d-1 ) );

figure(3); hold on
plot(linspace(-L/2, L/2, N_d), real(GFcut), "*-", "Linewidth", 1.5, "Markersize", 16)
legend("\boldmath{$G_{\lambda}(q)$}", "\boldmath{$G_{\lambda}(cut)$}", "box", "off", "Fontsize", 20, "Interpreter","latex")
set(gca, "Fontsize", 16); set(gca, 'TickLabelInterpreter', 'latex');
ax2 = axes('XAxisLocation','top', ...
    'Color','none', ...
    'XColor','k');
ax2.YAxis.Visible = 'off';
ax2.XLim = [1 N + 1];
xlabel(ax2, "Node Number","Interpreter","latex")
set(gca, "Fontsize", 16); set(gca, 'TickLabelInterpreter', 'latex');

K = toeplitz(GFcut);
end