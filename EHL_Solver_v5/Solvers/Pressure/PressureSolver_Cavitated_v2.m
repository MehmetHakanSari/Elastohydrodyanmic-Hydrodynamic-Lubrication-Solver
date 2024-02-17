function [pressure_dic,  theta] = PressureSolver_Cavitated_v2(h, mesh, fluid, DTxxDx, T_xy)


[M, N] = size(DTxxDx);
dydy = zeros(M, N);

dy = mesh.dy;
dx = mesh.dx;

for i = 1:M
    dydy(i,:) = dy;
end

TD = mesh.Y ./ h; % Transformed Grid 

zeta = h.^3;
W = zeros(1, N);
% E = zeros(1,N); 
W(2:end) = zeta(1:end-1);
% E(1:end-1) = zeta(2:end);
C = zeta(1:end);

I11 = dy .* trapz(TD(:,1) * (dy .* trapz(dydy .*cumtrapz(flip(DTxxDx))))) - dy .* trapz(dydy .* cumtrapz(dydy .* cumtrapz(flip(DTxxDx))));  % xx component of the pressure
I12 = dy .* (trapz(TD(:,1) * (dy .* trapz(flip(T_xy)))) -  trapz(dydy .* cumtrapz((flip(T_xy))))); % xy component of the pressure

theta = ones(1,N);
p = ones(1,N); 
p(1) = fluid.BC("p_i"); p(end) = fluid.BC("p_o");
p_xx = ones(1,N);  p_xy = ones(1,N);  p_xx(1) = 0; p_xx(end) = 0; p_xy(1) = 0; p_xy(end) = 0;
p_old = zeros(1,N); theta_old = zeros(1,N);
p_c =  fluid.BC("p_c"); 


%%Convergence
SS = 1e-12;
convergence = false;
counter = 0;
message = 5000;

while ~convergence
    counter = counter + 1;
    
    for i = 2:N-1
        if p(i) > 0 || theta(i) >= 1
%           Full backward Difference 
            term_1 = C(i) * theta(i) + W(i) * theta(i-1);
            term_2 = 6 * dx * fluid.beta * (h(i) * theta(i) - h(i-1) * theta(i-1)); 
            term_3 = C(i) * theta(i) * p(i+1) + W(i) * theta(i-1) * p(i-1);
            term_4 = 12 * dx * (theta(i) * (I11(i) + I12(i)) - theta(i-1) * (I11(i-1) + I12(i-1)));
            
            p(i) = (term_3 - term_2 - term_4) / (term_1); 
                       
            if p(i) >= 0 
                theta(i) = 1;
            else
                p(i) = p_c;
            end
            
        end
        
        if p(i) <= 0 || theta(i) < 1   
           %Full backward Difference  
            term_1 = C(i) * (p(i+1)-p(i)) - 6 * fluid.beta * dx * h(i) - 12 * dx * (I11(i) + I12(i));
            term_2 = W(i) * theta(i-1) * (p(i) - p(i-1));
            term_3 = 6 * fluid.beta * dx * h(i-1) * theta(i-1);
            term_4 = 12 * dx * theta(i-1) * (I11(i-1) + I12(i-1));

            theta(i) =  (term_2 - term_3 - term_4) / term_1;
                     
             if theta(i) < 1
                 p(i) = p_c;
             else
                 theta(i) = 1;
             end     
        end
    end
    
    error = norm(p_old-p) + norm(theta_old-theta);
    p_old = p;
    theta_old = theta;
   
    if error < SS
        convergence = true;
    end

%     if mod(counter,message) == 0
%         disp("Error: " + string(error))
%         plot(mesh.x, p,"LineWidth",1.9)
%         pause(0.1)
%     end
        
end

pressure_dic = containers.Map({char("p")}, {p});

end