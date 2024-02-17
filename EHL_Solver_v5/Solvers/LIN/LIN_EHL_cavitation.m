function [pressure_LIN,theta_LIN, flowrate_LIN] = ...
    LIN_EHL_cavitation(h, pressure_LIN, theta_LIN, fluid, mesh, simulation_proporties)

%%%WHAT ENTERS HERE SHOULD BE ALREADY NON_DIMENSIONAL PLEASE MAKE YOUR SOLVERS FULLY NON_DIMENSIONAL%%%%%

%Initial condiations 
p = pressure_LIN("p0"); pT = pressure_LIN("p"); pD = pressure_LIN("p1"); 
theta = theta_LIN("theta0"); thetaD = theta_LIN("theta1"); thetaT = theta_LIN("theta");

p(1) = fluid.BC("p_i"); p(end) = fluid.BC("p_o");
pT(1) = fluid.BC("p_i"); pT(end) = fluid.BC("p_o");

%Numeric Properties
N = mesh.N_x;
zeta = h.^3;
W = zeros(1,N);
C = zeros(1,N);
W(2:end) = zeta(1:end-1);
C(1:end) = zeta(1:end);
dx = mesh.dx;
% E = zeros(1,N);
% E(1:end-1) = zeta(2:end);   %unneccesary at the moment.

%%Convergence
convergence = false;
counter = 0;
message = 100;

p_old = zeros(1,N); theta_old = ones(1,N) * -1;  

while ~convergence
    counter = counter + 1;
    
    for i = 2:N-1
        if p(i) > 0 || theta(i) >= 1
            
            %Full backward Difference 
            term_1 = C(i) * theta(i) + W(i) * theta(i-1);
            term_2 = 6 * dx * (h(i) * theta(i) - h(i-1) * theta(i-1));
            term_3 = C(i) * theta(i) * p(i+1) + W(i) * theta(i-1) * p(i-1);
            
            p(i) = -(term_2 - term_3) / (term_1); 
            
            if p(i) >= 0 
                theta(i) = 1;
            else
                p(i) = 0;
            end
        end
        
        if p(i) <= 0 || theta(i) < 1
            
            term_1 = 6 * dx * h(i) - C(i) * (p(i+1)-p(i));
            term_2 = W(i) * theta(i-1) * (p(i) - p(i-1));
            term_3 = 6 * dx * h(i-1) * theta(i-1);
            
            theta(i) =  (-term_2 + term_3) / term_1;
               
             if theta(i) < 1
                 p(i) = 0;
             else
                 theta(i) = 1;
             end
        end
    end
    
    Error = norm(p_old - p);
    p_old = p;

    if Error < simulation_proporties.tollp
        convergence = true;
    end
    
    if counter > simulation_proporties.maxiterR
        error("Max itterration is reached!");
        break
    end

end

p_old = zeros(1,N); theta_old = zeros(1,N);

counter = 0;
convergence = false;

if fluid.De == 0 
    convergence = true;
    pT = p;
    thetaT = theta; 
end
hm = trapz(h.^-2) / trapz(h.^-3);  

while ~convergence
    counter = counter + 1;
    
    for i = 2:N-1
        if pT(i) > 0 || thetaT(i) >= 1
            
            term_1 = C(i) * theta(i) * pD(i+1) + W(i) * theta(i-1) * pD(i-1); 
            
            term_2 = C(i) * theta(i) + W(i) * theta(i-1);
            
            term_3 = C(i) * thetaD(i) * (p(i+1) - p(i)) - W(i) * thetaD(i-1) * (p(i) - p(i-1));
            
            term_4 = 6*dx*(thetaD(i)*h(i) - thetaD(i-1)*h(i-1));
            
            term_5 = 4 * (1-fluid.beta) * (theta(i)*(h(i+1)-h(i)) - theta(i-1)*(h(i)-h(i-1))); 
            term_6 = 18 * (1-fluid.beta) * hm * ( ...
                (theta(i)*(hm/h(i)-1)*(h(i+1)-h(i)) / h(i) ) - ...
                (theta(i-1)*(hm/h(i-1)-1)*(h(i) - h(i-1)) / h(i-1))...
                );
           
            pD(i) =  -(term_4 - term_5 - term_6 - term_3 - term_1) / (term_2);
            pT(i) = p(i) + fluid.De * pD(i);
                
            if pT(i) >= 0 
                thetaT(i) = 1;
                thetaD(i) = (1 - theta(i)) / fluid.De;
            else
                pT(i) = 0;
                pD(i) = -p(i) / fluid.De;
            end
        end
        
        if pT(i) <= 0 || thetaT(i) < 1
            
            term_1 = C(i) * (p(i+1) - p(i)) - 6*dx * h(i);
            
            term_2 = (W(i) * (p(i) - p(i-1)) - 6*dx * h(i-1))*thetaD(i-1);
           
            term_3 = C(i) * theta(i) * (pD(i+1)-pD(i)) - W(i) * theta(i-1) * (pD(i) - pD(i-1));
                  
            term_5 = 4 * (1-fluid.beta) * (theta(i)*(h(i+1)-h(i)) - theta(i-1)*(h(i)-h(i-1)));
            term_6 = 18 * (1-fluid.beta) * hm * ( ...
                (theta(i)*(hm/h(i)-1)*(h(i+1)-h(i)) / h(i) ) - ...
                (theta(i-1)*(hm/h(i-1)-1)*(h(i) - h(i-1)) / h(i-1))...
                );
            
            thetaD(i) =  (term_2 - term_3 - term_5 - term_6) / (term_1);
            thetaT(i) = theta(i) + fluid.De * thetaD(i);
               
             if thetaT(i) < 1
                 pT(i) = 0;
                 pD(i) = -p(i) / fluid.De;
             else
                 thetaT(i) = 1;
                 thetaD(i) = (1 - theta(i)) / fluid.De;
             end
             thetaT(i) = theta(i) + fluid.De * thetaD(i);
        end
    end
    
    Error = norm(p_old-pT) + norm(theta_old - thetaT);
    p_old = pT; theta_old = thetaT;

    if Error < simulation_proporties.tollp
        convergence = true;
    end
    
    if counter > simulation_proporties.maxiterR
        error("Max itteration is reached in non-newtonian part")
        break
    end
  
end

pressure_LIN = containers.Map({char('p'), char('p1'), char('p0')}, {pT, pD, p});
theta_LIN = containers.Map({char('theta'), char('theta1'), char('theta0')}, {thetaT, thetaD, theta});



%Mass Flow Rate
dpdx = OneDcentraldiff(p, mesh.dx);
dpDdx = OneDcentraldiff(pD, mesh.dx);
dhdx = OneDcentraldiff(h, mesh.dx);
flowrate = h/2 - h.^3/12 .* dpdx;  %u
flowrateD = ((1 - fluid.beta) ./ h) .* dhdx .* (1 - 3*hm./h) .* (2 - 3*hm./h) .* (-h/6) - h.^3/12 .* dpDdx; %uD
flowrateT = (flowrateD .* theta + flowrate .* thetaD)*fluid.De + flowrate .* theta;

flowrate_LIN = containers.Map({char('flowrateD'), char('flowrateT'), char('flowrate')}, {flowrateD, flowrateT, flowrate});


end 