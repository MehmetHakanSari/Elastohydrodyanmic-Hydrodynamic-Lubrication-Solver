function [press,flowrateT,flowrate,flowrateD, thetaT,theta,thetaD, pressure0, pressure1] = ...
    LIN_EHL_cavitation(nx, mu, h, pin, poldi, pout, pc, Ux, tollp, L, maxiterR, De, viscos, ...
    pressure0_old, pressure1_old, theta_old, thetaT_old, thetaD_old, h_ref)

%   nx      = number of nodes
%   mu      = viscosity of the lubricant              [use to convert dimensional pressure]
%   h       = heigt profile
%   dhdx    = derivative of height profile  (central difference)         
%   pin     = inlet pressure
%   poldi   = previous pressure 
%   pout    = pressure at outlet
%   pc      = pressure at cavitation 
%   Ux      = velocity profile                        [use to convert dimensional pressure]
%   tollp   = tolerance in pressure
%   maxiter = maxiteration

%non dimensional pressure 
% p  = (P - P_c) / (mu * Uref * L / H_max^2)
% x = x* / L
% h = h* / H_max


N = nx; %Grid number along x
dx = 1 / (N-1);      %non-dimensional step size. 

%surface profile defined by h

% h_max = max(h);
h = h / h_ref;   %non-dimensional height

pT_initial = poldi * h_ref^2 / (mu * Ux * L);   %non-dimensional old pressure
p_initial = pressure0_old * h_ref^2 / (mu * Ux * L); 
pD_inital = pressure1_old * h_ref^2 / (mu * Ux * L); 


%Initial condiations 
% theta = ones(1,N); thetaD = zeros(1,N); thetaT = ones(1,N); 
% p = ones(1,N); pD = zeros(1,N); pT = ones(1,N);   %Okay,Now try now initil condiations. 

%Initial condiations 
p = p_initial; pT = pT_initial; pD = pD_inital; 
theta = theta_old; thetaD = thetaD_old; thetaT = thetaT_old;

p(1) = pin; p(end) = pout;
pT(1) = pin; pT(end) = pout;


%Numeric Properties
hm = trapz(h.^-2) / trapz(h.^-3);  
zeta = h.^3;
W = zeros(1,N);
C = zeros(1,N);
W(2:end) = zeta(1:end-1);
C(1:end) = zeta(1:end);

% E = zeros(1,N);
% E(1:end-1) = zeta(2:end);   %unneccesary at the moment.

%%Convergence
convergence = false;
counter = 0;
message = 100;

p_old = zeros(1,N); theta_old = ones(1,N) * -1;   %is it non-dimensional?

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

    if Error < tollp
        convergence = true;
    end
    
    if counter > maxiterR
        error("Max itterration is reached!");
        break
    end

end

p_old = zeros(1,N); theta_old = zeros(1,N);

counter = 0;
convergence = false;

if De == 0 
    convergence = true;
    pT = p;
    thetaT = theta; 
end

while ~convergence
    counter = counter + 1;
    
    for i = 2:N-1
        if pT(i) > 0 || thetaT(i) >= 1
            
            term_1 = C(i) * theta(i) * pD(i+1) + W(i) * theta(i-1) * pD(i-1); 
            
            term_2 = C(i) * theta(i) + W(i) * theta(i-1);
            
            term_3 = C(i) * thetaD(i) * (p(i+1) - p(i)) - W(i) * thetaD(i-1) * (p(i) - p(i-1));
            
            term_4 = 6*dx*(thetaD(i)*h(i) - thetaD(i-1)*h(i-1));
            
            term_5 = 4 * (1-viscos) * (theta(i)*(h(i+1)-h(i)) - theta(i-1)*(h(i)-h(i-1))); 
            term_6 = 18 * (1-viscos) * hm * ( ...
                (theta(i)*(hm/h(i)-1)*(h(i+1)-h(i)) / h(i) ) - ...
                (theta(i-1)*(hm/h(i-1)-1)*(h(i) - h(i-1)) / h(i-1))...
                );
           
            pD(i) =  -(term_4 - term_5 - term_6 - term_3 - term_1) / (term_2);
            pT(i) = p(i) + De * pD(i);
                
            if pT(i) >= 0 
                thetaT(i) = 1;
                thetaD(i) = (1 - theta(i)) / De;
            else
                pT(i) = 0;
                pD(i) = -p(i) / De;
            end
        end
        
        if pT(i) <= 0 || thetaT(i) < 1
            
            term_1 = C(i) * (p(i+1) - p(i)) - 6*dx * h(i);
            
            term_2 = (W(i) * (p(i) - p(i-1)) - 6*dx * h(i-1))*thetaD(i-1);
           
            term_3 = C(i) * theta(i) * (pD(i+1)-pD(i)) - W(i) * theta(i-1) * (pD(i) - pD(i-1));
                  
            term_5 = 4 * (1-viscos) * (theta(i)*(h(i+1)-h(i)) - theta(i-1)*(h(i)-h(i-1)));
            term_6 = 18 * (1-viscos) * hm * ( ...
                (theta(i)*(hm/h(i)-1)*(h(i+1)-h(i)) / h(i) ) - ...
                (theta(i-1)*(hm/h(i-1)-1)*(h(i) - h(i-1)) / h(i-1))...
                );
            
            thetaD(i) =  (term_2 - term_3 - term_5 - term_6) / (term_1);
            thetaT(i) = theta(i) + De * thetaD(i);
               
             if thetaT(i) < 1
                 pT(i) = 0;
                 pD(i) = -p(i) / De;
             else
                 thetaT(i) = 1;
                 thetaD(i) = (1 - theta(i)) / De;
             end
             thetaT(i) = theta(i) + De * thetaD(i);
        end
    end
    
    Error = norm(p_old-pT) + norm(theta_old - thetaT);
    p_old = pT; theta_old = thetaT;

    if Error < tollp
        convergence = true;
    end
    
    if counter > maxiterR
        error("Max itteration is reached in non-newtonian part")
        break
    end
  
end

% figure(3); hold on
% plot(pT, "--")

press = pT * (mu * Ux * L) / (h_ref^2);  %convert non-dimensional pressure to dimensional
pressure0 = p * (mu * Ux * L) / (h_ref^2);
pressure1 = pD * (mu * Ux * L) / (h_ref^2);

%Mass Flow Rate
dpdx = OneDcentraldiff(p, dx);
dpDdx = OneDcentraldiff(pD, dx);
dhdx = OneDcentraldiff(h, dx);
flowrate = h/2 - h.^3/12 .* dpdx;  %u
flowrateD = ((1 - viscos) ./ h) .* dhdx .* (1 - 3*hm./h) .* (2 - 3*hm./h) .* (-h/6) - h.^3/12 .* dpDdx; %uD
flowrateT = (flowrateD .* theta + flowrate .* thetaD)*De + flowrate .* theta;


end 