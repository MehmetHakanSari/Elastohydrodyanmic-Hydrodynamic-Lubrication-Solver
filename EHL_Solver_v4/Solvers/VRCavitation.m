function [press, theta, stress_new, p_field] = VRCavitation(dt, N_x, N_y, Spatial_Scheme, Temporal_Scheme, h, mu, p_in, poldi, p_cavitation, ...
                                            Ux, tollp, L, maxiterR, De, viscos, stress_field, velocity_field, p_LIN, pD_LIN, iternumber, h_ref, cavitation_flag)

%   nx      = number of nodes               double
%   ny      = number of nodes               double
%   mu      = viscosity of the lubricant    double   [use to convert dimensional pressure]
%   h       =  heigt profile                 [1 x nx]
%   dhdx    = derivative of height profile  [1 x nx]        
%   p_in     = inlet pressure                double
%   poldi   = previous pressure             [1 x nx] 
%   pout    = pressure at outlet            double
%   p_cavitation = pressure at cavitation        double
%   Ux      = velocity profile              double          [use to convert dimensional pressure]
%   tollp   = tolerance in pressure         double
%   maxiterR = maxiteration                 double
%   viscos                                           double                
%   stress_field                                  cell {2x3}  T_xx,T_xy, T_
%   velocity_field                                cell {4x2}, u,v ; uD, vD; dudx, dvdx; dudy, dvdy  
%   href = reference height                double

%   itterationnumber = itteration numner of outer loop. 

%non dimensional pressure 
% p  = (P - P_c) / (fluid_viscosity * Uref * Total_length / H_H_reference^2)
% x = x* / Total_length
% h = h* / H_reference
p_field = {};


%Non-dimensionalizng
h = h / h_ref;   %non-dimensional height


% p_initial = poldi' * h_ref^2 / (mu * Ux * L);   %non-dimensional old pressure

%Grid Properties
N = N_x; 
J = N_y;

if length(h) ~= N
    h = interp1(linspace(0,1,length(h)), h, linspace(0,1,N));
end

x = linspace(0, 1, N);
dx = 1 / (N - 1);

X = zeros(J, N);
Y = zeros(J, N);

dhdx = OneDcentraldiff(h,dx);
d2hdx2 = OneDcentraldiff(dhdx,dx);

viscos = (1-viscos);   % viscos = 1:full poylmer. viscos = 0: full solvent. viscos is reverse in main code. (1 - viscos) is needed.  

for i = 1:J
    X(i,:) = x;
end
dy = h / (J - 1);
for i = 1:J
    Y(i,:) = h - dy .* (i-1);
end

% grady = TwoDcentraldiff(y_MAT,dx,dy);
dydx = -Y ./ h .* dhdx;             %Coordinate Transfrom

%Setting 
SSboundary = 10^(-7 - round(log10(dt^-1)));         %Steady-State Bound
Error_Boundary = 1000000;

% disp("-----------------------------------------")

% [u, v, DuDx, dudy, DvDx, dvdy, T_xx, T_xy, T_yy, p_LIN, pD_LIN] = ...
%     LinearSolver_v1({x, h, dhdx, dx, dy, d2hdx2, Y, dydx, viscos});  

u = velocity_field{1,1}; v = velocity_field{1,2};
DuDx = velocity_field{3,1}; DvDx = velocity_field{3,2};
dudy = velocity_field{4,1}; dvdy = velocity_field{4,2};
T_xx = stress_field{1,1}; T_xy = stress_field{1,2}; T_yy = stress_field{1,3};

coeffs_stress = {dx, De, T_xx, T_xy, T_yy, dy, u, v, DuDx, DvDx, dudy, dvdy,dydx,...
    Spatial_Scheme, Error_Boundary, SSboundary,dt, viscos, X, Y};

if Temporal_Scheme == "RK4"
    [T_xx, T_xy, T_yy, DTxxDx] = ...
        StressSolver_RK4(coeffs_stress); 
elseif Temporal_Scheme == "IM_EU"
    [T_xx, T_xy, T_yy, DTxxDx] = ...
        StressSolver_IE_EU(coeffs_stress); 
elseif Temporal_Scheme == "EX_EU"
    [T_xx, T_xy, T_yy, DTxxDx] = ...
        StressSolver_EX_EU(coeffs_stress); 
end
    
coeffs_pressure = {DTxxDx, x, dy, Y, h, dx, dhdx, T_xy, viscos, p_in, p_cavitation};

if cavitation_flag == "on" 
    [p, theta] = ...
        PressureSolver_Cavitated_v2(coeffs_pressure);
elseif cavitation_flag == "off"
    [p, p_s, p_xx, p_xy, p_xy_N, p_xy_NN, p_ve] = ...
        PressureSolver_v1(coeffs_pressure);
%     p = p_s + p_xx + p_xy;
    theta = ones(1, N);
    p_field = {p, p_s, p_xx, p_xy_NN, p_ve};
end



%---------OUTPUT MODIFICATION-------------%
% size(p);
press = p * (mu * Ux * L) / (h_ref^2);  %turning dimensional 
% press = p' * (mu * Ux * L) / (h_ref^2);  %turning dimensional 
stress_new = {T_xx, T_yy, T_yy};

end