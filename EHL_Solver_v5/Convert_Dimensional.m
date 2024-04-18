function  converted_res = Convert_Dimensional(pressure_VR, pressure_LIN, stress_VR, load_dic, friction_dic, DQ)
%Convert properties to dimensional cases
% DQ= dimensional quantaties


%pressure
% p^* = p * (eta0 * U * L) / h0^2 + p0 

keysArray = keys(pressure_VR);
valuesArray = values(pressure_VR);

pressure_VR_dim = containers.Map;

for i = 1:numel(keysArray)
    key = keysArray{i};
    p_ = valuesArray{i};
    
    p_dim = p_ * (DQ.eta0 * DQ.U * DQ.L) / DQ.h0^2 + DQ.p0;
    
    pressure_VR_dim(key) = p_dim;
end

keysArray = keys(pressure_LIN);
valuesArray = values(pressure_LIN);
pressure_LIN_dim = containers.Map;

for i = 1:numel(keysArray)
    key = keysArray{i};
    p_ = valuesArray{i};
    p_dim = p_ * (DQ.eta0 * DQ.U * DQ.L) / DQ.h0^2 + DQ.p0;
    pressure_LIN_dim(key) = p_dim;
end

%load
% W^* = int_L p dx = int_L^* p^* dx^* = int_L^*   (p * (eta0 * U * L) / h0^2 + p0) dx^*
% = int_L   (p * (eta0 * U * L) / h0^2 + p0) dx * L = int_L   p  dx  * (eta0 * U * L) / h0^2 + p0) * L 

keysArray = keys(load_dic);
valuesArray = values(load_dic);

load_dic_dim = containers.Map;

for i = 1:numel(keysArray)
    key = keysArray{i};
    W_ = valuesArray{i};
    
    W_dim = W_ * ((DQ.eta0 * DQ.U * DQ.L) / DQ.h0^2 + DQ.p0) * DQ.L;
    
    load_dic_dim(key) = W_dim;
end

%stress 
% T_xx ^* = T_xx *  (eta0 * U * L) / h0^2
% T_xy ^* = T_xy *  (eta0 * U) / h0
% T_yy ^* = T_yy *  (eta0 * U) / L

keysArray = keys(stress_VR);
valuesArray = values(stress_VR);

stress_VR_dim = containers.Map;

for i = 1:numel(keysArray)
    key = keysArray{i};
    tau_ = valuesArray{i};
    
    if key == "tau_xx"
        
       tau_dim = tau_ * ((DQ.eta0 * DQ.U * DQ.L) / DQ.h0^2);
        
    elseif key == "tau_xy"
       tau_dim = tau_ * ((DQ.eta0 * DQ.U) / DQ.h0);   
       
    elseif key == "tau_yy"
        tau_dim = tau_ * ((DQ.eta0 * DQ.U) / DQ.L);
    end
    
    stress_VR_dim(key) = tau_dim;
end

%fric
% f^* = f *  (eta0 * U) / h0 * L

keysArray = keys(friction_dic);
valuesArray = values(friction_dic);

friction_dic_dim = containers.Map;

for i = 1:numel(keysArray)
    key = keysArray{i};
    f_ = valuesArray{i};
    
    if key ~= "f_coeff"
    f_dim = f_ * ((DQ.eta0 * DQ.U) / DQ.h0)  * DQ.L;
    end
    
    friction_dic_dim(key) = f_dim;
end

converted_res = containers.Map;

converted_res(char("pressure_VR")) = pressure_VR_dim;
converted_res(char("pressure_LIN")) = pressure_LIN_dim;
converted_res(char("friction_dic")) = friction_dic_dim;
converted_res(char("load_dic")) = load_dic_dim;
converted_res(char("stress_VR")) = stress_VR_dim;

end






