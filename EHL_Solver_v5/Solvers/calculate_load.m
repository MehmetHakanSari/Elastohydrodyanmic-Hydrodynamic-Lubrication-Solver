function load_dic = calculate_load(solution, cav_flag)

if isempty(solution.pressure_VR)
    W_N= trapz(solution.mesh.x, solution.pressure_LIN("p0"));
    W_NN = trapz(solution.mesh.x, solution.pressure_LIN("p"));
else
    W_N= trapz(solution.mesh.x, solution.pressure_LIN("p0"));
    W_NN = trapz(solution.mesh.x, solution.pressure_LIN("p"));
    W_VR = trapz(solution.mesh.x, solution.pressure_VR("p"));
    if cav_flag == "off"
    W_VE = trapz(solution.mesh.x, solution.pressure_VR("p_ve"));
    elseif cav_flag == "on" 
    W_VE = 0;
    end
end
 load_dic = containers.Map({char('W_N'), char('W_NN'), char('W_VR'), char('W_VE')}, ...
     {W_N, W_NN, W_VR, W_VE}); 
end