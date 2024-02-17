function load_dic = calculate_load(solution)

if isempty(solution.pressure_VR)
    W_N= trapz(solution.mesh.x, solution.pressure_LIN("p0"));
    W_NN = trapz(solution.mesh.x, solution.pressure_LIN("p"));
else
    W_N= trapz(solution.mesh.x, solution.pressure_LIN("p0"));
    W_NN = trapz(solution.mesh.x, solution.pressure_LIN("p"));
    W_VR = trapz(solution.mesh.x, solution.pressure_VR("p"));
end
 load_dic = containers.Map({char('W_N'), char('W_NN'), char('W_VR')}, {W_N, W_NN, W_VR});           
end