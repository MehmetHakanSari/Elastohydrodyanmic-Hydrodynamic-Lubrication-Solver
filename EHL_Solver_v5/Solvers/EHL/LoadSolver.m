function [obj, loadflag] = LoadSolver(obj, F) 
    %F: given force
    %pressure: pressure output taken from the finished displacement
    %solver.
    pressure_D = obj.pressure * (obj.fluid_domain.eta_o * obj.U_o * obj.fliud_domain.length) / (obj.fluid_domain.h_ref^2)

    F_new = real(trapz(obj.mesh.dimensional("x"), pressure_D));
    
    % error on Load
    errL = abs(F_new - F) / abs(F);            
    
    if errL < obj.simulation_proporties.tollL
        loadflag = false;
    else % Load Control on gentral gap
        
        h00 = obj.h0 - obj.displacement(floor(domain.N_x / 2 + 1));  %displacement is column vector
        h00n = h00 * (F_new / F) ^ obj.simulation_proporties.relaxL;
        h0n = h00n + obj.displacement(floor(obj.mesh.N_x / 2 + 1));
        
        obj.h0 = h0n;
        obj.h = obj.h0 + domain.height_profile - obj.displacement;  
        loadflag = true;
    end
          
end