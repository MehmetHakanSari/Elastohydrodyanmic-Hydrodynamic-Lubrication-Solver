function [obj, dspflag] = displacementsolver(obj)
            
    d_old = obj.displacement;
    pressure_D = obj.pressure * (obj.fluid_domain.eta_o * obj.U_o * obj.fliud_domain.length) / (obj.fluid_domain.h_ref^2)
                        
    % Elastic displacement
    d = (obj.solid_domain.deformation_kernel * (pressure_D)')';    % (Nx1) = (NxN) (Nx1) changing Dps from (Nx1) to (1xN)
    
    % Aitken acceleration
    r_old = (obj.r)';                           %transpoze due to store preferences. 
    
    r = ( abs( (d - d_old) ./ d) )';
    
    relaxD = obj.simulation_proporties.relaxD * (r_old' * (r_old - r)) ./ ((r - r_old)' * (r - r_old));       
    relaxD = min(max(simulation_proporties.relaxD, obj.simulation_proporties.relaxDmin), ...
                                                                obj.simulation_proporties.relaxDmax);
     
    % error on displacament
    errD = max(abs((d - d_old) ./ d));
    % The displacement is underelaxated
    d = (1 - obj.simulation_proporties.relaxD) * d_old + obj.simulation_proporties.relaxD * d;
    h_new = obj.h0 + obj.domain.height_profile - d;      %the h0 changes when each load balance occurs. 
    
    obj.r = r';
    obj.displacement = d;
    
    if errD < obj.simulation_proporties.tollD
        dspflag = false;
    else
        obj.h = h_new;
        dspflag = true;
    end
    
end