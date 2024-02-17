function [p] = ReynoldsCondiation(h, p_previous, mesh, fluid, simulation_proporties)   
    dhdx = OneDcentraldiff(h, mesh.dx)
    c = 3 * dhdx ./ h;
    cn = 1/2 + c / (4 * dx);
    cs = 1/2 - c / (4 * dx);

    e = -(6 * dhdx) ./ h.^3;  
    g = e / (2 / dx^2);                     
    
    p = p_previous; 
    p(1) = fluid.BC("p_i"); p(end) = fluid.BC("p_o");
    relax = simulation_proporties.relaxp;
    
    convergence = false;
    counter = 0;
    while ~convergence
        counter = counter + 1; 
        p_old = p;
        
        for i = 2:mesh.N_x-1 
            p(i) = (1 - relax) * p_old(i) + relax * (cn(i)*p(i+1) + cs(i)*p(i-1) + g(i));
            if p(i) < fluid.BC("p_c")
                p = fluid.BC("p_c");
            end
        end
        % Error calculation
        err = max(abs((p - p_old) ./ p));
        if err < simulation_proporties.tollp
            convergence = true
        end
    end
    
%flowrate = h'/2 * Ux - h'.^3 / 12 / mu .* OneDcentraldiff(press, dx);

end



