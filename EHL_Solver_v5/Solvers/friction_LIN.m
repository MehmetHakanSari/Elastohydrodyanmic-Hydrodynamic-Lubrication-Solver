function friction_cell = friction_LIN()


Wax = trapz(x_phys, h .* (dpdx + De * dpDdx));


x_phys = solution.domain.x;

if length(x_phys) ~= N
    x_phys = interp1(linspace(0,1,length(x_phys)), x_phys, linspace(0,1,N));
end

%Dimensional friction
friction_top = - (trapz(x_phys,T_xy_full(1,:)) + De * trapz(x_phys,TD_xy(1,:))) * solution.velocity * solution.domain.mu / solution.domain.href;
friction_bottom = (trapz(x_phys,T_xy_full(end,:)) + De * trapz(x_phys,TD_xy(end,:))) * solution.velocity * solution.domain.mu / solution.domain.href;
% friction_viscoelastic = (- De * trapz(x_phys,TD_xy(1,:)) + De * trapz(x_phys,TD_xy(end,:))) * solution.velocity * solution.domain.mu / solution.domain.href;
friction_viscoelastic_bottom = De * trapz(x_phys,TD_xy(end,:)) * solution.velocity * solution.domain.mu / solution.domain.href;
friction_viscoelastic_top = - De * trapz(x_phys,TD_xy(1,:)) * solution.velocity * solution.domain.mu / solution.domain.href;


if isempty(solution.applied_load)
    generated_load = (x_phys(2) - x_phys(1)) * trapz(p);
    friction_coeff = - (friction_bottom) / generated_load;
else
    friction_coeff = - (friction_bottom) / solution.applied_load;
end



friction_cell = containers.Map({"f_t", "f_b", "f_ve_t","f_ve_b", "f_coeff"}, ...
        {friction_top, friction_bottom, friction_viscoelastic_top, friction_viscoelastic_bottom, friction_coeff});

end