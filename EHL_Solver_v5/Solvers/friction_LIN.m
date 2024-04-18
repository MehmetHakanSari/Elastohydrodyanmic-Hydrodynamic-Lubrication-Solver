function friction_cell = friction_LIN(load_cell, mesh, stress_VR, stress_LIN, De, p)


% Wax = trapz(x_phys, h .* (dpdx + De * dpDdx));

%default = VR
T_xy_full = stress_VR("tau_xy");
TD_xy = stress_LIN("tau_xy");


%Dimensional friction
friction_top = - trapz(mesh.x, T_xy_full(1,:));
friction_bottom = trapz(mesh.x,T_xy_full(end,:)); 
% friction_viscoelastic = (- De * trapz(x_phys,TD_xy(1,:)) + De * trapz(x_phys,TD_xy(end,:)));
friction_viscoelastic_bottom = De * trapz(mesh.x, TD_xy(end,:)) ;
friction_viscoelastic_top = - De * trapz(mesh.x,TD_xy(1,:));


if isempty(load_cell)
    generated_load = mesh.dx * trapz(p);
    friction_coeff = - (friction_bottom) / generated_load;
else
    friction_coeff = - (friction_bottom) / load_cell("W_VR");
end

friction_cell = containers.Map({char("f_t"), char("f_b"), char("f_ve_t"), char("f_ve_b"), char("f_coeff")}, ...
        {friction_top, friction_bottom, friction_viscoelastic_top, friction_viscoelastic_bottom, friction_coeff});

end