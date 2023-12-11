# EHL - HL Solver
A channel solver for low Re numbers fluid flow. Non-Newtonian effects is primaly caused by viscolelasticity caused by polymeric additives in fluid.The solver is capabale of solving EHL and HL cases by introduced variables. First option is to run EHL with defined Load, Speed, Viscosity Ratio and Deborah Number. For a fully Newtonian case, Deborah number should be zero. 
---------------------------------------
For HL conditions user should define a contact type such as parabolic slider, cosine slider, step etc. from domain script. After defining geometry, simple Reynolds solution can be obtained by Linearized Reynolds (LIN) solver via setting De = 0 for Newtonian solution. However, for non-Newtonian solvers the user should define De. There is two solves considering viscoelasticity, LIN and Viscoelastic Reynolds (VR). Meshing is required for solving VR. User can switch on and off cavitation by giving cavitation flag "on" or "off". 
---------------------------------------
This code is written for academic purposes and the output of the research is about to publish soon. Regarding this fact, the solver is not fulyy functional in terms of data processing. However, one solution class can be stored. The stored solution class contains:
---------------------------------------
solution:
pressure,
height,
friction,
velocity field,
stress field,
pressure field,
pressure_VR
stress_VR
pressure0_old
pressure1_old
Deborah number
---------------------------------------



