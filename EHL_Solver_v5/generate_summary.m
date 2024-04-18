function generate_summary(S)
%S: solution
close all
clc

if  ~isempty(S.PD)
    fprintf('----------------------------------------------\n');
    fprintf('Physical Quantaties of the Channel\n');
    fprintf('----------------------------------------------\n');

    fprintf('U         :      ' + string(S.PD.U)   +   '     m/s  \n');
    fprintf('L          :      ' + string(S.PD.L)   +   '     m  \n');
    fprintf('h0        :      ' + string(S.PD.h0)   +   '     m  \n');
    fprintf('eta0     :      ' + string(S.PD.eta0)   +   '   Pa . s  \n');
%     fprintf(' etas:      ' + string(S.PD.U)   +   '     m/s  \n');
    fprintf('lambda:      ' + string(S.PD.lambda)   +   '    s  \n');
    fprintf(' p0        :      ' + string(S.PD.p0)   +   '    Pa     \n');
else
    fprintf('----------------------------------------------\n');
    fprintf('Physical Quantaties Not Defined\n');
    fprintf('----------------------------------------------\n');
end
fprintf('\n');
fprintf('\n');



fprintf('      Fluid Properties              \n');
fprintf('----------------------------------------------\n');
if ~isempty(S.fluid.name)
    fprintf('Fluid:                       ' + string(S.fluid.name) + '\n');
end  
if ~isempty(S.fluid.beta)
      fprintf('Concentration:     ' + string(S.fluid.beta) + '\n');
end  
if ~isempty(S.fluid.De) 
    fprintf('Eq Debrorah     :    ' + string(S.fluid.De) + '\n');
    S.fluid = S.fluid.PipkinLink(S.PD.h0 / S.PD.L);
    fprintf('Eq Wiessenberg:    ' + string(S.fluid.Wi) + '\n');
end 
press_dic = S.dimensional_dic("pressure_VR");
if min(press_dic("p")) == 0
    fprintf('Cavitation:  ON \n');
else
    fprintf('Cavitation:  OFF \n');
end
fprintf('\n');
fprintf('\n');


fprintf('    Contact Information        \n');
fprintf('----------------------------------------------\n');
if ~isempty(S.dimensional_dic("load_dic"))
       press_dic = S.dimensional_dic("pressure_VR");
       fprintf('Peak Pressure (VR)   :      ' + string(max(press_dic("p")) / 10^6)  + '   MPa   |    ');
       fprintf('Bottom Pressure (VR)  :      ' + string(min(press_dic("p"))  / 10^6) + '   MPa  | \n  ');
end
if ~isempty(S.dimensional_dic("load_dic"))
      load_dic = S.dimensional_dic("load_dic");
      fprintf('Newtonian Load        :      ' + string(load_dic("W_N")) + '   N / m   |    ');
      fprintf('Non-newtonian Load:      ' + string(load_dic("W_VR")) + '  N  / m  \n');
      if load_dic("W_VE") == 0
      fprintf('Additional Load        :      ' + string(abs(load_dic("W_N") - load_dic("W_VR"))) + '  N / m       | ' );
      fprintf('  percentage change:  %% ' + string(abs(load_dic("W_N") - load_dic("W_VR")) / load_dic("W_N") * 100) + '  \n ');    
      else
      fprintf('Additional Load        :      ' + string(load_dic("W_VE")) + '  N  / m      | ' );
      fprintf('  percentage change:  %% ' + string(load_dic("W_VE") / load_dic("W_N") * 100) + '  \n ');
      end
end  
if ~isempty(S.dimensional_dic("friction_dic"))
    fric_dic = S.dimensional_dic("friction_dic");
    fprintf('Friction Force              :    ' + string(fric_dic('f_b')) + '  N / m   \n');
    fprintf('Cofficient of Friction :    ' + string(fric_dic('f_coeff')) + '   \n');
end  
fprintf('Estimate of Maximum Shear Rate U/h_0:   ' + string(S.PD.U / S.PD.h0) + "  s-1  \n");
fprintf('\n');
fprintf('\n');

%Plot pressure and channel

figure(1); 
subplot(1,2, 1); box on;
plot(S.mesh.x, S.h, "k-", "LineWidth", 2.5)
xlabel("x", "Interpreter", "latex")
title("Channel Profile", "Interpreter", "latex")
set(gca,'FontSize',23,'TickLabelInterpreter','latex');
subplot(1,2, 2); box on; hold on;
plot(S.mesh.x, S.pressure_LIN("p0"), "k-", "LineWidth", 2.3)
plot(S.mesh.x, S.pressure_VR("p"), "-.", "LineWidth", 2.3)
xlabel("x", "Interpreter", "latex")
title("Pressure Distibution", "Interpreter", "latex")
legend("Newtonian", "Non-Newtonian", "box", "off", "Interpreter", "latex")
set(gca,'FontSize',23,'TickLabelInterpreter','latex');
set(figure(1),'Position',[10 50 1300 700])
saveas(figure(1), "Summary_Fig" + ".svg" )
close 
    





end