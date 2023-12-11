function pressure_comparision(s1, s2)
%Compares s1 and s2 pressures in graphics
%s1 is the solution class variable
%s2 is the solution class variable



figure(1); hold on; box on;
plot(s1(2,1,2).domain.x / s1(2,1,2).domain.bH, s1(2,1,2).pressure,"k","LineWidth", 2)
plot(s1(2,1,5).domain.x / s1(2,1,5).domain.bH, s1(2,1,5).pressure,"k--","LineWidth", 2)
plot(s2(4,1,2).domain.x / s2(4,1,2).domain.bH, s2(4,1,2).pressure,"b","LineWidth", 2)
plot(s2(4,1,5).domain.x / s2(4,1,5).domain.bH, s2(4,1,5).pressure,"b--","LineWidth", 2)
xlim([-5 2])


s1_1 = find(s1(2,1,1).pressure(2:end) == 0, 1, 'first');
s1_2 = find(s1(2,1,5).pressure(2:end) == 0, 1, 'first');

s2_1 = find(s1(1,1,1).pressure(2:end) == 0, 1, 'first');
s2_2 = find(s1(1,1,5).pressure(2:end) == 0, 1, 'first');



disp("The x_cav shifted: " + string(100 * abs(s1(2,1,1).domain.x(s1_1) - s1(2,1,5).domain.x(s1_2)) / (s1(2,1,1).domain.x(s1_1) - s1(2,1,1).domain.x(1))) + " percentage")

disp("The x_cav shifted: " + string(100 * abs(s2(2,1,1).domain.x(s2_1) - s2(2,1,5).domain.x(s2_2)) / (s2(2,1,1).domain.x(s2_1) - s2(2,1,1).domain.x(1))) + " percentage")

figure(2); hold on; box on;
plot(s2(1,1,2).domain.x / s2(1,1,2).domain.bH, s2(1,1,2).displacement / min(s1(1,1,2).displacement),"k","LineWidth", 1.7)
plot(s2(1,1,5).domain.x / s2(1,1,5).domain.bH, s2(1,1,5).displacement / min(s1(1,1,5).displacement),"k--","LineWidth", 1.7)
plot(s1(2,1,2).domain.x / s1(2,1,2).domain.bH, s1(2,1,2).displacement / min(s1(2,1,2).displacement),"m","LineWidth", 1.7)
plot(s1(2,1,5).domain.x / s1(2,1,5).domain.bH, s1(2,1,5).displacement / min(s1(2,1,5).displacement),"m--","LineWidth", 1.7)
plot(s2(3,1,2).domain.x / s2(3,1,2).domain.bH, s2(3,1,2).displacement / min(s1(3,1,2).displacement),"r","LineWidth", 1.7)
plot(s2(3,1,5).domain.x / s2(3,1,5).domain.bH, s2(3,1,5).displacement / min(s1(3,1,5).displacement),"r--","LineWidth", 1.7)
plot(s2(4,1,2).domain.x / s2(4,1,2).domain.bH, s2(4,1,2).displacement / min(s1(4,1,2).displacement),"b","LineWidth", 1.7)
plot(s2(4,1,5).domain.x / s2(4,1,5).domain.bH, s2(4,1,5).displacement / min(s1(4,1,5).displacement),"b--","LineWidth", 1.7)
xlim([-5 2])


figure(3); hold on; box on;
plot(s1(2,1,2).domain.x / s1(2,1,2).domain.bH, s1(2,1,2).h / s1(2,1,2).domain.href ,"k","LineWidth", 2)
plot(s1(2,1,5).domain.x / s1(2,1,5).domain.bH, s1(2,1,5).h / s1(2,1,5).domain.href,"k--","LineWidth", 2)
plot(s2(4,1,2).domain.x / s2(4,1,2).domain.bH, s2(4,1,2).h / s1(4,1,2).domain.href,"b","LineWidth", 2)
plot(s2(4,1,5).domain.x / s2(4,1,5).domain.bH, s2(4,1,5).h / s1(4,1,5).domain.href,"b--","LineWidth", 2)
xlim([-2 2])

end