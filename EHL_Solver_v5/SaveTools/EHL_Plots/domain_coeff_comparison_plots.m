function wi_number = domain_coeff_comparison_plots(solution, call_index, savepath, wi_index)
% solution contains (k x 1 x m) sized solution class matrix. 

savepath = string(savepath);

[k, l, m] = size(solution);

line_type = ["-", "--", ":"," -."];
% legend_list = [];

wi_number = solution(1,1,wi_index).wiessenberg_Number;

for i = 1:k
    figure(1); hold on
    plot([solution(i,1,wi_index).domain.x] ./ [solution(i,1,wi_index).domain.bH], [solution(i,1,wi_index).pressure], line_type(call_index),"LineWidth",2)
    figure(2); hold on
    plot([solution(i,1,wi_index).domain.x] ./ [solution(i,1,wi_index).domain.bH], [solution(i,1,wi_index).h], line_type(call_index), "LineWidth",2)
end

figure(1); 
xlabel("\boldmath{$x/b_{H}$}", "Interpreter","latex")
ylabel("\bf{p (Pa)}", "Interpreter","latex")
title("\bf{Wi} = " + string(solution(i,1,wi_index).wiessenberg_Number) + ", \bf{U} = " + string(solution(i,1,wi_index).velocity) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)

% if velocity == 0.1 
%     xlim([-20 20])
% elseif velocity == 1
%     xlim([-20 20])
% end

figure(2); 
xlabel("\boldmath{$x/b_{H}$}","Interpreter","latex")
xlim([-3 3])
ylabel("\bf{h (m)}","Interpreter","latex")
title("\bf{Wi} = " + string(solution(i,1,wi_index).wiessenberg_Number) + ", \bf{U} = " + string(solution(i,1,wi_index).velocity) ,"Interpreter","latex")
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24)


pressure_plot = figure(1);
gap_plot = figure(2);

set(pressure_plot,'Position',[10 10 980 720])
set(gap_plot,'Position',[10 10 980 720])

% savepath = 'E:\Bilkent Dökümanları\Masterımsı\Haftalık\Dec 15 22\';

if ~isfolder(savepath + "/svg")
    mkdir(savepath + "/svg")
end
if ~isfolder(savepath + "/fig")
    mkdir(savepath + "/fig")
end
if ~isfolder(savepath + "/svg/gap")
    mkdir(savepath + "/svg/gap")
end
if ~isfolder(savepath + "/svg/pressure")
    mkdir(savepath + "/svg/pressure")
end
if ~isfolder(savepath + "/fig/gap")
    mkdir(savepath + "/fig/gap")
end
if ~isfolder(savepath + "/fig/pressure")
    mkdir(savepath + "/fig/pressure")
end



end