function [new_data] = data_organizer(data, command, no)
%seperates or combine data
%data type is ehl_solution, 

l_global = 0;

if length(data) == 1
    [m, l, k] = size(data{1});
    ex_data = data{1};
else
    for i = 1:length(data)
        [m, l, k] = size(data{i});
        l_global = l_global + l;
    end
     
end
    
l_global
if command == "merge_u"
    new_data(m, l_global, k) = solution;
    c = 1;
    for i = 1:length(data) 
        ex_data = data{i};
        new_data(:,c:c-1+length(ex_data(1,:,1)),:) = ex_data;  
        c = c + length(ex_data(1,:,1));
        c
    end
end

if command == "divide_de"
   
   new_data(m, l, k-1) = solution;
   index = [];
   for i = 1:k
      if i ~= no
         index = [index i]; 
      end
   end
   
   for j = 1:k-1
      new_data(:,:,j) = ex_data(:,:,index(j));  
   end
end

if command == "divide_u"
   
   new_data(m, l-1, k) = solution;
   index = [];
   for i = 1:l
      if i ~= no
         index = [index i]; 
      end
   end
   
   for j = 1:l-1
      new_data(:,j,:) = ex_data(:,index(j),:);  
   end
end


load_list = [new_data(:,1,1).applied_load];
beta_list = [new_data(1,1,1).viscocity_ratio];
de_list = [new_data(1,1,:).deborah_Number];
velocity_list = [new_data(1,:,1).velocity];
savenote = "";
% 
% ehl_savedata(load_list, beta_list, de_list, velocity_list, new_data, "De", savenote)


end