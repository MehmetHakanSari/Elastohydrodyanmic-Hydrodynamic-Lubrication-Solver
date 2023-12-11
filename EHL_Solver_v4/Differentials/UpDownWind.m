function Ans = UpDownWind(func, dx, dy, u, v)

%func is 2D Array. Each node has its own varabiles. 
%u is x-velocity, matrixs length should equal to func         
%v is v-velocity, matrixs length should equal to func
%dx is the spacing in x direction.
%dy is 1xN array as griding is not equally spaced each dy values has
    %special space w.r.t x coordinate.
    
%returns funcs gradient in a cell array. First index is d/dx where Second index is d/dy. 

[m, n] = size(func);

if n < 3
    error("Error: array size is too small!!");
end
if m < 3
    error("Error: array size is too small!!");
end

dfuncdx = zeros(m, n);
dfuncdy = zeros(m, n);

dfuncdx_uw = zeros(m, n); 
dfuncdy_uw = zeros(m, n);
dfuncdx_dw = zeros(m, n); 
dfuncdy_dw = zeros(m, n);

u_uw_logic = u > 0;
v_uw_logic = v > 0;
u_dw_logic = u <= 0;
v_dw_logic = v <= 0;
% 

dfuncdx_dw(:,2:end-1) = u_dw_logic(:, 2:end-1) .* (func(:, 3:end) - func(:,2:end-1)) ./ dx;   %forward
dfuncdy_dw(2:end-1,:) = v_dw_logic(2:end-1, :) .* (-func(3:end, :) + func(2:end-1,:)) ./ dy;  %check sign. should be reverse. 

dfuncdx_uw(:,2:end-1) = u_uw_logic(:, 2:end-1) .* (func(:, 2:end-1) - func(:,1:end-2)) ./ dx;   %backward   
dfuncdy_uw(2:end-1,:) = v_uw_logic(2:end-1, :) .* (-func(2:end-1, :) + func(1:end-2,:)) ./ dy;  %check sign. should be reverse. 


dfuncdx = dfuncdx_uw + dfuncdx_dw;
dfuncdy = dfuncdy_uw + dfuncdy_dw;

% %boundaries are FW and BW respectievly regardless of velocity. 
% 
dfuncdx(:,1) = (func(:,2) - func(:,1)) ./ dx;
dfuncdx(:,end) = (func(:,end) - func(:,end-1)) ./ dx;

dfuncdy(1,:) = (-func(2,:) + func(1,:)) ./ dy;
dfuncdy(end,:) = (func(end-1,:) - func(end,:)) ./ (dy);


% for j = 1:m
%     for i = 1:n                  %upwind
%         switch i
%             case 1
%                 dfuncdx(j, i) = (func(j,i+1) - func(j,i)) ./ (dx);     %Forward differance
%             case n
%                 dfuncdx(j,i) = ((func(j, i) - func(j,i-1))) ./ (dx);  %Backward differance
%             otherwise
%                 dfuncdx(j,i) = ((func(j, i) - func(j,i-1))) ./ (dx);  %Backward differance
%         end
%     end
% end
% 
% for i = 1:n
%     for j = 1:m                            %If (y(0)) != (y = 0). the case when matrix index is not coordinate. 
%         switch j
%             case 1
%                 dfuncdy(j,i) = (func(j+2,i) - 4 * func(j+1,i) + 3 * func(j,i)) ./ (2 * dy(i));  %Forward Differance
%             case m
%                 dfuncdy(j,i) = (-func(j-2,i) + 4 * func(j-1,i) - 3 * func(j,i)) ./ (2 * dy(i)); %Backward Differance
%             otherwise
%                 dfuncdy(j,i) = (func(j-1,i) - func(j+1,i)) ./ (2 * dy(i));                      %Central Differance
%         end
%     end
% end

func
dfuncdx
dfuncdy

Ans = {dfuncdx dfuncdy};
end




% for j = 1:n
%     for i = 1:m
%         if u(j,i) < 0      %downwind 
%             switch i
%                 case 1
%                     dfuncdx(j,i) = (func(j,i+1) - func(j,i)) / (dx);          %Forward differance
%                 case m
%                     dfuncdx(j,i) = (func(j, i) - func(j,i-1)) / (dx);         %Backward differance
%                 otherwise
%                     dfuncdx(j,i) = (func(j,i+1) - func(j,i)) / (dx);          %Forward differance
%             end
%         else                    %upwind
%             switch i
%                 case 1
%                     dfuncdx(j,i) = (func(j,i+1) - func(j,i)) / (dx);     %Forward differance
%                 case m
%                     dfuncdx(j,i) = ((func(j, i) - func(j,i-1))) / (dx);  %Backward differance
%                 otherwise
%                     dfuncdx(j,i) = ((func(j, i) - func(j,i-1))) / (dx);  %Backward differance
%             end
%         end
%     end
% end
% 
% for i = 1:m
%     for j = 1:n                            %If (y(0)) != (y = 0). the case when matrix index is not coordinate. 
%         if v(i,j) == 0
%             switch j
%                 case 1
%                     dfuncdy(j,i) = (func(j+2,i) - 4 * func(j+1,i) + 3 * func(j,i)) / (2 * dy(i));  %Forward Differance
%                 case n
%                     dfuncdy(j,i) = (-func(j-2,i) + 4 * func(j-1,i) - 3 * func(j,i)) / (2 * dy(i)); %Backward Differance
%                 otherwise
%                     dfuncdy(j,i) = (func(j-1,i) - func(j+1,i)) / (2 * dy(i));                      %Central Differance
%             end 
%         elseif v(i,j) < 0      %downwind
%             switch j
%                 case 1
%                     dfuncdy(j,i) = (-func(j+1,i) + func(j,i)) / (dy(i));    %Forward differance
%                 case n
%                     dfuncdy(j,i) = (func(j-1,i) - func(j,i)) / (dy(i));     %Backward differance
%                 otherwise
%                     dfuncdy(j,i) = (-func(j+1,i) + func(j,i)) / (dy(i));    %Forward differance
%             end
%         else                   %upwind
%             switch j
%                 case 1
%                     dfuncdy(j,i) = (-func(j+1,i) + func(j,i)) / (dy(i));    %Forward differance
%                 case n
%                     dfuncdy(j,i) = (func(j-1,i) - func(j,i)) / (dy(i));     %Backward differance
%                 otherwise
%                     dfuncdy(j,i) = (func(j-1,i) - func(j,i)) / (dy(i));     %Backward differance
%             end
%         end
%     end
% end
% 
% Ans = {dfuncdx dfuncdy};
% end

  
    
    
    
    
