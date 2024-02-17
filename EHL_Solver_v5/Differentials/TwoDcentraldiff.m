function Ans = TwoDcentraldiff(func, dx, dy)

%func is 2D Array. Each node has its own varabiles. 
%dx is the spacing in x direction.
%dy is 1xN array as griding is not equally spaced each dy values has
    %special space w.r.t x coordinate.
    
%returns funcs gradient in a cell array. First index is d/dx where Second index is d/dy. 

n = length(func(:,1));    %row    #    (mostly represent y in cartesian coordinate)
m = length(func(1,:));    %column #    (mostly represent x in cartesian coordinate)
dfuncdx = zeros(n,m);
dfuncdy = zeros(n,m);

if n < 3
    error("Error: array size is too small!!");
end
if m < 3
    error("Error: array size is too small!!");
end

% for j = 1:n
%     for i = 1:m
%         switch i 
%             case 1 
% %                 dfuncdx(j,i) = (-func(j,i+2) + 4 * func(j,i+1) - 3 * func(j,i)) / (2 * dx); %Forward 
% %                 func(j,i) = (func(j,i+2) - func(j,i+1)) / (2 * dx);
%                 dfuncdx(j,i) = 0 / (2 * dx);
%             case m
%                 dfuncdx(j,i) = (func(j,i-2) - 4 * func(j,i-1) + 3 * func(j, i)) / (2 * dx); %Backward
%             otherwise
%                 dfuncdx(j,i) = (func(j,i+1) - func(j,i-1)) / (2 * dx);                      %Central 
%         end
%     end
% end


dfuncdx(:,1) = (-func(:,3) + 4 * func(:,2) - 3 * func(:,1)) ./ (2 * dx);
dfuncdx(:,2:end-1) = (func(:,3:end) - func(:,1:end-2)) ./ (2 *dx);
dfuncdx(:,end) = (func(:,end-2) - 4 * func(:,end-1) + 3 * func(:, end)) ./ (2 * dx);

% for i = 1:m
%     for j = 1:n         %If (y(0)) != (y = 0). the case when matrix index is not coordinate. 
%         switch j 
%             case 1 
%                 dfuncdy(j,i) = (func(j+2,i) - 4 * func(j+1,i) + 3 * func(j,i)) / (2 * dy(i));  %Forward
%             case n
% %                 dfuncdy(j,i) = (-func(j-2,i) + 4 * func(j-1,i) - 3 * func(j,i)) / (2 * dy(i)); %Backward
%                 dfuncdy(j,i) = 0;
%             otherwise
%                 dfuncdy(j,i) = (func(j-1,i) - func(j+1,i)) / (2 * dy(i));                      %Central 
%         end
%     end
% end

%if dy is row vector the matrix elementwise division is done by each row
%element. 

% dy = dy';

dfuncdy(1,:) = (func(3,:) - 4 * func(2,:) + 3 * func(1,:)) ./ (2 * dy);
dfuncdy(2:end-1,:) = (func(1:end-2,:) - func(3:end,:)) ./ (2 * dy); 
dfuncdy(end,:) = (-func(end-2,:) + 4 * func(end-1,:) - 3 * func(end,:)) ./ (2 * dy);

%To check their performance
% figure(1)
% surf(dfuncdx - dfuncdx1, "Linestyle", "none")
% figure(2)
% surf(dfuncdy - dfuncdy1, "Linestyle", "none")

Ans = {dfuncdx dfuncdy};
end

  
    
    
    
    
