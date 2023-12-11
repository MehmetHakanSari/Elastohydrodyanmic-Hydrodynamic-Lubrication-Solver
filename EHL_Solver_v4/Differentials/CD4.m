function Ans = CD4(func, dx, dy)

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

for j = 1:n
    for i = 1:m
        switch i 
            case 1 
%                 dfuncdx(j,i) = (-func(j,i+2) + 4 * func(j,i+1) - 3 * func(j,i)) / (2 * dx); %Forward 
%                 func(j,i) = (func(j,i+2) - func(j,i+1)) / (2 * dx);
                dfuncdx(j,i) = 0 / (2 * dx);
            case 2
                dfuncdx(j,i) = (-func(j,i+2) + 4 * func(j,i+1) - 3 * func(j,i)) / (2 * dx); %Forward 
            case (m-1) 
                dfuncdx(j,i) = (func(j,i-2) - 4 * func(j,i-1) + 3 * func(j, i)) / (2 * dx);
            case m
                dfuncdx(j,i) = (func(j,i-2) - 4 * func(j,i-1) + 3 * func(j, i)) / (2 * dx); %Backward
            otherwise
                dfuncdx(j,i) = (-func(j,i+2) + 8*func(j,i+1) - 8*func(j,i-1) + func(j,i-2)) / (12 * dx);     %Central 
        end
    end
end

for i = 1:m
    for j = 1:n         %If (y(0)) != (y = 0). the case when matrix index is not coordinate. 
        switch j 
            case 1 
                dfuncdy(j,i) = (func(j+2,i) - 4 * func(j+1,i) + 3 * func(j,i)) / (2 * dy(i));  %Forward
            case 2
                dfuncdy(j,i) = (func(j+2,i) - 4 * func(j+1,i) + 3 * func(j,i)) / (2 * dy(i));                                 
            case n-1
                dfuncdy(j,i) = (-func(j-2,i) + 4 * func(j-1,i) - 3 * func(j,i)) / (2 * dy(i)); %Backward                               
            case n
%                 dfuncdy(j,i) = (-func(j-2,i) + 4 * func(j-1,i) - 3 * func(j,i)) / (2 * dy(i)); %Backward
                dfuncdy(j,i) = 0;
            otherwise
                dfuncdy(j,i) = (-func(j-2,i)+8*func(j-1,i) - 8*func(j+1,i)+func(j+2,i)) / (12 * dy(i)); %Central 
        end
    end
end
Ans = {dfuncdx dfuncdy};
end
