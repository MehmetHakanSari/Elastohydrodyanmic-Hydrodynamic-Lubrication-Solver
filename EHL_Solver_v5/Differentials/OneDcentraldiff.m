function Ans = OneDcentraldiff(func, dx, varargin)
%Takes 1D array and spacing array (or a single constant spacing)
%Taking three point forward and backward differancing in boundries

[m,n] = size(func);

if m > 1 && n == 1  % if it is a column vector
    func = func';
end

n = length(func);
dfuncdx = zeros(1, n);
Scheme = "CD2";
 if nargin > 1
            for i = 1:numel(varargin)
                Scheme = varargin{i};
            end
 end
 
 switch Scheme
     case "CD2" 
    %CD2 
    dfuncdx(1) = (-func(3) + 4 * func(2) - 3 * func(1)) ./ (2 * dx);
    dfuncdx(2:end-1) = (func(3:end) - func(1:end-2)) ./ (2 * dx);
    dfuncdx(end) = (func(end-2) - 4 * func(end-1) + 3 * func(end)) ./ (2 * dx);
     case "CD4"
    %CD4
    dfuncdx(1:2) =  (-25/12 * func(1:2) + 4 * func(2:3)	- 3	* func(3:4) + 4/3 * func(4:5) - 1/4 * func(5:6) )	./ (dx);
    dfuncdx(3:end-2) = (-1/12 * func(5:end) + 2/3 * func(4:end-1) - 2/3 * func(2:end-3) + 1/12 * func(1:end-4)) ./ (dx);
    dfuncdx(end-1:end) =  (25/12 * func(end-1:end) - 4 * func(end-2:end-1) + 3	* func(end-3:end-2)...
                                                                                                                            - 4/3 * func(end-4:end-3) + 1/4 * func(end-5:end-4)) ./ (dx);
                                                                                                                        
     case "CD6"
     dfuncdx(1:3) =  (-49/20 * func(1:3) + 6 * func(2:4) -15/2 * func(3:5)	+  20/3	* func(4:6) ...
                                                                                -  15/4 * func(5:7) + 6/5 * func(6:8) - 1/6 * func(7:9) ) ./ (dx);
     dfuncdx(4:end-3) = (-func(1:end-6) + 9 * func(2:end-5) - 45 * func(3:end-4)  ...
                                                                                                                + 45 * func(5:end-2) -  9 * func(6:end-1) +  func(7:end)) ./ (60 * dx);
     dfuncdx(end-2:end) =   (49/20 * func(end-2:end) - 6 * func(end-3:end-1) + 15/2 * func(end-4:end-2)	- 20/3 * func(end-5:end-3) ...
                                           + 15/4 * func(end-6:end-4) - 6/5 * func(end-7:end-5) + 1/6 * func(end-8:end-6)  ) / ( dx); 
     case "UW"
%     dfuncdx(1) = (-func(3) + 4 * func(2) - 3 * func(1)) ./ (2 * dx);
%     dfuncdx(2) = (-func(4) + 4 * func(3) - 3 * func(2)) ./ (2 * dx);
%     dfuncdx(3:end) = (func(1:end-2) - 4 * func(2:end-1) + 3 * func(3:end)) ./ (2 * dx);
    dfuncdx(1) = (func(2) - func(1)) ./ dx;
    dfuncdx(2:end) = (func(2:end) - func(1:end-1)) ./ dx;
     case "DW"
%     dfuncdx(1:end-2) = (func(3:end) - 4 * func(2:end-1) + 3 * func(1:end-2)) ./ (2 * dx);
%     dfuncdx(end-1) = (func(end-3) - 4 * func(end-2) + 3 * func(end-1)) ./ (2 * dx);
%     dfuncdx(end) = (func(end-2) - 4 * func(end-1) + 3 * func(end)) ./ (2 * dx);
    dfuncdx(1:end-1) = (func(2:end) - func(1:end-1)) ./ dx;
    dfuncdx(end) = (func(end) - func(end-1)) ./ dx;
     otherwise
    %CD2 
    dfuncdx(1) = (-func(3) + 4 * func(2) - 3 * func(1)) ./ (2 * dx);
    dfuncdx(2:end-1) = (func(3:end) - func(1:end-2)) ./ (2 * dx);
    dfuncdx(end) = (func(end-2) - 4 * func(end-1) + 3 * func(end)) ./ (2 * dx);
 end


Ans = dfuncdx;
end


