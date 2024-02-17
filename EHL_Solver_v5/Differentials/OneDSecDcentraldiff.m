function Ans = OneDSecDcentraldiff(func, dx, varargin)
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

     case "CD4"
    %CD4
    dfuncdx(1:2) = (2 * func(1:2) - 5 * func(2:3) + 4 * func(3:4) - func(4:5) ) ./ (dx.^2);
    dfuncdx(3:end-2) = (-func(1:end-4) + 16 * func(2:end-3) - 30 * func(3:end-2) + 16 * func(4:end-1) - func(5:end)) ./ (12 * dx.^2);
    dfuncdx(end-1:end) = (2 * func(end-4:end-3) - 5 * func(end-3:end-2) + 4 * func(end-2:end-1) - func(end-1:end) ) ./ (dx.^2);
     case "CD6"
     dfuncdx(1:3) =  (469/90 * func(1:3) - 223/10 * func(2:4) + 	879/20 * func(3:5)	- 949/18	* func(4:6) ...
                                                                                + 41 * func(5:7) - 201/10 * func(6:8) + 1019/180 * func(7:9) - 7/10 *  func(8:10) ) / ( dx.^2);
     dfuncdx(4:end-3) = (2 * func(1:end-6) - 27 * func(2:end-5) + 270 * func(3:end-4) - 490 * func(4:end-3) ...
                                                                                                                + 270 * func(5:end-2) -  27 * func(6:end-1) + 2 * func(7:end)) ./ (180 * dx.^2);
     dfuncdx(end-2:end) =   (469/90 * func(end-2:end) - 223/10 * func(end-3:end-1) + 	879/20 * func(end-4:end-2)	- 949/18	* func(end-5:end-3) ...
                                           + 41 * func(end-6:end-4) - 201/10 * func(end-7:end-5) + 1019/180 * func(end-8:end-6) - 7/10 *  func(end-9:end-7) ) / ( dx.^2);                                                    
     case "UW"

     case "DW"

     otherwise
    %CD2 
    dfuncdx(1) = (func(3) - 2 * func(2) + func(1)) ./ (dx^2);
    dfuncdx(2:end-1) = ( func(3:end) - 2 * func(2:end-1)  + func(1:end-2) ) ./ (dx^2);
    dfuncdx(end) = ( func(end-2) - 2 * func(end-1) + func(end) ) ./ (dx^2);
 end


Ans = dfuncdx;
end


