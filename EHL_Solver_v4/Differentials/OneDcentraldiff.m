function Ans = OneDcentraldiff(func, dx)
%Takes 1D array and spacing array (or a single constant spacing)
%Taking three point forward and backward differancing in boundries

[m,n] = size(func);

if m > 1 && n == 1  % if it is a column vector
    func = func';
end

n = length(func);
dfuncdx = zeros(1, n);

dfuncdx(1) = (-func(3) + 4 * func(2) - 3 * func(1)) ./ (2 * dx);
dfuncdx(2:end-1) = (func(3:end) - func(1:end-2)) ./ (2 * dx);
dfuncdx(end) = (func(end-2) - 4 * func(end-1) + 3 * func(end)) ./ (2 * dx);

Ans = dfuncdx;
end


