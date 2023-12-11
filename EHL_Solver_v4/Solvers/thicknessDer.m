function [dhdx,a,c,u,cn,cs] = thicknessDer(h, dx)

ddx = 2 * dx;
d2x = dx^2;

dhdx = OneDcentraldiff(h, dx);

% Coeffs. discretized equation (How it come to me)
% for i = 1:nx    
%         a = 1;
%         c(i) = (3 / h(i))*(dhdx(i));
%         u = 2*(a/d2x);
%         cn(i) = ( (a/d2x) + (c(i)/ddx)) / u;
%         cs(i) = ( (a/d2x) - (c(i)/ddx)) / u; 
% end 

% Coeffs. discretized equation (How I fixed it)
% a = 1;
% for i = 1:nx    
%         
%         c(i) = 3 * dhdx(i) / h(i);
%         u = 2 / d2x;
%         cn(i) = ( 1/d2x + c(i)/ddx ) / u;
%         cs(i) = ( 1/d2x - c(i)/ddx ) / u; 
% end 

%%OR 

a = 1;
c = 3 * dhdx ./ h;
u = 2 / d2x;
cn = 1/2 + c / (u * ddx);
cs = 1/2 - c / (u * ddx);

end