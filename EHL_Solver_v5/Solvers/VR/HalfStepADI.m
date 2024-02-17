function T = HalfStepADI(T, Inclusive, Constants, u, v, dx, dy, dydx, dt, De, scheme)

%First handle y derivative. 
%           t + dt                          t 
% T / De + T + dt v dTdy - (2 T dudx) = dt C + T  

[Ny, Nx] = size(T);

switch scheme
    case "UW"
        
u_uw_logic = u > 0;
v_uw_logic = v > 0;
u_dw_logic = u <= 0;
v_dw_logic = v <= 0;

%%%%%%
gradT = TwoDcentraldiff(T, dx, dy);
DTDx = gradT{1} + gradT{2} .* dydx;

%First Half step: Regard (v dTdy)^{t + dt/2}; (u dTdx)^{t}
dt = dt / 2;

Q = dt * (Constants - u .* DTDx) + T ; 
C = ( dt * (1 / De +  Inclusive) + 1 ); [m, n] = size(C); isSizeOneByOne = (m == 1) && (n == 1);
if isSizeOneByOne
    C = C * ones(Ny, Nx);
end
S = dt * v ./ (2 * dy) ;
N = -dt * v ./ (2 * dy) ; %bacause the matrix is reverse order in y, north has the negatve sign.

for i = 1:Nx
%      T(2:end-1,i) = TDMA(S(1:end-2,i), C(2:end-1, i), N(3:end, i), Q(2:end-1, i));
     T(:, i) = TDMA(S(:, i), C(:, i), N(:, i), Q(:, i));
end


% gradT = UpDownWind(T, dx, dy, u, v);
% DTDx = gradT{1} + gradT{2} .* dydx;
% 
% %First Half step: Regard (v dTdy)^{t + dt/2}; (u dTdx)^{t}
% dt = dt / 2;
% 
% Q = dt * (Constants - u .* DTDx) + T ; 
% C = ( dt * (1 / De +  Inclusive + v / dy .* v_uw_logic - v / dy .* v_dw_logic) + 1 ); [m, n] = size(C); isSizeOneByOne = (m == 1) && (n == 1);
% if isSizeOneByOne
%     C = C * ones(Ny, Nx);
% end
% S = -dt * v ./ dy .* v_uw_logic;
% N = dt * v ./  dy .* v_dw_logic; %bacause the matrix is reverse order in y, north has the negatve sign.
% 
% for i = 1:Nx
% %      T(2:end-1,i) = TDMA(S(1:end-2,i), C(2:end-1, i), N(3:end, i), Q(2:end-1, i));
%      T(:, i) = TDMA(N(:, i), C(:, i), S(:, i), Q(:, i));
% end


%Second Half step: Regard (u dTdx)^{t + dt}; (v dTdv)^{t+dt/2}

gradT = UpDownWind(T, dx, dy, u, v);
dTdy = gradT{2};

Q = dt * (Constants - v .* dTdy) + T;

C = (dt * (1 / De + Inclusive + u / dx .* u_uw_logic - u / dx .* u_dw_logic) + 1); [m, n] = size(C); isSizeOneByOne = (m == 1) && (n == 1);
if isSizeOneByOne
    C = C * ones(Ny, Nx);
end
E = dt * (u .* (1 / dx  + dydx .* dTdy)) .* u_dw_logic;
W = -dt * (u .* (1 /  dx + dydx .* dTdy)) .* u_uw_logic;

for i = 1:Ny
%     T(i, 2:end-1) = TDMA(W(i, 1:end-2), C(i, 2:end-1), E(i, 3:end), Q(i, 2:end-1));
    T(i, :) = TDMA(W(i, :), C(i, :), E(i, :), Q(i, :));
end
        
        
        
    otherwise 
        
gradT = TwoDcentraldiff(T, dx, dy);
DTDx = gradT{1} + gradT{2} .* dydx;

%First Half step: Regard (v dTdy)^{t + dt/2}; (u dTdx)^{t}
dt = dt / 2;

Q = dt * (Constants - u .* DTDx) + T ; 
C = ( dt * (1 / De +  Inclusive) + 1 ); [m, n] = size(C); isSizeOneByOne = (m == 1) && (n == 1);
if isSizeOneByOne
    C = C * ones(Ny, Nx);
end
S = dt * v ./ (2 * dy) ;
N = -dt * v ./ (2 * dy) ; %bacause the matrix is reverse order in y, north has the negatve sign.

for i = 1:Nx
%      T(2:end-1,i) = TDMA(S(1:end-2,i), C(2:end-1, i), N(3:end, i), Q(2:end-1, i));
     T(:, i) = TDMA(S(:, i), C(:, i), N(:, i), Q(:, i));
end


%Second Half step: Regard (u dTdx)^{t + dt}; (v dTdv)^{t+dt/2}

gradT = TwoDcentraldiff(T, dx, dy);
dTdy = gradT{2};

Q = dt * (Constants - v .* dTdy) + T;
C = (dt * (1 / De + Inclusive) + 1); [m, n] = size(C); isSizeOneByOne = (m == 1) && (n == 1);
if isSizeOneByOne
    C = C * ones(Ny, Nx);
end
E = dt * (u .* (1 / (2* dx)  + dydx .* dTdy));
W = -dt * (u .* (1 / (2* dx)  + dydx .* dTdy));

for i = 1:Ny
%     T(i, 2:end-1) = TDMA(W(i, 1:end-2), C(i, 2:end-1), E(i, 3:end), Q(i, 2:end-1));
    T(i, :) = TDMA(W(i, :), C(i, :), E(i, :), Q(i, :));
end
        
    
end





end