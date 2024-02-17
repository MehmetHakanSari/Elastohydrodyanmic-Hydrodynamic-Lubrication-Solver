function T_xx_nextstep = T_xx_ADI(T_xx, T_xy, u, v, DuDx, dudy, dx, dy, dydx, dt, De, viscos)
    %below this function create operators for T_xx terms.

    [ny, nx] = size(T_xx);

    T_xx_U = T_xx_nextimestepterm(T_xx, DuDx, De);

    Const = 2 * T_xy .* dudy; %semi-implicit part

    gradTxx = TwoDcentraldiff(T_xx, dx, dy);
    dTxxdy = gradTxx{2};

    Q = (- v .* dTxxdy + Const) * 1/2 * dt + T_xx; %Forcing Term of TDMA.  
    C = 1 - T_xx_U * dt / 2; %Center Diagonal of TDMA.
    E = dt * dydx .* u / (4 * dx); %Upper Diagonal of TDMA.
    W = -dt * dydx .* u / (4 * dx); %Lower Diagonal of TDMA.

    for i = 1:ny
        T_xx(i,2:end-1) = TDMA(W(i,1:end), C(i,2:end-1), E(i,3:end), Q(i,2:end-1));
    end
    
    T_xx_halfstep = T_xx;
    T_xx_U_halfstep = T_xx_nextimestepterm(T_xx_halfstep, DuDx, De);
    gradTxx_halfstep = TwoDcentraldiff(T_xx_halfstep, dx, dy);
    DTxxDx_halfstep = gradTxx_halfstep{1}+ gradTxx_halfstep{2} .* dydx;

    Q = (-u .* DTxxDx_halfstep + Const) * 1/2 * dt + T_xx_halfstep; %Forcing Term of TDMA.
    C = 1 - T_xx_U_halfstep * dt / 2; %Center Diagonal of TDMA.
    N = dt * v ./ (4 * dy); %Upper Diagonal of TDMA.
    S = -dt * v ./ (4 * dy); %Lower Diagonal of TDMA.

    for i = 1:nx
        T_xx(2:end-1,i) = TDMA(S(1:end-2,i), C(2:end-1,i), N(3:end,i), Q(2:end-1,i));
    end

    T_xx_nextstep = T_xx;
end