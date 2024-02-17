function T_xy_nextstep = T_xy_ADI(T_xx, T_xy, T_yy, u, v, DuDx, dudy, DvDx, dvdy, dydx, dt, De, viscos, dx, dy)
    %below this function create operators for T_xy terms.

    [ny, nx] = size(T_xy);

    T_xy_U = T_xy_nextimestepterm(T_xy, De);

    Const = T_xx .* dvdy + T_yy .* dudy + viscos * dudy; %semi-implicit part

    gradTxy = TwoDcentraldiff(T_xy, dx, dy);
    dTxydx = gradTxy{1};
    dTxydy = gradTxy{2};

    Q = (- v .* dTxydy + Const) * 1/2 * dt + T_xy; %Forcing Term of TDMA.
    C = 1 - T_xy_U * dt / 2; %Center Diagonal of TDMA.
    E = dt * dydx .* u / (4 * dx); %Upper Diagonal of TDMA.
    W = -dt * dydx .* u / (4 * dx); %Lower Diagonal of TDMA.

    for i = 1:ny
        T_xy(i,2:end-1) = TDMA(W(i,1:end), C(i,2:end-1), E(i,3:end), Q(i,2:end-1));
    end
    T_xy_halfstep = T_xy;
    T_xy_U_halfstep = T_xy_nextimestepterm(T_xy_halfstep, De);
    gradTxy_halfstep = TwoDcentraldiff(T_xy_halfstep, dx, dy);
    DTxyDx_halfstep = gradTxy_halfstep{1} + gradTxy_halfstep{2} .* dydx;

    Q = (-u .* DTxyDx_halfstep + Const) * 1/2 * dt + T_xy_halfstep; %Forcing Term of TDMA.
    C = 1 - T_xy_U_halfstep * dt / 2; %Center Diagonal of TDMA.
    N = dt * v ./ (4 * dy); %Upper Diagonal of TDMA.
    S = -dt * v ./ (4 * dy); %Lower Diagonal of TD

    for i = 1:nx
        T_xy(2:end-1,i) = TDMA(S(1:end-2,i), C(2:end-1,i), N(3:end,i), Q(2:end-1,i));
    end
    T_xy_nextstep = T_xy;
end