function T_yy_nextstep = T_yy_ADI(T_yy, T_xy, u, v, dx, dy, DvDx, dvdy, dydx, dt, De, viscos)
    %below this function create operators for T_yy terms.

    [ny, nx] = size(T_yy);

    T_yy_U = T_yy_nextimestepterm(T_yy, dvdy, De);

    Const = 2 * T_xy .* DvDx + 2 * viscos * dvdy; %semi-implicit part

    gradTyy = TwoDcentraldiff(T_yy, dx, dy);
    dTyydy = gradTyy{2};


    Q = (- v .* dTyydy + Const) * 1/2 * dt + T_yy; %Forcing Term of TDMA.
    C = 1 - T_yy_U * dt / 2; %Center Diagonal of TDMA.
    E = dt * dydx .* u / (4 * dx); %Upper Diagonal of TDMA.
    W = -dt * dydx .* u / (4 * dx); %Lower Diagonal of TDMA.
    
    for i = 1:ny
        T_yy(i,2:end-1) = TDMA(W(i,1:end), C(i,2:end-1), E(i,3:end), Q(i,2:end-1));
    end
    T_yy_halfstep = T_yy;
    T_yy_U_halfstep = T_yy_nextimestepterm(T_yy_halfstep, dvdy, De);
    gradTyy_halfstep = TwoDcentraldiff(T_yy_halfstep, dx, dy);
    DTyyDx_halfstep = gradTyy_halfstep{1} + gradTyy_halfstep{2} .* dydx;

    Q = (-u .* DTyyDx_halfstep + Const) * 1/2 * dt + T_yy_halfstep; %Forcing Term of TDMA.
    C = 1 - T_yy_U_halfstep * dt / 2; %Center Diagonal of TDMA.
    N = dt * v ./ (4 * dy); %Upper Diagonal of TDMA.
    S = -dt * v ./ (4 * dy); %Lower Diagonal of TDMA.

    for i = 1:nx
        T_yy(2:end-1,i) = TDMA(S(1:end-2,i), C(2:end-1,i), N(3:end,i), Q(2:end-1,i));
    end
    T_yy_nextstep = T_yy;
end