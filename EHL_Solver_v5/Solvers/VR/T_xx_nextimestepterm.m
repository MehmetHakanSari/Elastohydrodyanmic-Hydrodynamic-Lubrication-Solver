function T_xx_U = T_xx_nextimestepterm(T_xx, DuDx, De)
    T_xx_U = 2 * T_xx .* DuDx - T_xx / De;
end