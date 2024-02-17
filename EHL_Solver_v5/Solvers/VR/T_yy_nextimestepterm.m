function T_yy_U = T_yy_nextimestepterm(T_yy, dvdy, De)
    T_yy_U = 2 * T_yy .* dvdy - T_yy / De;
end
