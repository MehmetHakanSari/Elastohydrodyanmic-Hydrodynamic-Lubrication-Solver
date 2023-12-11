function coeff = domain_coeff_function(u)

    if u == 0.0001
        coeff = 10;
    elseif u == 0.00056234
        coeff = 12;
%         coeff = 12;
    elseif u == 0.0031623
%         coeff = 10;
        coeff = 20;
    elseif u == 0.017783
        coeff = 22;
%         coeff = 40;
    elseif u == 0.1
%         coeff = 10;
        coeff = 25;
    elseif u == 1
%         coeff = 10;
        coeff = 30;
    elseif u == 2.5
        coeff = 10;
%         coeff = 30;
    end

end