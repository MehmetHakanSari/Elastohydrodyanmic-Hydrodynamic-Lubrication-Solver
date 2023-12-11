function [press] = FluidPressISOSliding(nx,mu,h,dhdx,cn,cs,u,pin,poldi,pout,pc,Ux,relax,tollp,maxiter, dx)
        % Inizialization Coeffs
        e = -(6*mu*Ux*dhdx) ./ h.^3;  
        g = e / u;                     
        
        % Inizialization
        press = zeros(1,nx);
        press(1) = pin;
        press(nx) = pout;

        press(2:nx-1) = poldi(2:nx-1); 
        
        for iter = 1:maxiter
            pold = press;
            %Ciclo for per il calcolo della pressione nel corpo
            for i = 2:nx-1
                press(i) = (1-relax)*pold(i) + relax * (cn(i)*press(i+1) + cs(i)*press(i-1) + g(i));
                if press(i) < pc
                    press(i) = pc;
                end
            end
            % Error calculation
            err = max(abs((press-pold)./press));
            if err < tollp
                break;
            end
        end
        
%         flowrate = h'/2 * Ux - h'.^3 / 12 / mu .* OneDcentraldiff(press, dx);

end



