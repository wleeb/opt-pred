        function val = svshr_mp_eval(t,gam,sig)
%
%        evaluates standard Marchenko-Pastur density with 
%        aspect ratio gam at t
%
        x0 = sig^2*(1-sqrt(gam))^2;
        x1 = sig^2*(1+sqrt(gam))^2;
        val = (x1 - t) * (t - x0);
        val = sqrt(val) / (2*pi*t*gam) / sig^2;

        end
%
%
%
%
%
