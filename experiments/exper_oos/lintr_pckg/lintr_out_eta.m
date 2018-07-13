        function [eta,err_hat] = lintr_out_eta(ell,gam,sig)
%
%        out-of-sample denoising coefficient for spike strength ell,
%        for white noise with variance sig^2. Also returns estimated
%        error.
%
        ell = ell / sig^2;
        done=1;
        [rlam,cos_out,cos_inn] = svshr_pop2emp_white2(ell,gam,done);
        eta = ell*cos_out^2 / (ell*cos_out^2 + 1);

        aa = 1 + ell*cos_out^2;
        bb = -2*ell*cos_out^2;
        err_hat = ell + eta^2*aa + eta*bb;

        ell = ell*sig^2;
        err_hat = err_hat * sig^2;

        end
%
%
%
%
%
