        function [rlam,cos_out,cos_inn] = svshr_pop2emp_white2(ell,gam,sig)
%
        ell = ell / sig^2;
        [rlam,cos_out,cos_inn] = svshr_pop2emp_white(ell,gam);
        rlam = sig^2 * rlam;

        end
%
%
%
%
%
