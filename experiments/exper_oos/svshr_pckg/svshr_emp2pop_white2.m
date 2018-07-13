        function [ell,cos_out,cos_inn] = svshr_emp2pop_white2(rlam,gam,sig)
%
        rlam = rlam / sig^2;
        [ell,cos_out,cos_inn] = svshr_emp2pop_white(rlam,gam);
        ell = sig^2 * ell;



        end
%
%
%
%
%
