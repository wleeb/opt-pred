        function [ell,cos_out,cos_inn] = svshr_emp2pop_white_flip(rlam,gam)
%
%       takes in empirical eigenvalue of \hat{\Sigma}_{Y/sqrt{delta}}, and 
%       returns the corresponding population eigenvalue if possible (if 
%       rlam exceeds the bulk edge)
%
%       Critically, this code assumes gam > 1
%
        [ell,cos_inn,cos_out] = svshr_emp2pop_white(rlam/gam,1/gam);
        ell = ell*gam;

        end
%
%
%
%
%
