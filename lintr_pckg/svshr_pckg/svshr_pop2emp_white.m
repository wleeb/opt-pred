        function [rlam,cos_out,cos_inn] = svshr_pop2emp_white(ell,gam)
%
%       takes the population eigenvalue ell in unreduced noise model, and 
%       returns the asymptotic empirical eigenvalue of the unbiased
%       covariance estimator \hat{\Sigma}_s
%

        if (ell <= sqrt(gam))
%
        rlam = (1 + sqrt(gam))^2;
        cos_out = 0;
        cos_inn = 0;
        return
    end


        rlam = (ell + 1) * (1 + gam/ell);

        cos_out = (1 - gam / ell^2) / (1 + gam / ell);
        cos_out = sqrt(cos_out);
%
        cos_inn = (1 - gam / ell^2) / (1 + 1 / ell);
        cos_inn = sqrt(cos_inn);

        end
%
%
%
%
%
