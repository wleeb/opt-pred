        function err_fr = svshr_err_fmla(ells,cos_inn,...
           cos_out,k)
%
        err_fr = 0;
        for i=1:k
%
        err_fr = err_fr + ells(i) * (1 - cos_inn(i)^2 * cos_out(i)^2);
    end
        end
%
%
%
%
%
