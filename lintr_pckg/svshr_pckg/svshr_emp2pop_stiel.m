        function [si_hat,cos_out,cos_inn] = svshr_emp2pop_stiel(s,m,n,k,i,bedge)
%
        if (m <= n)
%
        [si_hat,cos_out,cos_inn] = svshr_emp2pop_fat(s,m,n,k,i,bedge);
    end
%
        if (m > n)
%
        [si_hat,cos_inn,cos_out] = svshr_emp2pop_fat(s,n,m,k,i,bedge);
    end

        end
%
%
%
%
%
