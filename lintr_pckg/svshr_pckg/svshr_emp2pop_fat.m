        function [si_hat,cos_out,cos_inn] = svshr_emp2pop_fat(s,m,n,k,i,bedge)
%
%        estimate i^th singular value using the Stieltjes transform
%
        si = s(i);
        rlam = si^2
%
        if (rlam < bedge)
%
        si_hat = 0;
        cos_out=0;
        cos_inn=0;
        return
    end


        [d_hat,d_der,stra,stra_der,sbar,sbar_der] = svshr_integrs_stiel(s,...
           k,m,n,rlam)

%
%        estimate i^th singular value
%
        si_hat = 1/sqrt(d_hat);
%
%
%        estimate i^th inner and outer angles between pop and emp singular
%        vectors
%
        cos_inn = sqrt(sbar / (d_der * si_hat^2))
        cos_out = sqrt(stra / (d_der * si_hat^2))


        end
%
%
%
%
%
