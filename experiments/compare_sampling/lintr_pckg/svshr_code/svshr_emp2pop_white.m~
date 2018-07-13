        function [ell,cos_out,cos_inn] = svshr_emp2pop_white(rlam,gam)
%
%        computes the population spike and angles, given empirical spike
%        where noise is assumed white
%
        btop = (1+sqrt(gam))^2;
        if (rlam <= btop);
%
        ell=0;
        cos_out=0;
        cos_inn=0;
        return;
    end

        [d_hat,d_der,stra,stra_der,sbar,sbar_der] = ...
           svshr_integrs_white(rlam,gam);
%
        ell = 1/d_hat;
        cos_out = sqrt(stra / (d_der * ell));
        cos_inn = sqrt(sbar / (d_der * ell));


        return
%
%       alternative formulas, explicit in terms of ell and gam
%
        yy = rlam - 1 - gam;
        ell2 = .5*(yy + sqrt(yy^2 - 4*gam));
        chk0 = ell-ell2

        cos_out2 = (1 - gam / ell^2) / (1 + gam / ell);
        cos_out2 = sqrt(cos_out2);
        chk0 = cos_out-cos_out2

        cos_inn2 = (1 - gam / ell^2) / (1 + 1 / ell);
        cos_inn2 = sqrt(cos_inn2);
        chk0 = cos_inn-cos_inn2

        end
%
%
%
%
%
