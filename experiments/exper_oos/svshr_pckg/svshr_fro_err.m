        function [err_fr,cos_out,cos_inn] = svshr_fro_err(ells,...
           u,u_hat,v,v_hat,k)
%
        cos_out = zeros(1,k);
        cos_inn = zeros(1,k);

        for i=1:k
%
        cos_out(i) = sum(u(:,i) .* u_hat(:,i));
        cos_inn(i) = sum(v(:,i) .* v_hat(:,i));
    end

        err_fr = svshr_err_fmla(ells,cos_inn,cos_out,k);
        end
%
%
%
%
%
