        function [d_hat,d_der,stra,stra_der,sbar,sbar_der] = ...
           svshr_integrs_stiel(s,k,m,n,rlam)
%
%        computes empirical estimates of the stieljes and D-transforms 
%        at the value rlam, with empirical singular values s (already
%        normalized by dimension). Critically, code assumes that m < n
%
        stra = 0
        stra_der = 0

        for i=k+1:m
%
        stra = stra + 1 / (s(i)^2 - rlam);
        stra_der = stra_der + 1 / (rlam - s(i)^2)^2;        
    end

        stra = stra/(m-k);
        stra_der = stra_der/(m-k);

        sbar = (m-k)*stra/(n-k) - (n-m)/rlam/(n-k);
        sbar_der = (m-k)*stra_der/(n-k) + (n-m)/rlam^2/(n-k);

        d_hat = stra*sbar*rlam;
        d_der = stra_der*sbar*rlam + stra*sbar_der*rlam + stra*sbar;


        end
%
%
%
%
%
