        function rlam = lsnuc_find_rlam(var_ep,ps,m,n,nsims)
%
%
%        approximates rlam as the operator norm of noise term,
%        depending on variances and sampling probabilities
%
        rlam = 0;
        rnorms=zeros(1,nsims);
        for i=1:nsims
%
        inds00 = lsnuc_rand_inds(m,n,ps);
        ep00 = diag(sqrt(var_ep))*randn(m,n);
        ep00 = ep00.*inds00;
        rnorms(i) = norm(ep00);
    end
        rlam=median(rnorms);

        end
%
%
%
%
%
