        function [rlam,rs,cs] = lsnuc_find_rlam(var_ep,pvec,m,n,nsims)
%
%
%        approximates rlam as the operator norm of noise term,
%        depending on variances and sampling probabilities
%

        rs = sqrt(pvec / min(pvec));
        cs = ones(1,n);
        ps = repmat(pvec,1,n);

        rnorms=zeros(1,nsims);
        for i=1:nsims
%
        inds00 = lsnuc_rand_inds(m,n,ps);
        ep00 = diag(sqrt(var_ep))*randn(m,n);
        ep00 = ep00.*inds00;
        ep00 = diag(1./rs)*ep00;
%
        rnorms(i) = norm(ep00);
    end
        rlam=median(rnorms);

        end
%
%
%
%
%
