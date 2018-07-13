        function [u,s,v] = lsnuc_svdsmart(a,m,n,k)
%
        if (min(m,n)/k > 20)
%
        [u,s,v] = svds(a,k);
        s=diag(s);
    end

        if (min(m,n)/k <= 20)
%
        [u,s,v] = svd(a,'econ');
        s=diag(s(1:k,1:k));
%%%        s=diag(s);
        u = u(1:m,1:k);
        v = v(1:n,1:k);
    end

        end
%
%
%
%
%
