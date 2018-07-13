        function [u,s,v] = svshr_svdsmartc(a,m,n,k)
%
        if (k==1)
%
        [u,s,v] = svds(a,2);
        s=s(1);
        u=u(:,1);
        v=v(:,1);
        return;
    end

        if (m/k > 20)
%
        [u,s,v] = svds(a,k);
        s=diag(s);
        return;
    end

        if (m/k <= 20)
%
        [u,s,v] = svd(a,'econ');
        s=diag(s(1:k,1:k));
%%%        s=diag(s);
        u = u(1:m,1:k);
        v = v(1:n,1:k);
        return;
    end

        end
%
%
%
%
%
