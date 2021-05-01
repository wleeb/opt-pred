        function bedge = lintr_edge_diag(pvec,var_ep,m,n,nsims)
%
%        estimate bulk edge of eigenvalue distribution of
%        the backprojected noise 
%
        rnorms = zeros(1,nsims);

        pvec = repmat(pvec,1,n);
        for i=1:nsims
%
        as = lintr_rand_inds(m,n,pvec);

        ep0 = diag(sqrt(var_ep))*randn(m,n);
        ep0_back = conj(as) .* ep0;
        rnorms(i)=norm(ep0_back)^2 / n;
    end
        bedge = median(rnorms)

        end
%
