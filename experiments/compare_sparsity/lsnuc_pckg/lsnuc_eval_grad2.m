        function grad = lsnuc_eval_grad2(a,iis,y,len,m,n,whts)
%
%        evaluates gradient of little f (fsmall), the smooth term
%        in the objective function (the Frobenius norm)
%        y is the observed data, a is the variable matrix
%
        grad = zeros(m,n);

        dd=whts(iis).*a(iis);
        grad(iis) = (dd - y(iis)) .*  whts(iis);

%%%        grad(iis) = 2*(a(iis) - y(iis));

        end
%
%
%
%
%
