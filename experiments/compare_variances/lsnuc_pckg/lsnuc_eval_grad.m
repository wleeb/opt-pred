        function grad = lsnuc_eval_grad(a,iis,y,len,m,n)
%
%        evaluates gradient of little f (fsmall), the smooth term
%        in the objective function (the Frobenius norm)
%        y is the observed data, a is the variable matrix
%
        grad = zeros(m,n);
        grad(iis) = a(iis) - y(iis);

        end
%
%
%
%
%
