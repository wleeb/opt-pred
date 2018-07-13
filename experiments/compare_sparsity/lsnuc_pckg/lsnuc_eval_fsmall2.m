        function fval = lsnuc_eval_fsmall2(a,iis,y,len,m,n,whts)
%
%        evaluates the smooth term in the objective function.
%        y is the observed data, a is the variable matrix
%
        dd = a(iis) .* whts(iis);
        fval = norm(dd - y(iis))^2 / 2;

        end
%
%
%
%
%
