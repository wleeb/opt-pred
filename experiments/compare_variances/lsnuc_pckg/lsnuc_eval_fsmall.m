        function fval = lsnuc_eval_fsmall(a,iis,y,len,m,n,rlam)
%
%        evaluates the smooth term in the objective function.
%        y is the observed data, a is the variable matrix
%
        fval = norm(a(iis) - y(iis))^2 / 2;

        end
%
%
%
%
%
