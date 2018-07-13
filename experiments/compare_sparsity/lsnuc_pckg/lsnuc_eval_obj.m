        function [val,fval,rnuc] = lsnuc_eval_obj(a,iis,len,y,m,n,rlam)
%
%        evaluates the objective function of the optimization; squared frobenius
%        loss on observed entries plus rlam times nuclear norm.
%
        val=0;
        rnuc=0;
        fval = lsnuc_eval_fsmall(a,iis,y,len,m,n,rlam);

        k=min(m,n);
        [ua,sa,va] = lsnuc_svdsmart(a,m,n,k);

%%%        chk0 = norm(ua*diag(sa)*va' - a,'fro')

        rnuc = sum(sa);
        val = fval + rlam*rnuc;

        end
%
%
%
%
%
