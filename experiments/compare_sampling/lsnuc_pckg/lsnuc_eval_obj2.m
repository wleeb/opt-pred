        function [val,fval,rnuc] = lsnuc_eval_obj2(a,iis,len,y,m,n,rlam,whts)
%
%        evaluates the objective function of the optimization; squared frobenius
%        loss on observed entries plus rlam times nuclear norm.
%
        val=0;
        rnuc=0;
        fval = lsnuc_eval_fsmall2(a,iis,y,len,m,n,whts);

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
