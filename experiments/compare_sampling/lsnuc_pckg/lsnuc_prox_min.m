        function wmin = lsnuc_prox_min(w1,tt,rlam,y,iis,len,m,n)
%
        grad = lsnuc_eval_grad(w1,iis,y,len,m,n);
        bb = w1 - grad / tt;
        thr = rlam / tt;
        wmin = lsnuc_softhresh(bb,m,n,thr);


        end
%
%
%
%
%
