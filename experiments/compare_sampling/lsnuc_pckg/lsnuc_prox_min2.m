        function wmin = lsnuc_prox_min2(w1,tt,rlam,y,iis,len,m,n,whts)
%
        grad = lsnuc_eval_grad2(w1,iis,y,len,m,n,whts);
        bb = w1 - grad / tt;
        thr = rlam / tt;
        wmin = lsnuc_softhresh(bb,m,n,thr);


        end
%
%
%
%
%
