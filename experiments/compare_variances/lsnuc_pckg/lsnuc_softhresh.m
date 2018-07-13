        function a_thr = lsnuc_softhresh(a,m,n,thr)
%
        k_thr = min(m,n);
        k=min(m,n);
        [ua,sa,va] = lsnuc_svdsmart(a,m,n,k);
%%%        chk0 = norm(ua*diag(sa)*va' - a,'fro')

        s_thr = 0*sa;

        s_thr = max(sa-thr,0);
        k_thr = sum(s_thr>0);

        s_thr=s_thr(1:k_thr);
        a_thr = ua(:,1:k_thr) * diag(s_thr) * va(:,1:k_thr)';



        return


        for i=1:k
%
        s_thr(i) = sa(i) - thr;

        if (s_thr(i) < 0)
           s_thr(i) = 0;
           k_thr = i - 1;
           break;
        end

    end
%
        s_thr = s_thr(1:k_thr);
        a_thr = ua(:,1:k_thr) * diag(s_thr) * va(:,1:k_thr)';

        end
