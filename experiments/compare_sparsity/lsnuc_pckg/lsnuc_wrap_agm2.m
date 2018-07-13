        function [ww,nmax,whts] = lsnuc_wrap_agm2(y,inds,m,n,rlam,w0,niter,thresh)
%
        iis = find(inds == 1);
        len=length(iis);

        rs=mean(inds,2);
        rs = rs / min(rs);

%%%        cs=mean(inds,1);
%%%        cs = cs / min(cs);
        cs = ones(n,1);

        [ws,nmax,vals,whts] = lsnuc_agm2(y,m,n,iis,len,rlam,w0,niter,thresh^2,rs,cs);
        ww=ws(:,:,nmax);

        end
%
%
%
%
%
