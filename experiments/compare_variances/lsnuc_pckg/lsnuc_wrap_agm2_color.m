        function [ww,nmax,whts] = lsnuc_wrap_agm2_color(y,inds,...
            m,n,w0,niter,thresh,rlam,rs,cs)
%
        iis = find(inds == 1);
        len=length(iis);
%
        [ws,nmax,vals,whts] = lsnuc_agm2(y,m,n,iis,len,rlam,w0,niter,thresh,rs,cs);
        ww=ws(:,:,nmax);

        ww = diag(1./rs) * ww;

        end
%
%
%
%
%
