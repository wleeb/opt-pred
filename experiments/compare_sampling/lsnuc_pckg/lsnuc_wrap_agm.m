        function [ww,nmax,rlam] = lsnuc_wrap_agm(y,inds,m,n,w0,niter,thresh,sig)
%
        iis = find(inds == 1);
        len=length(iis);

        pp = len / (m*n);

        rlam = sig*(sqrt(n) + sqrt(m))*sqrt(pp);
        [ws,nmax,vals] = lsnuc_agm(y,m,n,iis,len,rlam,w0,niter,thresh);
        ww=ws(:,:,nmax);


        end
%
%
%
%
%
