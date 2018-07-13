        function [ww,nmax,rlam,whts] = lsnuc_wrap_agm2(y,inds,m,n,w0,niter,thresh,pvec,sig)
%
%        code for convenience: assumes sampling varies only across rows, with
%        probabilities given by vector pvec; and noise is white with standard 
%        deviation sig. Correct weight rlam is determined automatically, as the operator
%        norm of the pure noise term.
%
        iis = find(inds == 1);
        len=length(iis);
%
        rs = sqrt(pvec);
        rs = rs / min(sqrt(pvec));
        rrr=min(sqrt(pvec));
        rlam = sig*(sqrt(n) + sqrt(m)) * rrr;
%
        cs = ones(n,1);
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
