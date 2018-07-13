        function [ws,nmax,vals] = lsnuc_egm(y,m,n,iis,len,rlam,...
           w0,niter,thresh)
%
%
%                            description:
%
%   This code minimizes the convex objective
%
%                     || P(A - Y) ||_F^2 + L * ||A||_*
%
%   where P project an m-by-n matrix onto the coordinates iis, defined columnwise.
%   It uses the extended gradient method defined in the paper ``An Accelerated 
%   Gradient Method for Trace Norm Minimization'' by Ji and Ye. The error on the
%   k^th iteration decays like 1/k.
%
%
%                           input parameters:
%
%   y - the m-by-n data matrix of observations, with zeros filled in for
%      the unobserved entries
%
%   m,n - the dimensions
%
%   iis - a vector of length len, containing the indices of the observed values,
%      stored columnwise.
%
%   len - the number of observed entries
%
%   rlam - the weight place on the nuclear norm term in the optimization
%
%   ell - a parameter defining the Lipschitz constant of the smooth part. Should
%      be greater than or equal to 2.
%
%   w0 - an m-by-n matrix. The initial guess of the optimization.
%
%   niter - the maximum number of iterations
%
%   thresh - the precision. When the norm between consecutive terms is below
%      thresh, the procedure terminates
%
%
%                        output parameters:
%
%   ws - the m-by-n-by-nmax cube of matrices in the iteration, starting with w0
%
%   nmax - the number of iterations until convergence
%
%   vals - the values of the objective function (squared error plus nuclear norm) at
%      the matrices in ws
%
%

        ws=zeros(m,n,niter);
        vals=zeros(niter,1);
        fvals=zeros(niter,1);
        rnucs=zeros(niter,1);
%
        w=w0;
        ws(:,:,1) = w;
        [vals(1),fvals(1),rnucs(1)] = lsnuc_eval_obj(w,iis,len,y,m,n,rlam);
        
        ell=1;

        nmax=niter;
        for ijk=2:niter
%
        w = lsnuc_prox_min(w,ell,rlam,y,iis,len,m,n);
        ws(:,:,ijk) = w;
        [vals(ijk),fvals(ijk),rnucs(ijk)] = lsnuc_eval_obj(w,iis,len,y,m,n,rlam);

        dif = norm(ws(:,:,ijk) - ws(:,:,ijk-1),'fro');
        if (dif < thresh)
           nmax=ijk;
           break;
        end

    end

        ws=ws(:,:,1:nmax);
        vals = vals(1:nmax);

        end
%
%
%
%
%
