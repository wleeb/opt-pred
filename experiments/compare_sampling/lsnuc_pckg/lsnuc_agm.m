        function [ws,nmax,vals] = lsnuc_agm(y,m,n,iis,len,rlam,w0,niter,thresh)
%
%
%                            description:
%
%   This code minimizes the convex objective
%
%                     || P(A - Y) ||_F^2 + L * ||A||_*
%
%   where P projects an m-by-n matrix onto the coordinates iis, defined columnwise.
%   It uses the accelerated gradient method defined in the paper ``An Accelerated 
%   Gradient Method for Trace Norm Minimization'' by Ji and Ye.  The error on the
%   k^th iteration decays like 1/k^2.
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
%   rlam - the weight placed on the nuclear norm term in the optimization
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
        alps=zeros(niter,1);
        alps(1)=1;

        ws=zeros(m,n,niter);
        zs=zeros(m,n,niter);

        z=w0;
        w=w0;

        [vals(1),fvals(1),rnucs(1)] = lsnuc_eval_obj(w0,iis,len,y,m,n,rlam);

        ws(:,:,1)=w0;
        zs(:,:,1)=w0;

        ell=1
        nmax=niter;

        vvs=zeros(1,niter);
        uus=zeros(1,niter);

        for ijk=2:niter
%
        alps(ijk) = 1 + sqrt(1 + 4*alps(ijk-1)^2);
        alps(ijk) = alps(ijk)/2;

        w = lsnuc_prox_min(z,ell,rlam,y,iis,len,m,n);
        ws(:,:,ijk) = w;

        coef = (alps(ijk-1) - 1) / alps(ijk);
        z = w + coef*(w - ws(:,:,ijk-1));
        zs(:,:,ijk) = z;
        coefs(ijk) = coef;
        [vals(ijk),fvals(ijk),rnucs(ijk)] = lsnuc_eval_obj(w,iis,len,y,m,n,rlam);

%%%        dif = norm(ws(:,:,ijk) - ws(:,:,ijk-1),'fro') / norm(ws(:,:,ijk),'fro');
        dif2 = abs(vals(ijk) - vals(ijk-1)) / vals(ijk);

        if (dif2 < thresh)
           nmax=ijk;
           break;
        end

    end

        ws=ws(:,:,1:nmax);
        vals=vals(1:nmax);


        end
%
%
%
%
%
