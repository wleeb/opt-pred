        function [xs_in,xs_out,errs_hat] = lintr_inout_spike(ys_in,ys_out,...
           m,nin,nout,k,sig)
%
%        computes the optimal in-sample (shrinkage) denoiser and 
%        out-of-sample denoiser for data from spiked model with
%        white noise of specified variance
%
        [xs_in,s_op,s_fr,cos_out,cos_inn,uy,sy,vy,errs] = ...
           svshr_white2(ys_in,m,nin,k,sig);

        ells = s_op.^2;
        gam_in = m / nin;
        [xs_out,errs_hat] = lintr_out_spike(ys_out,uy,ells,m,nout,k,gam_in,sig);

%%%        chk0 = norm(errs_hat - errs')

        end
%
%
%
%
%
