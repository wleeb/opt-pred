        function [xs_in,xs_out,errs_hat,wvals] = lintr_inout(ys_in,ys_out,...
           as_in,as_out,m,nin,nout,k,var_ep)
%
%
%        whiten in-sample data 
%
        as2_in = conj(as_in).*as_in;
        as2_mean = mean(as2_in,2);
        [ys4_in,wvals] = lintr_gen2whi2(ys_in,as_in,m,nin,k,var_ep,as2_mean);

        winv = 1./wvals;
%
%        whiten out-of-sample data
%
        [ys4_out,wvals2] = lintr_gen2whi2(ys_out,as_out,m,nout,...
           k,var_ep,as2_mean);

%
%        apply predictor to both in-sample and out-of-sample, and unwhiten
%
        done=1;
        [xs_in,xs_out,errs_hat] = lintr_inout_spike(ys4_in,ys4_out,m,nin,nout,...
            k,done);

        xs_in = diag(winv) * xs_in;
        xs_out = diag(winv) * xs_out;

        end
%
