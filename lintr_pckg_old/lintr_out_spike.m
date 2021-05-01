        function [xs_hat,errs_hat] = lintr_out_spike(ys_out,uy_in,ells,...
           m,n_out,k,gam,sig)
%
%        uses PCs from in-sample data to denoise out-of-sample data in
%        standard spiked model with istropic noise
%
        etas = zeros(k,1);
        errs_hat = zeros(k,1);

        gam_out = m / n_out;

        for i=1:k
%
        [etas(i),errs_hat(i)] = lintr_out_eta(ells(i),gam,sig);
    end

        prods = ys_out'*uy_in;
        xs_hat = uy_in*diag(etas)*prods';


        return

%
%        explicit formula, with loops around out-of-sample data and
%        in-sample singular vectors
%
        xs_hat2 = zeros(m,n_out);
        for i=1:n_out

        for j=1:k
%
        p_ij = sum(ys_out(:,i) .* conj(uy_in(:,j)));
        xs_hat2(:,i) = xs_hat2(:,i) + p_ij * etas(j) * uy_in(:,j);
    end
    end

%%%        xs_hat(:,1:2)

        chk0 = norm(xs_hat2 - xs_hat)

%%%        stopnow


        end
%
%
%
%
%
