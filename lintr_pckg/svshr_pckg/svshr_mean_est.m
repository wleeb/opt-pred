        function [est_mean,ds_mean] = svshr_mean_est(ys,ds,m,n)
%
        est_mean = zeros(m,1);
        ds_mean = zeros(m,1);
%
        for i = 1:m
%
        ds_sum=0;
        for j=1:n
%
        est_mean(i) = est_mean(i) + ys(i,j);
        ds_sum = ds_sum + ds(i,j);
    end

        est_mean(i) = est_mean(i)/ds_sum;
        ds_mean(i) = ds_sum / n;
    end

        return

        est_mean2 = sum(ys,2) ./ sum(inds,2);
        chk0 = norm(est_mean - est_mean2)

        end
%
%
%
%
%
