        function errs = svshr_errs(ells,cos_inn,...
           cos_out,k)
%
        errs = zeros(1,k);
        for i=1:k
%
        errs(i) = ells(i) * (1 - cos_inn(i)^2 * cos_out(i)^2);
    end
        end
%
%
%
%
%
