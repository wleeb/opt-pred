        function dinv = svshr_dinv_white(x,gam)
%
%       computes inverse D-transform at x, for aspect ratio gam,
%       for white noise (unity variance)
%
        dinv = (1 + (1/x)) * (1 + gam*x);
        end
%
%
%
%
%
