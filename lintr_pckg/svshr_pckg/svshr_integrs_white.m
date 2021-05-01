        function [d_hat,d_der,stra,stra_der,sbar,sbar_der] = ...
           svshr_integrs_white(z,gam)
%
%       evaluates the Stieltjes transform of white noise with unit
%       covariance with aspect ratio gam at z -- only tested for real z
%
        top = (1 - gam- z) + sqrt( (z - 1 - gam)^2 - 4*gam);
        bot = 2*gam*z;
        stra = top/bot;
%
        top_der = -1 + (z-1-gam) / sqrt( (z - 1 - gam)^2 - 4*gam);
        bot_der = 2*gam;
        stra_der = (bot*top_der - top*bot_der) / bot^2;
%
        sbar = gam*stra - (1-gam)/z;
        sbar_der = gam*stra_der + (1-gam)/z^2;
%
        d_hat = stra*sbar*z;
        d_der = stra_der*sbar*z + stra*sbar_der*z + stra*sbar;


        end
%
%
%
%
%
