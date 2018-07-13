        function w = svshr_joint_basis(u,v,m,k)
%
%
%       make the joint basis w
%
        g = zeros(m,m);
        g(1:m, 1:k) = u(1:m, 1:k);
        g(1:m, k+1:2*k) = v(1:m, 1:k);
        g(1:m, 2*k+1:m) = v(1:m, k+1:m-k);

        w = svshr_gs_dumb(g,m);

        chk0 = norm(w(1:m,1:k) - u(1:m, 1:k))

        iperm = 1:m;
        for j=1:k
%
        iperm(2*j-1) = j;
        iperm(2*j) = k+j;
    end


%%%        iperm
        w = w(1:m,iperm);


        end
%
%
%
%
%
