        function u = svshr_gs_dumb(g,m)
%
        u = zeros(m,m);
%
        u(1:m,1) = g(1:m,1) / norm(g(1:m,1));

        for ijk = 2:m
%
        vec = g(1:m,ijk);

        for ii = 1:ijk-1
%
        vec2 = u(1:m,ii);

        prod = sum(vec2 .* vec);
        vec = vec - prod*vec2;

        prod = sum(vec2 .* vec);
        vec = vec - prod*vec2;

%%%        chk0 = sum(vec .* vec2)

        u(1:m,ijk) = vec/norm(vec);
    end
    end
%
        end
%
%
%
%
%
