        function s_hat = svshr_emp2frob(s,m,n,k,i)
%
%        estimates i^th singular value using RRN's formula
%

        if (m > n)
%
        m2=m;
        m=n;
        n=m2;
    end

        s_tail = s(k+1:m);
        si = s(i);
%
        n1=n-k;
        m1=m-k;

        t1 = sum(1./(si^2 - s_tail.^2))*si/n1 + (n1-m1)/si/n1;
        t2 = sum(1./(si^2 - s_tail.^2))*si/m1;
%
        t3 = sum(1./(si^2 - s_tail.^2) - 2*si^2./(si^2 - s_tail.^2).^2)/m1;
        t4 = sum(1./(si^2 - s_tail.^2) - 2*si^2./(si^2 - s_tail.^2).^2)/n1;
        t4 = t4 - 2*(n1-m1)/si^2/n1 + (n1-m1)/si^2/n1;
%
        dbot = t1*t3 + t2*t4;
        dtop = t1*t2;
%
        s_hat = -2*dtop/dbot;



        end
%
%
%
%
%
