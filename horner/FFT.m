function [b] = FFT(F,n)

N = 2^n;
b = zeros(N,1);

for j = 0:N-1
    b(j+1) = F(rho(j,n)+1);
end

for m = 1:n
    for j = 0:2^(m-1)-1
        temp = complex(0,-2*pi*j/(2^m));
        e = exp(temp);
        
        for q = 0:2^m:N-1
            u = b(q+j+1);
            v = b(q+j+2^(m-1)+1)*e;
            b(q+j+1) = u+v;
            b(q+j+2^(m-1)+1) = u-v;
        end
    end
end

end