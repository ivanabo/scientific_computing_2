function [A,B] = trig_FFT(F,n)

N = 2^n;
M = N/2;

A = zeros(M+1,1);
B = zeros(M,1);

beta = FFT(F,n)./N;

A(1) = 2*beta(1);

for k = 1:M-1
    A(k+1) = beta(k+1) + beta(N-k+1);
    B(k+1) = i*( beta(k+1) - beta(N-k+1) );
end

A(M+1) = 2*beta(M+1);
A = real(A);
B = real(B);

end