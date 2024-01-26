% primjena FFT

f = @(x) exp(-x.*x/4);
n = 4;
N = 16;

x = zeros(N,1);
F = zeros(N,1);
y = zeros(N,1);

for k = 0:N-1
    x(k+1,1) = 2*k*pi/N;
    F(k+1,1) = f(x(k+1,1));
end

beta = FFT(F,n)/N;
%{
fprintf('interpolacijski\n');
sb = max(size(beta));
for k = 1:sb
    fprintf('beta(k): %f\n',beta(k));
end
%}
[Ah, Bh] = trig_FFT(F,n);

for k = 0:N-1
    y(k+1,1) = horner(x(k+1,1),beta);
end

err = norm(F-y,'inf');
fprintf('%e\n',err);

fplot(f, [0 6], 'b');
hold on;
plot(x,y,'r-o');
legend('f','y');
hold off;