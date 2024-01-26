% primjena generaliziranog Hornerovog algoritma

f = @(x) exp(-x.*x/4);
n = 4;
N = 16;

x = zeros(N,1);
F = zeros(N,1);
y = zeros(N,1);

for k = 0:N-1
    x(k+1) = 2*k*pi/N;
    F(k+1) = f(x(k+1));
end

beta = FFT(F,n)./N;
[Ah,Bh] = trig_FFT(F,n);

for k = 0:N-1
	y(k+1) = horner(x(k+1),beta);
end

Ah(1) = Ah(1)/2;
Ah(end) = Ah(end)/2;
Bh = [Bh; 0];

tp = @(x) gen_horner_trig(x, Ah, Bh);

fplot(f, [0 2*pi], 'm');
hold on;
plot(x, y, 'b*');
fplot(tp, [0 2*pi], 'g');
legend('funkcija f','interpolacijske tocke','trigonometrijski polinom');
hold off;