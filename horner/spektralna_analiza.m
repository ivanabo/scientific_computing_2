% spektralna analiza

phi1 = 20;
phi2 = 50;
phi3 = 70;

g = @(t) 0.2*cos(2*pi*phi1*t) + 0.35*sin(2*pi*phi2*t) + 0.3*sin(2*pi*phi3*t);

N = 512;
n = 9;
T = 2;
M = N/2;
delt = T/N;
delphi = 1/T;
epsilon = normrnd(0,0.5,N,1);

for k = 1:N
    t(k) = delt * (k-1);
    y(k) = g(t(k));
    y_ss(k) = y(k) + epsilon(k);
end

for k = 1:M+1
    phi(k) = (k-1)/T;
end

beta = FFT(y,n);
beta_ss = FFT(y_ss,n);

mod_beta = abs(beta)/N;
mod_beta_ss = abs(beta_ss)/N;

[A,B] = trig_FFT(y,n);
B(end+1) = 0;
[A_ss,B_ss] = trig_FFT(y_ss,n);
B_ss(end+1) = 0;

figure(1)
plot(phi,mod_beta(1:M+1),'b');
title('toèke (phi(k),|beta(k)|) bez šuma');

figure(2)
plot(phi,mod_beta_ss(1:M+1),'b');
title('toèke (phi(k),|beta(k)|) sa šumom');

figure(3)
plot(phi,A(1:M+1));
title('toèke (phi(k),A(k)) bez šuma');

figure(4)
plot(phi,B(1:M+1));
title('toèke (phi(k),B(k)) bez šuma');

figure(5)
plot(phi,A_ss(1:M+1));
title('toèke (phi(k),A(k)) sa šumom');

figure(6)
plot(phi,B_ss(1:M+1));
title('toèke (phi(k),B(k)) sa šumom');