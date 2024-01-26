function [y] = gen_horner_trig(x, a, b)

N = max(size(a));

B1 = 0;
C1 = 0;
B0 = a(N);
C0 = b(N);
alpha = 2*cos(x);

for k = (N-1):-1:1
    B2 = B1;
    C2 = C1;
    B1 = B0;
    C1 = C0;
    B0 = a(k) + alpha*B1 - B2;
    C0 = b(k) + alpha*C1 - C2;
end

y = B0 - B1*0.5*alpha + C1*sin(x);

end