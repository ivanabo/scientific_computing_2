function [v,s] = opt_isko_fd_ode_s1

% broj tocaka diskretizacije
N = 11; % uz dvije rubne æe ih biti N+2
% domena [-L,L]
Ld = 10;

% zadani parametri za jednadzbu; epsilon proizvoljan, gamma in <0,1>
epsilon = 0.1;
gamma = 0.4;

h = 1/(N+1);

t = zeros(N+2,1);
for i = 1:N+2
    t(i) = (i-1)*h;
end

% z = zeros((N+2)*6,1);
L = zeros((N+2)*6,(N+2)*6);
b = zeros((N+2)*6,1);

% rubni uvjeti
L(1,end-5) = 1;
L(2,1) = -gamma;
L(2,2) = f(t(1)-y(1));
    b(2) = n(t(1))*y(1);
L(3,5) = 1;
L(4,end-2) = f(t(end)-y(2));
L(4,end-3) = -gamma;
    b(4) = n(t(end))*y(2);
L(5,4) = -1;
L(5,end-6) = 1;
L(6,5) = -1;
L(6,end-5) = 1;

% pomocne varijable za vrijednosti elemenata koji se ponavljaju
kvad = -h*gamma/epsilon;
zvj = -(1+h*y(1)/epsilon);
zvjc = -(1+h*y(2)/epsilon);

% ---------------------------------------------------
for j = 1:N+1
    L(j*6+1,j*6+1) = -1;    L(j*6+1,j*6+2) = -h;    
    L(j*6+2,j*6+1) = kvad;  L(j*6+2,j*6+2) = zvj;   
    L(j*6+3,j*6+3) = -1;                            
    L(j*6+4,j*6+4) = -1;    L(j*6+4,j*6+5) = -h;    
    L(j*6+5,j*6+4) = kvad;  L(j*6+5,j*6+5) = zvjc;  
    L(j*6+6,j*6+6) = -1;                            
    if j < N+1
        L(j*6+1,(j+1)*6+1) = 1;
        L(j*6+2,(j+1)*6+2) = 1;
        L(j*6+3,(j+1)*6+3) = 1;
        L(j*6+4,(j+1)*6+4) = 1;
        L(j*6+5,(j+1)*6+5) = 1;
        L(j*6+6,(j+1)*6+6) = 1;
    end
end
% -----------------------------------------------------

R = @nR;
A = @nA;
J = @Jacobian;

% --------------- iteracije Newtonove metode ----------------------------
rj0 = ones((N+2)*6,1);

for m = 1:2
  rj = rj0 - J(rj0)\A(rj0);
  rj0 = rj;
end

% ------------- premjestanje rjesenja u pogodan oblik -------------------
v = zeros((N+2)*2,1);
s = zeros(N+2,1);
for j = 1:N+2
    v(j) = rj(1+(j-1)*6);
    v(N+2+j) = rj(4+(j-1)*6);
    s(j) = rj(3+(j-1)*6);
end

% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% pomocne funkcije kojima su definirane zadane funkcije f, n i y
% -----------------------------------------------------------------------
%
    function [z] = f(x)
        z = 1-exp(x);
    end % f
% -----------------------------------------------------------------------
%
    function [z] = n(x)
        z = exp(x);
    end % n
% -----------------------------------------------------------------------
%
    function [z] = y(x)
        switch x
            case 1 % y-
                z = 1;
            case 2 % y+
                z = 4;
            otherwise
                display('Error in value of y(x)');
        end
    end % y
% ----------- derivacija funkcije f -------------------------------------
%
    function [z] = dfdz(x)
        z = -exp(x);
    end % dfdz
% ----------- derivacija funkcije n -------------------------------------
%
    function [z] = dndz(x)
        z = exp(x);
    end % dndz
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% pomocna funkcija kojom je zadan R dio sustava
    function [z] = nR(x)
        z = zeros((N+2)*6,1);
        
        z(1) = -n(x(end-4));
        z(3) = -n(x(6));
        
        for k = 1:N+1
            z(k*6+2) = h*((f(Ld*t(k)-Ld+x(k*6+3)))*x(k*6+2)+n(Ld*t(k)-Ld+x(k*6+3))*y(1))/epsilon;
            z(k*6+5) = h*((f(Ld*t(k)-Ld+x(k*6+6)))*x(k*6+5)+n(Ld*t(k)-Ld+x(k*6+6))*y(2))/epsilon;
        end
    end % nR
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% pomocna funkcija kojom je zadan cijeli sustav A(z)
    function [u] = nA(z)
        u = L*z + R(z) + b;
    end % nA
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% pomocna funkcija kojom je zadan Jacobian
    function [M] = Jacobian(x)
        dRdz = zeros((N+2)*6,(N+2)*6);
        dRdz(1,end-3) = -dndz(x(end-3));
        dRdz(2,6) = -dndz(x(6));
        for l = 1:N+1
            dRdz(l*6+2,l*6+2) = h*(f(Ld*t(l)-Ld+x(l*6+3)))/epsilon;
            dRdz(l*6+2,l*6+3) = h*(dfdz(Ld*t(l)-Ld+x(l*6+3))+y(1)*dndz(Ld*t(l)-Ld+x(l*6+3)))/epsilon;
            dRdz(l*6+5,l*6+5) = h*(f(Ld*t(l)-Ld+x(l*6+6)))/epsilon;
            dRdz(l*6+5,l*6+6) = h*(dfdz(Ld*t(l)-Ld+x(l*6+6))+y(2)*dndz(Ld*t(l)-Ld+x(l*6+6)))/epsilon;
        end
        M = L + dRdz;
    end % Jacobian
% -----------------------------------------------------------------------

end % opt_isko_fd_ode_s1