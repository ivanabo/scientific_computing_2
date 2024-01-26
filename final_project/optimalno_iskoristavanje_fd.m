% AMR 1.22 Optimal Harvesting
% Cilj je za odabranu funkciju y odrediti tocke s t.d. je maksimizirana
% vrijednost fukcije v koja zadovoljava jednadzbu
% \epsilon v'' + (f(x)-y(x))v' - \gamma v + n(x)y(x) = 0 .............(1a)
% za x in <-\infty, \infty>
% parametri \epsilon i \gamma su zadani
% funkcije f, n i y su zadane, s time da za je 
% y(x) = y- za n(x)<v'(x)
% y(x) = y+ za n(x)>v'(x)
% te vrijedi lim_{abs(x)->\infty}v(x) < \infty
% -------------------------------------------------------------------------
% Pretpostavljamo da postoji samo jedna tocka skoka s te je domena [-L,L]
% Tada imamo rubne uvjete:
% (f(-L)-y-)v'(-L) - \gamma v(-L) + n(-L)y- = 0 ........................(2)
% (f(L)-y+)v'(L) - \gamma v(L) + n(L)y+ = 0 ............................(3)
% Prebacujemo tocku s u tocku x = 0, transformacijom x = x-s
% Tada (1a) prelazi u 
% \epsilon v'' + (f(x+s)-y)v' - \gamma v + n(x+s)y = 0 ................(1b)
% rjesavanje konacnim diferencijama
% parametar s dobivamo jedostavno postupkom pogadanje <-> usporedivanje
%--------------------------------------------------------------------------

function [x,tocka,rjesenje] = optimalno_iskoristavanje_fd


% broj tocaka diskretizacije
N = 11; % uz dvije rubne æe ih biti N+2
% domena [-L,L]
L = 10;
% korak
h = 2*L/(N+1);
% zadani parametri za jednadzbu, epsilon proizvoljan, gamma in <0,1>
epsilon = 0.1;
gamma = 0.4;

% diskretizacija domene sa gore odabranim parametrima
x = zeros(N+2,1);
for i = 1:N+2
    x(i) = -L + (i-1)*h;
end

maks = 0; % pomocna za odabir tocke s

for s = -L:0.01:L;
    % problem cemo svesti na rjesavanje sustava Av=b
    % gdje je A matrica dimenzija (N+2)x(N+2) dobivena aproksimiranjem v'', v'
    % v vektor nepoznanica v(i) = v(x(i)), b vektor desne strane duljine N+2
    A = zeros(N+2,N+2);
    b = zeros(N+2,1);
    % prvi redak dobivamo iz rubnog uvjeta (2) uz v'(-L) = (v(2)-v(1))/2
    A(1,1) = f(x(1))-y(1) - h*gamma;
    A(1,2) = f(x(1))-y(1);
    b(1) = -n(x(1))*y(1);

    for j = 2:N+1
        if j < (N+2)/2
            A(j,j-1) = 1+h*(y(1)-f(x(j)+s)/(2*epsilon));
            A(j,j) = -2-h*h*gamma;
            A(j,j+1) = 1-h*(y(1)+f(x(j)+s))/(2*epsilon);
            b(j) = n(x(j)+s)*y(1);
        elseif j == (N+2)/2
            A(j,j-1) = -1;
            A(j,j) = 2;
            A(j,j+1) = -1;
        else % j > (N+2)/2
            A(j,j-1) = 1+h*(y(2)-f(x(j)+s)/(2*epsilon));
            A(j,j) = -2-h*h*gamma;
            A(j,j+1) = 1-h*(y(2)+f(x(j)+s))/(2*epsilon);
            b(j) = n(x(j)+s)*y(2);
        end
    end
    
    A(N+2,N+1) = -f(x(N+2))+y(2);
    A(N+2,N+2) = f(x(N+2))-y(2)-h*gamma;
    b(N+2) = -n(x(N+2))*y(2);
    
    v = A\b;
    
    if v(end) > maks
        tocka = s;
        rjesenje = v;
    end
end

% -----------------------------------------------------------------------
% pomocne funkcije kojima su definirane zadane funkcije f, n i y
%
    function [z] = f(x)
        z = 1-exp(x);
    end % f
% -----------------------------------------------------------------------
    function [z] = n(x)
        z = exp(x);
    end % n
% -----------------------------------------------------------------------
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
% -----------------------------------------------------------------------
end % opi_fd