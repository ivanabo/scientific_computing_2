function [sol] = opt_isko

% globalni parametri
epsilon = 0.1;
gama = 0.2;

% diskretizacija, uz dupliciranje tocke x = 0
L = 10;
right = 0:1:L;
left = fliplr(-right);
xinit = [left, right];

% inicijalna aproksimacija rjesenja
yinit = [1; 1; 1];

% postavljanje starta za bvp4c
sol = bvpinit(xinit,yinit);

% poziv funkcije koja vraca aproksimaciju rjesenja
sol = bvp4c(@f,@bc,sol);

len = length(sol.y(3,:));
tocka = sol.y(3,floor(len/2));
figure(1)
plot(sol.x,sol.y(1,:));
legend('v(x)');
title(sprintf('funkcija v, tocka s = %d', tocka));


% -----------------------------------------------------------------------
% pomocne funkcije
%

    function z = F(x)
        % Funkcija 'koeficijenta' brzine rasta dobra - na papiru oznaka f
        z = 1-exp(x);
    end % F

% -----------------------------------------------------------------------

    function z = N(x)
        % Funkcija 'ucinka' iskoristavanja - na papiru oznaka n
        z = exp(x);
    end % N

% -----------------------------------------------------------------------

    function z = w(podrucje)
        % Funkcija 'intenziteta iskoristavanja' - na papiru oznaka y
        switch podrucje
            case 1 % x in [-L 0]
                z = 1;
            case 2 % x in [0 L]
                z = 4;
        end
    end % z

% -----------------------------------------------------------------------

    function dydx = f(x,y,region)
      % Sustav
      dydx = zeros(3,1);
      
      dydx(1) = y(2);
      
      % w ima prekid u tocki x = 0 pa je stoga v'' razlicito definirana
      % u pojednom podrucju domene
      switch region
         case 1    % x in [-L 0]
            dydx(2) = -((F(x+y(3))-w(1))*y(2)-gama*y(1)+N(x+y(3))*w(1))/epsilon;
         case 2    % x in [0 L]
            dydx(2) = -((F(x+y(3))-w(2))*y(2)-gama*y(1)+N(x+y(3))*w(2))/epsilon;
         otherwise
            error('MATLAB:threebvp:BadRegionIndex','Incorrect region index: %d',region);
      end
      
      dydx(3) = 0;
    end % f
% -----------------------------------------------------------------------

   function res = bc(YL,YR)
      % Rubni (i u nutrasnji) uvjeti
      res = [ YR(2,1) - N(YR(3,1))             % v'(0) = n(s(0)) in region 1
              YL(2,2) - N(YR(3,2))             % v'(0) = n(s(0)) in region 2
              YR(1,1) - YL(1,2)                % neprekidnost v(x) at x = 0
              YR(2,1) - YL(2,2)                % neprekidnost v'(x) at x = 0
              (F(-L)-w(1))*YL(2,1)-gama*YL(1,1)+N(-L)*w(1)   % uvjet 5 na papiru
              (F(L)-w(2))*YR(2,2)-gama*YR(1,2)+N(L)*w(2)];   % uvjet 6 na papiru
   end % bc
% -----------------------------------------------------------------------


end  % opt_isko

