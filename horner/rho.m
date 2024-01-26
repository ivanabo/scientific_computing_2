function [z] = rho(x, n)

z = dec2bin(x,n);
z = fliplr(z);
z = bin2dec(z);

end