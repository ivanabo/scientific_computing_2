function [y] = horner(x, beta)

size_beta = max(size(beta));

epsilon = exp(i*x);
y = beta(end);

for j = size_beta-2:-1:0
    y = y*epsilon + beta(j+1);
end

end