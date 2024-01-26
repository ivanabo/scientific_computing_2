function [z] = odj_newton(g,df,z0,k)

if abs(g(z0)) < 1e-8
    z = z0;
    return
end

for i = 1:k
    
    z = z0 - g(z0)/df(z0);
    z0 = z;    

end
