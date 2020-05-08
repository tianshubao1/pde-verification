syms x1 x2 real
x = [x1;x2]; % column vector of symbolic variables
f = log(1 + 3*(x2 - (x1^3 - x1))^2 + (x1 - 4/3)^2);
% gradf = jacobian(f,x).'; % column gradf
% hessf = jacobian(gradf,x);
fh = matlabFunction(f,'vars',{x});
options = optimoptions('fminunc', ...
    'SpecifyObjectiveGradient', true, ...
    'HessianFcn', 'objective', ...
    'Algorithm','trust-region', ...
    'Display','final');
% [xfinal,fval,exitflag,output] = fminunc(fh,[-1,2],options);