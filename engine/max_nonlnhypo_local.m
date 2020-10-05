function [maxlist, fmax] = max_nonlnhypo_local(u_min, u_max, deltat, deltax, xlist)

    dx = deltax;
    dt = deltat;
    mesh = size(xlist);
    m = mesh(2);     
    maxlist = zeros(1,m);   
    
    function f = fun(x) %flux function
        f = x * x;
    end 

    tgt = 0;
    
    maxlist(1) = 0;
    maxlist(m) = 0;      
    for j = 2 : m - 1
        
        str_1 = strcat('u_', int8(num2str(j-1)));
        str_2 = strcat('u_', int8(num2str(j)));
        str_3 = strcat('u_', int8(num2str(j+1)));
        
        k1 = str2sym(str_1);
        k2 = str2sym(str_2);        
        k3 = str2sym(str_3);

        u = [k1, k2, k3];

        
        u_1 = (k3 + k2)/2 - dt/(2 * dx) * (fun(k3) - fun(k2));  %u_+1/2
        u_2 = (k2 + k1)/2 - dt/(2 * dx) * (fun(k2) - fun(k1));  %u_-1/2

        tgt = k2 - dt/dx * (fun(u_1) - fun(u_2));     %f(u_i+1/2)-f(u_i-1/2)
        
        tgt = -tgt;    
        ht = matlabFunction(tgt,'vars',{u});
        x0 = [u_min(j - 1), u_min(j), u_min(j + 1)];
        lb = [u_min(j - 1), u_min(j), u_min(j + 1)];
        ub = [u_max(j - 1), u_max(j), u_max(j + 1)];  

        [list, fmin] = fmincon(ht, x0,[],[],[],[], lb, ub);
        fmax = -fmin;
        maxlist(j) = fmax;
    end    
    

    
%     for j = 2 : m - 1   %calculate value for next step
%         
%         maxlist(1) = 0;
%         maxlist(m) = 0;
%         
%         k1 = list(j - 1);
%         k2 = list(j);        
%         k3 = list(j + 1);
% 
%         u_1 = (k3 + k2)/2 - dt/(2 * dx) * (fun(k3) - fun(k2));  %u_+1/2
%         u_2 = (k2 + k1)/2 - dt/(2 * dx) * (fun(k2) - fun(k1));  %u_-1/2
% 
%         maxlist(j) = k2 - dt/dx * (fun(u_1) - fun(u_2));     %f(u_i+1/2)-f(u_i-1/2)    
%         
%     end        
    
% dx = deltax;
% dt = deltat;
% 
%     function tgt = target(x)
%         
%         function f = fun(x) %flux function
%             f = x * x;
%         end
%         
%         u_1 = (x(3) + x(2))/2 - dt/(2 * dx) * (fun(x(3)) - fun(x(2)));  %u_+1/2
%         u_2 = (x(2) + x(1))/2 - dt/(2 * dx) * (fun(x(2)) - fun(x(1)));  %u_-1/2
% 
%         tgt = x(2) - dt/dx * (fun(u_1) - fun(u_2));     %f(u_i+1/2)-f(u_i-1/2)
%         tgt = -tgt;
%     end
% 
% 
% x0 = [(x_min + x_max)/2, (y_min + y_max)/2, (z_min + z_max)/2];
% lb = [x_min, y_min, z_min];
% ub = [x_max, y_max, z_max];
% 
% [~, fmin] = fmincon(@target,x0,[],[],[],[],lb,ub);
% fmax = -fmin;


end