
    dx = 1;
    dt = 1;
    
    tgt = 0;
    u_min = [0, 0.3, 0.3, 1, 0];
    u_max = [0, 1.5, 1.7, 1.3, 0];
    for j = 2 : 5 - 1
        
        str_1 = strcat('u_', int8(num2str(j-1)));
        str_2 = strcat('u_', int8(num2str(j)));
        str_3 = strcat('u_', int8(num2str(j+1)));
        
        k1 = str2sym(str_1);
        k2 = str2sym(str_2);        
        k3 = str2sym(str_3);

        if j == 2
            u = [k1, k2, k3];
        else
            u = [u, k3];
        end

        
%         u_1 = (k3 + k2)/2 - dt/(2 * dx) * (k3^2 - k2^2);  %u_+1/2
%         u_2 = (k2 + k1)/2 - dt/(2 * dx) * (k2^2 - k1^2);  %u_-1/2
        
        u_1 = (k3 + k2)/2 - dt/(2 * dx) * (fun(k3) - fun(k2));  %u_+1/2
        u_2 = (k2 + k1)/2 - dt/(2 * dx) * (fun(k2) - fun(k1));  %u_-1/2

        tgt_new = k2 - dt/dx * (fun(u_1) - fun(u_2));     %f(u_i+1/2)-f(u_i-1/2)    
%         tgt_new = k1^2 - k2 - k3;
        tgt = tgt_new + tgt;
    end 
    
    ht = matlabFunction(tgt,'vars',{u});
    x0 = (u_min + u_max)/2;
    lb = u_min;
    ub = u_max;

    [umin, fmin] = fmincon(ht, x0,[],[],[],[], lb, ub);
    
    function f = fun(x) %flux function
        f = x * x;
    end    
% end