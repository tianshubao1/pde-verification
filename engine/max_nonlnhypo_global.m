function [maxlist, fmax] = max_nonlnhypo_global(u_min, u_max, deltat, deltax, xlist, h_min, h_max)

    dx = deltax;
    dt = deltat;
    mesh = size(xlist);
    m = mesh(2);     
    maxlist = zeros(1,m);   
    
    function f = fun(x) %flux function
        f = 1/2 * x * x;
    end 

    tgt = 0;
    
    maxlist(1) = 0; %boundary is 0
    maxlist(m) = 0;      
    for j = 2 : m - 1
        
        str_1 = strcat('u_', int8(num2str(j-1)));
        str_2 = strcat('u_', int8(num2str(j)));
        str_3 = strcat('u_', int8(num2str(j+1)));
        
        k1 = str2sym(str_1);        %u_1
        k2 = str2sym(str_2);        %u_2
        k3 = str2sym(str_3);        %u_3

        str_4 = strcat('h_', int8(num2str(j-1)));
        str_5 = strcat('h_', int8(num2str(j)));
        str_6 = strcat('h_', int8(num2str(j+1)));
        
        h1 = str2sym(str_4);        %h_1
        h2 = str2sym(str_5);        %h_2
        h3 = str2sym(str_6);        %h_3
        
        if j == 2
            u1 = [k1, k2, k3];
            u2 = [h1, h2, h3]; 
        else
            u1 = [u1, k3];
            u2 = [u2, h3];
        end
        
        u_1 = (k1 + k2)/2 - dt/(2 * dx) * (fun(k2) - fun(k1)) + dt/2* (h1 + h2)/2;  %u_-1/2
        u_2 = (k2 + k3)/2 - dt/(2 * dx) * (fun(k3) - fun(k2)) + dt/2* (h2 + h3)/2;  %u_+1/2
        tgt_new = k2 - dt/dx * (fun(u_2) - fun(u_1))+ dt*h2;  
        
        tgt = tgt_new^2 + tgt;
        
    end
    
    
    tgt = -tgt;    
    u = [u1, u2];       %[1, 42]    
    ht = matlabFunction(tgt,'vars',{u});
    lb = [u_min, h_min'];
    ub = [u_max, h_max'];
    x0 = (lb + ub)/2;
    [list, fmin] = fmincon(ht, x0,[],[],[],[], lb, ub);
    fmax = -fmin;
        
    for j = 2 : m - 1   %calculate value for next step
        
        maxlist(1) = 0;
        maxlist(m) = 0;
        
        k1 = list(j - 1);
        k2 = list(j);        
        k3 = list(j + 1);

        u_1 = (k3 + k2)/2 - dt/(2 * dx) * (fun(k3) - fun(k2));  %u_+1/2
        u_2 = (k2 + k1)/2 - dt/(2 * dx) * (fun(k2) - fun(k1));  %u_-1/2

        maxlist(j) = k2 - dt/dx * (fun(u_1) - fun(u_2));     %f(u_i+1/2)-f(u_i-1/2)    
        
    end       

end