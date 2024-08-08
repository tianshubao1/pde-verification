function [maxlist, fmax] = max_nonlnhypo_local(u_min, u_max, deltat, deltax, xlist, h_min, h_max)

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
        
        u = [k1, k2, k3, h1, h2, h3]; 

        u_1 = (k1 + k2)/2 - dt/(2 * dx) * (fun(k2) - fun(k1)) + dt/2* (h1 + h2)/2;  %u_-1/2
        u_2 = (k2 + k3)/2 - dt/(2 * dx) * (fun(k3) - fun(k2)) + dt/2* (h2 + h3)/2;  %u_+1/2
        tgt = k2 - dt/dx * (fun(u_2) - fun(u_1))+ dt*h2;              
        tgt = -tgt;    
        
        
        ht = matlabFunction(tgt,'vars',{u});
        lb = [u_min(j - 1), u_min(j), u_min(j + 1), h_min(j - 1), h_min(j), h_min(j + 1)];
        ub = [u_max(j - 1), u_max(j), u_max(j + 1), h_max(j - 1), h_max(j), h_max(j + 1)];
        x0 = (lb + ub)/2;
        [~, fmin] = fmincon(ht, x0,[],[],[],[], lb, ub);
        fmax = -fmin;
        maxlist(j) = fmax;
    end    
    

end