function [minlist, fmin] = min_nonlnhypo(u_min, u_max, deltat, deltax, xlist)
%ingnore this

    dx = deltax;
    dt = deltat;
    mesh = size(xlist);
    m = mesh(2);     
    minlist = zeros(1,m);
    
    function f = fun(x) %flux function
        f = x * x;
    end

    tgt = 0;
    
    for j = 2 : m - 1
        
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
        
        u_1 = (k3 + k2)/2 - dt/(2 * dx) * (fun(k3) - fun(k2));  %u_+1/2
        u_2 = (k2 + k1)/2 - dt/(2 * dx) * (fun(k2) - fun(k1));  %u_-1/2

        tgt_new = k2 - dt/dx * (fun(u_1) - fun(u_2));     %f(u_i+1/2)-f(u_i-1/2)    
        
        tgt = tgt_new + tgt;
    end        
        
    ht = matlabFunction(tgt,'vars',{u});
    x0 = (u_min + u_max)/2;
    lb = u_min;
    ub = u_max;

    [list, fmin] = fmincon(ht, x0,[],[],[],[], lb, ub);
    
    
    for j = 2 : m - 1   %calculate value for next step
        
        minlist(1) = 0;
        minlist(m) = 0;
        
        k1 = list(j - 1);
        k2 = list(j);        
        k3 = list(j + 1);

        

        u_1 = (k3 + k2)/2 - dt/(2 * dx) * (fun(k3) - fun(k2));  %u_+1/2
        u_2 = (k2 + k1)/2 - dt/(2 * dx) * (fun(k2) - fun(k1));  %u_-1/2

        minlist(j) = k2 - dt/dx * (fun(u_1) - fun(u_2));     %f(u_i+1/2)-f(u_i-1/2)    
        
    end        
    
    
    
    
    
% %             u_min_new(j) = min_nonlnhypo(0, u_min(j), u_min(j + 1), 0, u_max(j), u_max(j + 1), deltat, deltax);
% %             u_max_new(j) = max_nonlnhypo(0, u_min(j), u_min(j + 1), 0, u_max(j), u_max(j + 1), deltat, deltax);          
%             
% %         elseif j == m   %right boundary
% 
% %             u_min_new(j) = min_nonlnhypo(u_min(j - 1), u_min(j), 0, u_max(j - 1), u_max(j), 0, deltat, deltax);
% %             u_max_new(j) = max_nonlnhypo(u_min(j - 1), u_min(j), 0, u_max(j - 1), u_max(j), 0, deltat, deltax);
%             
% 
%         
%         else    %in the middle
%             x_j
%             tgt = tgt .* @loc;            
% %             u_min_new(j) = min_nonlnhypo(u_min(j - 1), u_min(j), u_min(j + 1), u_max(j - 1), u_max(j), u_max(j + 1), deltat, deltax);
% %             u_max_new(j) = max_nonlnhypo(u_min(j - 1), u_min(j), u_min(j + 1), u_max(j - 1), u_max(j), u_max(j + 1), deltat, deltax);
% %         end        
%     
    


% x0 = [(x_min + x_max)/2, (y_min + y_max)/2, (z_min + z_max)/2];
% lb = [x_min, y_min, z_min];
% ub = [x_max, y_max, z_max];

% ht = matlabFunction(tgt);
% x0 = (u_min + u_max)/2;
% lb = u_min;
% ub = u_max;
% 
% [~, fmin] = fmincon(ht, x0,[],[],[],[], lb, ub);
% disp('min is');
% disp(fmin);


end