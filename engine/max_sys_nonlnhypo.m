function [maxlist, fmin] = max_sys_nonlnhypo(u_min, u_max, deltat, deltax, xlist, gamma)

    dx = deltax;
    dt = deltat;
    mesh = size(xlist);
    m = mesh(2);     
    maxlist_1 = zeros(1,m);     %for rho
    maxlist_2 = zeros(1,m);     %for momemtum
    maxlist_3 = zeros(1,m);     %for energy
    maxlist = [maxlist_1; maxlist_2; maxlist_3];
    maxlist(1, 1) = 1; %u_max(1, 1);
    maxlist(1, m) = 0.7; %u_max(1, m);
    
    function result = fun(rho, mmt, engy, gamma) %flux function
        f1 = mmt;
        f2 = 1/2*(3 - gamma)*mmt^2/rho + (gamma - 1)*engy;
        f3 = gamma*mmt*engy/rho - 1/2*(gamma - 1)*mmt^2/(rho^2);
        result = [f1; f2; f3];
    end

    tgt = 0;
    
    for i = 1 : 3       % 3 variables for euler equations
        
        for j = 2 : m - 1  %construct min function

            rho_1 = strcat('rho_', int8(num2str(j-1)));
            rho_2 = strcat('rho_', int8(num2str(j)));
            rho_3 = strcat('rho_', int8(num2str(j+1)));
            r1 = str2sym(rho_1);
            r2 = str2sym(rho_2);        
            r3 = str2sym(rho_3);

            mmt_1 = strcat('mmt_', int8(num2str(j-1)));
            mmt_2 = strcat('mmt_', int8(num2str(j)));
            mmt_3 = strcat('mmt_', int8(num2str(j+1)));
            m1 = str2sym(mmt_1);
            m2 = str2sym(mmt_2);        
            m3 = str2sym(mmt_3);
            
            
            engy_1 = strcat('engy_', int8(num2str(j-1)));
            engy_2 = strcat('engy_', int8(num2str(j)));
            engy_3 = strcat('engy_', int8(num2str(j+1)));
            e1 = str2sym(engy_1);
            e2 = str2sym(engy_2);        
            e3 = str2sym(engy_3);
            
            k1 = [r1; m1; e1];
            k2 = [r2; m2; e2];
            k3 = [r3; m3; e3];
            u = [r1, r2, r3, m1, m2, m3, e1, e2, e3];   %r1 = rho at x = 1
            
            
%             if j == 2
%                 u = [k1, k2, k3];
%             else
%                 u = [u, k3];
%             end

            u_1 = (k3 + k2)/2 - dt/(2 * dx) * (fun(k3(1), k3(2), k3(3), gamma) - fun(k2(1), k2(2), k2(3), gamma));  %u_+1/2
            u_2 = (k2 + k1)/2 - dt/(2 * dx) * (fun(k2(1), k2(2), k2(3), gamma) - fun(k1(1), k1(2), k1(3), gamma));  %u_-1/2

            tgt = k2 - dt/dx * (fun(u_1(1), u_1(2), u_1(3), gamma) - fun(u_2(1), u_2(2), u_2(3), gamma));     %f(u_i+1/2)-f(u_i-1/2)    

%             tgt = tgt_new + tgt;

            tgt = -tgt;
                %solve 3 optimization problems at each stage
            ht_1 = matlabFunction(tgt(1),'vars',{u});
            ht_2 = matlabFunction(tgt(2),'vars',{u});
            ht_3 = matlabFunction(tgt(3),'vars',{u});    

            
            lb = [u_min(1, j - 1), u_min(1, j), u_min(1, j + 1), ...
                u_min(2, j - 1), u_min(2, j), u_min(2, j + 1), ...
                u_min(3, j - 1), u_min(3, j), u_min(3, j + 1)];
            ub = [u_max(1, j - 1), u_max(1, j), u_max(1, j + 1),...
                u_max(2, j - 1), u_max(2, j), u_max(2, j + 1),...
                u_max(3, j - 1), u_max(3, j), u_max(3, j + 1)];
            x0 = (lb + ub)/2;

            [list_1, fmin1] = fmincon(ht_1, x0,[],[],[],[], lb, ub);
            [list_2, fmin2] = fmincon(ht_2, x0,[],[],[],[], lb, ub);
            [list_3, fmin3] = fmincon(ht_3, x0,[],[],[],[], lb, ub);
            
            maxlist(1,j) = -fmin1;  %max value
            maxlist(2,j) = -fmin2;
            maxlist(3,j) = -fmin3;
            
%             disp(-fmin1);
%             disp(-fmin2);
%             disp(-fmin3);    
        end 
        
    end
    
    
    
%     for j = 2 : m - 1  %calculate rho for next step
%         
%         maxlist_1(:, 1) = [0.5; 0; 0.7];      %BC conditions
%         maxlist_1(:, m) = [1.5; 0; 1.5];  
%         
%         k1 = list_1(:, j - 1);
%         k2 = list_1(:, j);        
%         k3 = list_1(:, j + 1);
% 
%         
% 
%         u_1 = (k3 + k2)/2 - dt/(2 * dx) * (fun(k3(1), k3(2), k3(3), gamma) - fun(k2(1), k2(2), k2(3), gamma));  %u_+1/2
%         u_2 = (k2 + k1)/2 - dt/(2 * dx) * (fun(k2(1), k2(2), k2(3), gamma) - fun(k1(1), k1(2), k1(3), gamma));  %u_-1/2
% 
%         maxlist_1(:, j) = k2 - dt/dx * (fun(u_1(1), u_1(2), u_1(3), gamma) - fun(u_2(1), u_2(2), u_2(3), gamma));     %f(u_i+1/2)-f(u_i-1/2)    
%         
%     end        
%     
%     for j = 2 : m - 1  %calculate momentum for next step
%         
%         maxlist_2(:, 1) = [0.5; 0; 0.7];      %BC conditions
%         maxlist_2(:, m) = [1.5; 0; 1.5];       
%         
%         k1 = list_2(:, j - 1);
%         k2 = list_2(:, j);        
%         k3 = list_2(:, j + 1);
% 
%         
% 
%         u_1 = (k3 + k2)/2 - dt/(2 * dx) * (fun(k3(1), k3(2), k3(3), gamma) - fun(k2(1), k2(2), k2(3), gamma));  %u_+1/2
%         u_2 = (k2 + k1)/2 - dt/(2 * dx) * (fun(k2(1), k2(2), k2(3), gamma) - fun(k1(1), k1(2), k1(3), gamma));  %u_-1/2
% 
%         maxlist_2(:, j) = k2 - dt/dx * (fun(u_1(1), u_1(2), u_1(3), gamma) - fun(u_2(1), u_2(2), u_2(3), gamma));     %f(u_i+1/2)-f(u_i-1/2)    
%         
%     end    
%     
%     for j = 2 : m - 1  %calculate energy for next step
%         
%         maxlist_3(:, 1) = [0.5; 0; 0.7];      %BC conditions
%         maxlist_3(:, m) = [1.5; 0; 1.5]; 
%         
%         k1 = list_3(:, j - 1);
%         k2 = list_3(:, j);        
%         k3 = list_3(:, j + 1);
% 
%         
% 
%         u_1 = (k3 + k2)/2 - dt/(2 * dx) * (fun(k3(1), k3(2), k3(3), gamma) - fun(k2(1), k2(2), k2(3), gamma));  %u_+1/2
%         u_2 = (k2 + k1)/2 - dt/(2 * dx) * (fun(k2(1), k2(2), k2(3), gamma) - fun(k1(1), k1(2), k1(3), gamma));  %u_-1/2
% 
%         maxlist_3(:, j) = k2 - dt/dx * (fun(u_1(1), u_1(2), u_1(3), gamma) - fun(u_2(1), u_2(2), u_2(3), gamma));  %f(u_i+1/2)-f(u_i-1/2)    
%         
%     end   
%     
%     maxlist = [maxlist_1(1, :); maxlist_2(2, :); maxlist_3(3, :)];
end