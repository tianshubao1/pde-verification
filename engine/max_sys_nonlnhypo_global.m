function [maxlist, fmin] = max_sys_nonlnhypo_global(u_min, u_max, deltat, deltax, xlist, gamma, h_min, h_max)

    dx = deltax;
    dt = deltat;
    mesh = size(xlist);
    m = mesh(2);   
    
    maxlist_1 = zeros(1,m);     %for rho
    maxlist_2 = zeros(1,m);     %for momemtum
    maxlist_3 = zeros(1,m);     %for energy
    
    maxlist = [maxlist_1; maxlist_2; maxlist_3];

    
    tgt = 0;
    
    
    function result = fun(rho, mmt, engy, gamma) %flux function
        f1 = mmt;
        f2 = 1/2*(3 - gamma)*mmt^2/(rho + 0.1) + (gamma - 1)*engy;
        f3 = gamma*mmt*engy/(rho + 0.1) - 1/2*(gamma - 1)*mmt^2/(rho*rho + 0.1);
        result = [f1; f2; f3];
    end

            
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
            
            k1 = [r1; m1; e1];      %combine them to a 3*1 vector
            k2 = [r2; m2; e2];
            k3 = [r3; m3; e3];
            
            hin_1 = strcat('input_', int8(num2str(j-1)));       % only consider 1 h input
            hin_2 = strcat('input_', int8(num2str(j)));
            hin_3 = strcat('input_', int8(num2str(j+1)));
            h1 = str2sym(hin_1);
            h2 = str2sym(hin_2);        
            h3 = str2sym(hin_3);
            
            
            if j == 2
                u1 = [k1, k2, k3];      %this is a matrix, 3*3      
                u2 = [h1, h2, h3];      %this is just for h, 1*3    
            else
                u1 = [u1, k3];  % 3* m     [r1,r2,r3,..; m1,m2,m3,...;e1,e2,e3,...]
                u2 = [u2, h3];  % 1* m     [h1,h2,h3,...]
            end
            
            
            u_1 = (k2 + k1)/2 - dt/(2 * dx) * (fun(k2(1), k2(2), k2(3), gamma) - fun(k1(1), k1(2), k1(3), gamma)) + dt/2* ([h1 + h2; 0; 0])/2;  %u_-1/2
            u_2 = (k3 + k2)/2 - dt/(2 * dx) * (fun(k3(1), k3(2), k3(3), gamma) - fun(k2(1), k2(2), k2(3), gamma)) + dt/2* ([h2 + h3; 0; 0])/2;  %u_+1/2     3*1 vector

            tgt_new = k2 - dt/dx * (fun(u_2(1), u_2(2), u_2(3), gamma) - fun(u_1(1), u_1(2), u_1(3), gamma))+ dt*[h2; 0; 0];     %f(u_i+1/2)-f(u_i-1/2)    3*1 vector

            tgt = tgt_new.^2 + tgt;
        end
        
        tgt = -tgt;
        u = [u1; u2];   %[4, m]
        ht_1 = matlabFunction(tgt(1),'vars',{u});
        ht_2 = matlabFunction(tgt(2),'vars',{u});
        ht_3 = matlabFunction(tgt(3),'vars',{u}); 

            

        lb = [u_min; h_min'];   %[4, m]
        ub = [u_max; h_max'];   %[4, m]
        x0 = (lb + ub)/2;
        

        [list1, ~] = fmincon(ht_1, x0,[],[],[],[], lb, ub);
        [list2, ~] = fmincon(ht_2, x0,[],[],[],[], lb, ub);
        [list3, ~] = fmincon(ht_3, x0,[],[],[],[], lb, ub);

%         if fmin1 <= 0
%             fmin1 = 0.01;
%         end
%             
%         if fmin2 <= 0
%             fmin2 = 0.01;
%         end
%             
%         if fmin3 <= 0
%             fmin3 = 0.01;
%         end
        
        list = [list1; list2; list3];    %[4*3, m]
            
        
        
    for i = 1 : 3    % for 3 variables, pick minimum for them seperately
        for j = 2 : m - 1   %calculate value for next step

            maxlist(1) = 0;     %[3, m]
            maxlist(m) = 0;
            
            currlist = list(i*4-3 : i*4 -1, :);   %[3, m]
            
            k1 = currlist(:,j - 1);   %[3, 1]
            k2 = currlist(:,j);        %[3, 1]
            k3 = currlist(:,j + 1);       %[3, 1]
            
            h1 = list(i*4, j - 1);
            h2 = list(i*4, j);
            h3 = list(i*4, j + 1);
 
            u_1 = (k2 + k1)/2 - dt/(2 * dx) * (fun(k2(1), k2(2), k2(3), gamma) - fun(k1(1), k1(2), k1(3), gamma)) + dt/2* ([h1 + h2; 0; 0])/2;  %u_-1/2
            u_2 = (k3 + k2)/2 - dt/(2 * dx) * (fun(k3(1), k3(2), k3(3), gamma) - fun(k2(1), k2(2), k2(3), gamma)) + dt/2* ([h2 + h3; 0; 0])/2;  %u_+1/2     3*1 vector
            result = k2 - dt/dx * (fun(u_2(1), u_2(2), u_2(3), gamma) - fun(u_1(1), u_1(2), u_1(3), gamma))+ dt*[h2; 0; 0];  %[3, 1]
            
            if(result(i) <= 0)
                result(i) = 0.02;
            end
            
            maxlist(i,j) = result(i);
%             maxlist(:,j) = result;
            

        
        end         
        
    end   

        

    
    
end