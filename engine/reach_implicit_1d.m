function [sol_min, sol_max] = reach_implicit_1d(deltat, deltax, init_min, init_max, lambda_min, lambda_max, time, xlist, tlist)
    % u(0) = 0 and u_x(L) = 0
    
    tic
    u_min = init_min;
    u_max = init_max;    
    
    sol_min = zeros(length(xlist),length(tlist));
    sol_max = zeros(length(xlist),length(tlist));    

    size2 = size(xlist);
    m = size2(2);
    sol_min(:,1) = init_min;
    sol_max(:,1) = init_max;
        
    r = deltat/(deltax*deltax);     %alpha = 1
    A = zeros(m, m);
    
    for i = 1 : m
          
            if i == 1                   %dirichlet bd condition
                A(i, i) =  1 + 2*r;
                A(i, i + 1) = -r;
                
            elseif i == m               %neumann bd condition
                A(i, i - 1) = -2*r;
                A(i, i) =  1 + 2*r;
                
            else
                A(i, i - 1) = -r;
                A(i, i) =  1 + 2*r;
                A(i, i + 1) = -r;
            end
                  
    end
    
%-------------------pick the maximum of thomax algo and gminres    


   
    a = ones(m, 1) * (1 + 2*r);
    b = ones(m, 1) * (-r);
    b(m, 1) = -2*r;
    c = ones(m, 1) * (-r);
    
    for i = 2 : time    % max value
        
        u_min = u_min *(1 + lambda_min * deltat);
        
        [temp, ~, relres]= gmres(A, u_min, 10, 0.001);
    
        %temp_thomas_min = tridiag(a, b, c, u_min);
        
        for j = 1 : m
            
            temp(j) = temp(j) - relres * abs(u_min(j)); %   bloating
            
        end
        
        
        u_min = temp;
        sol_min(:,i) = u_min; 

    end
        
    for i = 2 : time     % max value
        
        u_max = u_max *(1 + lambda_max * deltat);  
        
        [temp, ~, relres]= gmres(A, u_max, 10, 0.001);

        %temp_thomas_max = tridiag(a, b, c, u_max);
        
        for j = 1 : m
            temp(j) = temp(j) + relres * abs(u_max(j));     %   bloating
        end
        
        u_max = temp;
        sol_max(:,i) = u_max; 

    end        
    toc;
    
    plot_2dbox(sol_min(:,4/deltat + 1), sol_max(:,4/deltat + 1), deltax, xlist);   
    axis([-1 11 0 inf])
%     
%     figure;
%     surf(xlist,tlist,sol_min') 
%     title('Minimum Numerical solution')
%     xlabel('Distance x')
%     ylabel('Time t')
% 
%     figure;
%     surf(xlist,tlist,sol_max') 
%     title('Maximum Numerical solution')
%     xlabel('Distance x')
%     ylabel('Time t')
    
%     plot_3dbox(sol_min', sol_max', deltax, deltat, xlist, tlist);
end
    