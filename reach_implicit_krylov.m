function [sol1,sol2] = reach_implicit_krylov(deltat, deltax, init_min, init_max, lambda_min, lambda_max, time, xlist, tlist)
    % u(0) = 0 and u_x(L) = 0

    u_min = init_min;
    u_max = init_max;    
    
    sol1 = zeros(length(xlist),length(tlist));
    sol2 = zeros(length(xlist),length(tlist));    

    size2 = size(xlist);
    m = size2(2);
    sol1(:,1) = init_min;
    sol2(:,1) = init_max;
    
    
    r = deltat/(deltax*deltax);
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


    tic;    
    a = ones(m, 1) * (1 + 2*r);
    b = ones(m, 1) * (-r);
    c = ones(m, 1) * (-r);
    
    for i = 2 : time
        
        u_min = u_min *(1 + lambda_min * deltat);
        
        [temp, ~, relres]= gmres(A, u_min, 10, 0.001);
    
        temp_thomas_min = tridiag(a, b, c, u_min);
        
        for j = 1 : m
            
            temp(j) = temp(j) - relres * abs(u_min(j));
            temp(j) = min(temp_thomas_min(j), temp(j));
            
        end
        
        
        u_min = temp;
        sol1(:,i) = u_min; 

    end
        
    for i = 2 : time
        
        u_max = u_max *(1 + lambda_max * deltat);  
        
        [temp, ~, relres]= gmres(A, u_max, 10, 0.001);

        temp_thomas_max = tridiag(a, b, c, u_max);
        
        for j = 1 : m
            temp(j) = temp(j) + relres * abs(u_max(j));
            temp(j) = max(temp_thomas_max(j), temp(j));
        end
        
        u_max = temp;
        sol2(:,i) = u_max; 

    end        
    toc;
    
    figure;
    surf(xlist,tlist,sol1') 
    title('Minimum Numerical solution')
    xlabel('Distance x')
    ylabel('Time t')

    figure;
    surf(xlist,tlist,sol2') 
    title('Maximum Numerical solution')
    xlabel('Distance x')
    ylabel('Time t')
    
    %plot_3dbox(sol1', sol2', deltax, deltat, xlist, tlist);
end
    