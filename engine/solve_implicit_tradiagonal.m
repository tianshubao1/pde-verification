function sol = solve_implicit_tradiagonal(deltat, deltax, init, time, xlist, tlist)
    % u(0) = 0 and u_x(L) = 0

    u = init;
    sol = zeros(length(xlist),length(tlist));
    size(sol);
    size2 = size(xlist);
    m = size2(2);
    sol(:,1) = init;
    r = deltat/(deltax*deltax);
    A = zeros(m, m);
    
    for i = 1 : m
          
            if i == 1
                A(i, i) =  1 + 2*r;
                A(i, i + 1) = -r;
                
            elseif i == m
                A(i, i - 1) = -r;
                A(i, i) =  1 + r;
                
            else
                A(i, i - 1) = -r;
                A(i, i) =  1 + 2*r;
                A(i, i + 1) = -r;
            end
                  
    end
    
        
    for i = 2 : time
        
        u_new = gmres(A, u, 10, 0.001);
        
        u_max *(1 + lambda_max * deltat);
        y = tridiag(1 + 2*r, -r, -r, u );
        
%         u = A\u;
        sol(:,i) = u_new; 

    end
        
        
    
    figure;
    surf(xlist,tlist,sol') 
    title('Numerical solution')
    xlabel('Distance x')
    ylabel('Time t')


end
    