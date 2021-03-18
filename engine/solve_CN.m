function sol_print = solve_CN(deltat, deltax, deltay, lambda, alpha, init, time, xlist, ylist, tlist, operator)
    % u(0) = 0 and u_x(L) = 0

    u = init;
       
    
    size1 = size(xlist);
    m1 = size1(2);
    size2 = size(ylist);
    m2 = size2(2);
    m = m1 * m2;
    
    sol = zeros(m,length(tlist));
    sol(:,1) = init;
    Q1 = zeros(m, m);
    Q2 = zeros(m, m);
    r1 = alpha * deltat/(2 * deltax * deltax);
    r2 = alpha * deltat/(2 * deltay * deltay);

    f = zeros(m, 1);
% update boundary condition    
%-----------------Q1 * u^{n + 1} = Q2 * u^{n}------------------------------ 
%-----------------build matrix Q1, blocked tridiagonal---------------------

    for i = 1 : m           
        
            Q1(i, i) =  1 + 2*r1 + 2*r2 - 1/2 * deltat * lambda;
            
            
            if rem(i, m1) == 1              %left boundary, u_x = 0
                  
                Q1(i, i + 1) = -r1;
                Q1(i, i) =  1 + r1 + 2*r2 - 1/2 * deltat * lambda;
                
            elseif rem(i, m1) == 0         %right boundary, u = 0.8
              
                Q1(i, i - 1) = -r1;
  
                           
            else                                %interior points
                Q1(i, i - 1) = -r1;
                
                Q1(i, i + 1) = -r1;                 
               
            end
            
            if i >= 1 && i <= m1         %lower boundary, u_y = 0.5   
                
                Q1(i, i + m1) = -2*r2;  
                Q1(i, i) =  1 + 2*r1 + 2*r2 - 1/2 * deltat * lambda;
                
                
            elseif i >= m - m1 + 1        %upper boundary, u_y = 1 - u 
                
                Q1(i, i - m1) = -2*r2;
                Q1(i, i) =  1 + 2*r1 + 2*r2*(1 + deltax) - 1/2 * deltat * lambda;
                         
            else
                
                Q1(i, i + m1) = -r2;  
                Q1(i, i - m1) = -r2;
            end
                  
    end    
    
    %disp(Q1);
    

%-----------------build matrix Q2, blocked tridiagonal-------------------

    for i = 1 : m           
        
            Q2(i, i) =  1 - 2*r1 - 2*r2 + 1/2 * deltat * lambda;
            
            
            if rem(i, m1) == 1              %left boundary, u_x = 0
                
 
                Q2(i, i + 1) = r1;
                Q2(i, i) =  1 - r1 - 2*r2 + 1/2 * deltat * lambda;
                
            elseif rem(i, m1) == 0         %right boundary, u = 0.8
                
                Q2(i, i - 1) = r1;

                f(i) = 0.8 * r1 + 0.8 * r1;
                
            else                                %interior points
                Q2(i, i - 1) = r1;
                
                Q2(i, i + 1) = r1;

                
            end
            
            if i >= 1 && i <= m1         %lower boundary, u_y = 0.5  
                
                Q2(i, i + m1) = r2;  
                Q2(i, i) =  1 - r1 - 2*r2 + 1/2 * deltat * lambda;
%                 Q2(i, i + m1) = 2*r2;  
%                 Q2(i, i) =  1 - 2*r1 - 2*r2 + 1/2 * deltat * lambda;                
                
                f(i) = - 0.5 * 2 * alpha * deltat/deltay;
                
            elseif i >= m - m1 + 1       %upper boundary, u_y = 1 - u 
                
                Q2(i, i - m1) = 2*r2;
                Q2(i, i) =  1 - 2*r1 - 2*r2*(1 + deltay) + 1/2 * deltat * lambda;
                f(i) =   1 * 2 * alpha* deltat/deltay;
                
            else
                Q2(i, i + m1) = r2;  
                Q2(i, i - m1) = r2;
            end
                  
    end    
    
    %disp(Q2);
     
   
    %bloat this solution later
    for i = 2 : time
        
        u_curr = Q2 * u + f;   
        [u, ~, relres] = gmres(Q1, u_curr, 10, 0.001);
        
        for j = 1 : m      % bloat
            
            
            if operator == "min"
                u(j) = u(j) - relres * abs(u(j));   
            elseif operator == "max"
                u(j) = u(j) + relres * abs(u(j));   
            end
            
            
        end
        
        sol(:,i) = u; 

    end
    
    sol9 = sol(:,91); %deltat = 0.1, t = 9s
    
    sol_print = zeros(m2, m1);      %m2 index for y, number of rows, m1 index for x, number of column
    
    for i = 1 : m2
        for j = 1 : m1
            
            sol_print(i, j) = sol9(m1 * (i - 1) + j, 1);
            
        end
    end

    
    figure;
    surf(xlist,ylist,sol_print); 
    title('Numerical solution')
    xlabel('Distance x')
    ylabel('Distance y')


end
    