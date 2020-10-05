function sol_print = solve_peaceman(alpha, deltat, deltax, deltay, lambda, init, time, xlist, ylist, tlist)
    % u(0) = 0 and u_x(L) = 0
    % solving 2d heat eq using ADI
    %add neumann and robin B.C. 
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
    Q3 = zeros(m, m);
    Q4 = zeros(m, m);    
    r1 = alpha * deltat/(2 * deltax * deltax);
    r2 = alpha * deltat/(2 * deltay * deltay);   
    
    f1 = zeros(m, 1);
    f2 = zeros(m, 1);
%----------------------step 1---- Q_1 u = Q_2 u + f--------------------

%Q1, x direction difference
    for i = 1 : m           %build matrix Q_1
        

            if rem(i, m1) == 1              %left boundary, u_x = 0
                

                Q1(i, i) =  1 + 2*r1;
                Q1(i, i + 1) = -2*r1;

                
            elseif rem(i, m1) == 0         %right boundary, u = 0.8
                

                Q1(i, i - 1) = -r1;
                Q1(i, i) =  1 + 2*r1;
                
                f1(i) = r1*0.8;
            
                
            else                                %interior points
                Q1(i, i - 1) = -r1;
                Q1(i, i) =  1 + 2*r1;
                Q1(i, i + 1) = -r1;

            end
                  
    end
    
    a1 = zeros(m, 1);
    b1 = zeros(m, 1);
    c1 = zeros(m, 1);
    
    for i = 1 : m
        
        a1(i) = Q1(i, i);
        
    end
    
    for i = 2 : m
        
        b1(i) = Q1(i, i - 1);    %subdiag        
        
    end
    
    
    for i = 1 : m - 1
        
        c1(i) = Q1(i, i + 1);    %upper diag           
        
    end   
    
%Q2, y direction difference 

     for i = 1 : m           %build matrix Q_2
        
            if i >= 1 && i <= m1         %lower boundary, u_y = 0.5

                Q2(i, i) =  1 - 2*r2;

                Q2(i, i + m1) = 2*r2;   
                
                f1(i) = - 0.5 * alpha * deltat/deltay;
                
            elseif i >= m - m1 + 1            %upper boundary, u_y = 1 - u
                
                Q2(i, i - m1) = r2;

                Q2(i, i) =  1 - 2*r2*(1 + deltay);

                f1(i) = 1 * alpha * deltat/deltay;
                               
            else                            %interior points
                Q2(i, i - m1) = r2;

                Q2(i, i) =  1 - 2*r2;
  
                Q2(i, i + m1) = r2;
            end
                  
    end  
    
    
    
    
%---------------------step 2--------- Q_3 u = Q_4 u + f--------------------   
    
    
    
    
    
    for i = 1 : m           %build matrix Q3
        
            if i >= 1 && i <= m1         %lower boundary u_y = 0.5

                Q3(i, i) =  1 + 2*r2;

                Q3(i, i + m1) = -r2;   
                
                f2(i) = - 0.5 * alpha * deltat/deltay;   
                
            elseif i >= m - m1 + 1            %upper boundary, u_y = 1 - u
                
                Q3(i, i - m1) = -r2;

                Q3(i, i) =  1 + 2*r2*(1 + deltay);

                f2(i) = 1 * alpha * deltat/deltay;     
                
            else                            %interior points
                Q3(i, i - m1) = -r2;

                Q3(i, i) =  1 + 2*r2;
  
                Q3(i, i + m1) = -r2;
            end
                  
    end    
    
    
     
    a2 = zeros(m, 1);
    b2 = zeros(m, 1);
    c2 = zeros(m, 1);
    
    for i = 1 : m
        
        a2(i) = Q3(i, i);
        
    end
    
    for i = 1 + m1 : m
        
        b2(i) = Q3(i, i - m1);    %subdiag        
        
    end
    
    
    for i = 1 : m - m1
        
        c2(i) = Q3(i, i + m1);    %upper diag           
        
    end 
    
    
    for i = 1 : m           %build matrix Q4
        

            if rem(i, m1) == 1              %left boundary, u_x = 0
                

                Q4(i, i) =  1 - 2*r1;
                Q4(i, i + 1) = 2*r1;

                
            elseif rem(i, m1) == 0         %right boundary, u = 0.8
                

                Q4(i, i - 1) = r1;
                Q4(i, i) =  1 - 2*r1;
            
                f2(i) = r1*0.8;   
                
            else                            %interior points
                
                Q4(i, i - 1) = r1;
                Q4(i, i) =  1 - 2*r1;
                Q4(i, i + 1) = r1;

            end
                  
    end
    
   
    
    for i = 2 : time
        
        %here solving matrxi at each step!!!
        u = Q2*u + lambda*deltat/2*u + f1;
        u_star = tridiag(a1, b1, c1, u);
        

        
        u_star = Q4*u_star + lambda*deltat/2*u_star + f2;
        u = block_tridiag(a2, b2, c2, m1, u_star);       

        sol(:,i) = u; 

    end
    
    
    
    
    
    sol9 = sol(:,91); %deltat = 0.1
    
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
    