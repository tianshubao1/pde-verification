function sol = solve_peaceman(deltat, deltax, deltay, lambda, init, time, xlist, ylist, tlist)
    % u(0) = 0 and u_x(L) = 0

    u = init;
    sol = zeros(length(xlist),length(tlist));
    size(sol);
    
    
    size1 = size(xlist);
    m1 = size1(2);
    size2 = size(ylist);
    m2 = size2(2);
    m = m1 * m2;
    
    sol(:,1) = init;
    Q1 = zeros(m, m);
    Q2 = zeros(m, m);
    Q3 = zeros(m, m);
    Q4 = zeros(m, m);
    r1 = deltat/(2 * deltax * deltax);
    r2 = deltat/(2 * deltay * deltay);

    
    
    
%--------------------------------step 1------------------------------------


    for i = 1 : m           %build matrix Q
        

            if rem(i, m2) == 1              %left boundary
                

                Q1(i, i) =  1 + 2*r1 + 2*r2;
                Q1(i, i + 1) = -r1;

                
            elseif rem(i, m2) == m2 - 1         %right boundary
                

                Q1(i, i - 1) = -r1;
                Q1(i, i) =  1 + 2*r1;
            
                
            else                                %interior points
                Q1(i, i - 1) = -r1;
                Q1(i, i) =  1 + 2*r1;
                Q1(i, i + 1) = -r1;

            end
                  
    end
    
 
     for i = 1 : m           %build matrix Q
        
            if i > 1 && i < m         %upper boundary

                Q2(i, i) =  1 + 2*r2;

                Q2(i, i + m1) = -r2;   
                
            elseif i > m - m1 + 1            %lower boundary
                
                Q2(i, i - m1) = -r2;

                Q2(i, i) =  1 + 2*r2;

                               
            else                            %interior points
                Q2(i, i - m1) = -r2;

                Q2(i, i) =  1 + 2*r2;
  
                Q2(i, i + m1) = -r2;
            end
                  
    end  
    
    for i = 2 : time
        
        u = u + deltat * lambda * u;   
        u = gmres(Q1, u, 10, 0.001);
        sol(:,i) = u; 

    end
        
        
    
    
    
    
    
%--------------------------------step 2------------------------------------    
    
    
    
    
    
      for i = 1 : m           %build matrix Q
        
            if i > 1 && i < m         %upper boundary

                Q3(i, i) =  1 + 2*r2;

                Q3(i, i + m1) = -r2;   
                
            elseif i > m - m1 + 1            %lower boundary
                
                Q3(i, i - m1) = -r2;

                Q3(i, i) =  1 + 2*r2;

                               
            else                            %interior points
                Q3(i, i - m1) = -r2;

                Q3(i, i) =  1 + 2*r2;
  
                Q3(i, i + m1) = -r2;
            end
                  
    end    
    
    
    
    for i = 1 : m           %build matrix Q
        

            if rem(i, m2) == 1              %left boundary
                

                Q4(i, i) =  1 + 2*r1 + 2*r2;
                Q4(i, i + 1) = -r1;

                
            elseif rem(i, m2) == m2 - 1         %right boundary
                

                Q4(i, i - 1) = -r1;
                Q4(i, i) =  1 + 2*r1;
            
                
            else                            %interior points
                
                Q4(i, i - 1) = -r1;
                Q4(i, i) =  1 + 2*r1;
                Q4(i, i + 1) = -r1;

            end
                  
    end
    
 
 
    
 
   
    
    for i = 2 : time
        
        u = u + deltat * lambda * u;   
        u = gmres(Q1, u, 10, 0.001);
        sol(:,i) = u; 

    end
    
    
    figure;
    surf(xlist,tlist,sol') 
    title('Numerical solution')
    xlabel('Distance x')
    ylabel('Time t')


end
    