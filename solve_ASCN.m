function sol = solve_ASCN(deltat, deltax, alpha, init, time, xlist, tlist, lambda)
    
    % solve a 1d heat equation using Alternating Segment Crank–Nicolson scheme
    % u(0) = 0 and u_x(L) = 0
    % consider 3 parallel threads, the size of 1st is m, the 2nd and 3rd are 2*m
    
    % todos add BCs and lambda
    u = init;
       
    size1 = size(xlist);
    m = size1(2);
    l = m / 5;
    
    sol = zeros(m,length(tlist));
    sol(:, 1) = init;
    
    
    Gm1 = zeros(l, l);
    Gm2 = zeros(l, l);
    G2m = zeros(2*l, 2*l);
    G2m1 = zeros(2*l, 2*l);
    G2m2 = zeros(2*l, 2*l);
    
    r = alpha * deltat/(deltax * deltax);

    %f = zeros(m, 1);
    



%-----------------construct G_l^1---------------------

    for i = 1 : l
        
        if i == 1
            
            Gm1(i, i) = 2;
            Gm1(i, i + 1) = -1;
              
        elseif i == l
            
            Gm1(i, i) = 1;
            Gm1(i, i - 1) = -1;
            
        else            

            Gm1(i, i) = 2;
            Gm1(i, i - 1) = -1;        
            Gm1(i, i + 1) = -1;  
            
        end
        
    end
    
    Gm1 = 0.5*Gm1;

    
%-----------------construct G_l^2---------------------    
    
    for i = 1 : l
        
        if i == 1
            
            Gm2(i, i) = 1;
            Gm2(i, i + 1) = -1;
              
        elseif i == l
            
            Gm2(i, i) = 1;                  %neumann boundary condition on right bondary 
            Gm2(i, i - 1) = -1;
            
        else            

            Gm2(i, i) = 2;
            Gm2(i, i - 1) = -1;        
            Gm2(i, i + 1) = -1;  
            
        end
        
    end
    Gm2 = 0.5*Gm2;

%-----------------construct G_2l---------------------    
    
    for i = 1 : 2*l
        
        if i == 1
            
            G2m(i, i) = 1;
            G2m(i, i + 1) = -1;
              
        elseif i == l
            
            G2m(i, i) = 3;
            G2m(i, i - 1) = -1;      
            G2m(i, i + 1) = -2; 
            
        elseif i == l + 1            

            G2m(i, i) = 3;
            G2m(i, i - 1) = -2;        
            G2m(i, i + 1) = -1; 

        elseif i == 2*l            

            G2m(i, i) = 1;
            G2m(i, i - 1) = -1;  
            
        else            

            G2m(i, i) = 2;
            G2m(i, i - 1) = -1;        
            G2m(i, i + 1) = -1; 
            
        end
        
    end
    G2m = 0.5*G2m;
    
%-----------------construct G_2l^1---------------------    
    
    for i = 1 : 2*l
        
        if i == 1
            
            G2m1(i, i) = 1;
            G2m1(i, i + 1) = -1;
              
        elseif i == l
            
            G2m1(i, i) = 3;
            G2m1(i, i - 1) = -1;      
            G2m1(i, i + 1) = -2; 
            
        elseif i == l + 1            

            G2m1(i, i) = 3;
            G2m1(i, i - 1) = -2;        
            G2m1(i, i + 1) = -1; 

        elseif i == 2*l            

            G2m1(i, i) = 1;                 %neumann boundary condition on right bondary                  
            G2m1(i, i - 1) = -1;  
            
        else            

            G2m1(i, i) = 2;
            G2m1(i, i - 1) = -1;        
            G2m1(i, i + 1) = -1; 
            
        end
        
    end
    G2m1 = 0.5*G2m1;
    
    
%-----------------construct G_2l^2---------------------    
    
    for i = 1 : 2*l
        
        if i == 1
            
            G2m2(i, i) = 2;
            G2m2(i, i + 1) = -1;
              
        elseif i == l
            
            G2m2(i, i) = 3;
            G2m2(i, i - 1) = -1;      
            G2m2(i, i + 1) = -2; 
            
        elseif i == l + 1            

            G2m2(i, i) = 3;
            G2m2(i, i - 1) = -2;        
            G2m2(i, i + 1) = -1; 

        elseif i == 2*l           

            G2m2(i, i) = 1;
            G2m2(i, i - 1) = -1;  
            
        else            

            G2m2(i, i) = 2;
            G2m2(i, i - 1) = -1;        
            G2m2(i, i + 1) = -1; 
            
        end
        
    end
    G2m2 = 0.5*G2m2;
   
%-------------------------apply ASC-N scheme--------------------------     
%-----------------(I + rG2) * u^{n + 1} = (I - rG1) * u^{n} for n is odd------------------------------     
%-----------------(I + rG1) * u^{n + 1} = (I - rG2) * u^{n} for n is even------------------------------    

    
    pool = parpool(3);
    
    tic
    
    for i = 2 : time
        
        if mod(i,2) == 0    %we want to update t = i * delta t
%-------------------------1st thread-------------------------- 
            
            temp{1} = u(1:l, 1);
            temp{2} = u(l+1:3*l, 1);
            temp{3} = u(3*l+1:5*l, 1);
            
            parfor j = 1:3
                
                if j == 1
                    
                    temp{j} = (eye(l) - r*Gm1) * temp{j} + deltat*lambda*temp{j};

                end
                
                
%-------------------------2nd thread-------------------------- 

                if j == 2
                    
                    temp{j} = (eye(2*l) - r*G2m) * temp{j}+ deltat*lambda*temp{j};

                end
                
%-------------------------3rd thread--------------------------
                if j == 3
                    
                    temp{j} = (eye(2*l) - r*G2m1) * temp{j}+ deltat*lambda*temp{j};


                end
            end
            
            u(1:l, 1) = temp{1};            
            u(l+1:3*l, 1) = temp{2};   
            u(3*l+1:5*l, 1) = temp{3};            

            
            update{1} = u(1:2*l, 1);            
            update{2} = u(2*l+1:4*l, 1);   
            update{3} = u(4*l+1:5*l, 1);  
            
  
            parfor j = 1:3
                
                if j == 1
                    %update{j} = (eye(2*l) + r*G2m2) \ update{j};         
                    [update{j}, ~, ~]= gmres((eye(2*l) + r*G2m2), update{j}, 10, 0.001);

                end
                
                
%-------------------------2nd thread-------------------------- 

                if j == 2
                    
                    %update{j} = (eye(2*l) + r*G2m) \ update{j};
                    [update{j}, ~, ~]= gmres((eye(2*l) + r*G2m), update{j}, 10, 0.001);

                end
                
%-------------------------3rd thread--------------------------
                if j == 3
                    
                    %update{j} = (eye(l) + r*Gm2) \ update{j};
                    [update{j}, ~, ~]= gmres((eye(l) + r*Gm2), update{j}, 10, 0.001);

                end
            end   
            
            u(1:2*l, 1) = update{1};  
            u(2*l+1:4*l, 1) = update{2};    
            u(4*l+1:5*l, 1) = update{3};
            
            sol(:, i) = u;
            
        else
            
%-------------------------1st thread-------------------------- 
            temp{1} = u(1:2*l, 1);
            temp{2} = u(2*l+1:4*l, 1);
            temp{3} = u(4*l+1 : 5*l, 1);

            parfor j = 1:3
                
                if j == 1
                    
                    temp{j} = (eye(2*l) - r*G2m2) * temp{j} + deltat*lambda*temp{j};

                end
                
                
%-------------------------2nd thread-------------------------- 

                if j == 2
                    
                    temp{j} = (eye(2*l) - r*G2m) * temp{j} + deltat*lambda*temp{j};

                end
                
%-------------------------3rd thread--------------------------
                if j == 3
                    
                    temp{j} = (eye(l) - r*Gm2) * temp{j} + deltat*lambda*temp{j};


                end
            end


            u(1:2*l, 1) = temp{1};            
            u(2*l+1:4*l, 1) = temp{2};   
            u(4*l+1 : 5*l, 1) = temp{3};            

            
            update{1} = u(1:l, 1);            
            update{2} = u(l+1:3*l, 1);   
            update{3} = u(3*l+1 : 5*l, 1);
            
            
            parfor j = 1:3
                
                if j == 1
                    
                    %update{j} = (eye(l) + r*Gm1) \ update{j};
                    [update{j}, ~, ~]= gmres((eye(l) + r*Gm1), update{j}, 10, 0.001);

                end
                                
%-------------------------2nd thread-------------------------- 

                if j == 2
                    
                    %update{j} = (eye(2*l) + r*G2m) \ update{j};
                    [update{j}, ~, ~]= gmres((eye(2*l) + r*G2m), update{j}, 10, 0.001);

                end
                
%-------------------------3rd thread--------------------------
                if j == 3
                    
                    %update{j} = (eye(2*l) + r*G2m1) \ update{j};
                    [update{j}, ~, ~]= gmres((eye(2*l) + r*G2m1), update{j}, 10, 0.001);

                end
            end 
            
            u(1:l, 1) = update{1};  
            u(l+1:3*l, 1) = update{2};    
            u(3*l+1 : 5*l, 1) = update{3};
            
            sol(:, i) = u;
                                                       
            
        end
        
    end  
    
    
    toc
    delete(pool);
    
    figure;
    surf(xlist,tlist,sol'); 
    title('ASC-N Numerical solution')
    xlabel('Distance x')
    ylabel('Time t')
    
    
end
    