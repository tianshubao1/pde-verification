function sol = solve_nonlnhypo(deltat, deltax, init, time, xlist, tlist, bdcnd)

    function f = fun(x) %flux function
        f = 1/2* x * x;
    end


    if strcmp(bdcnd,'Dirichlet')
        u = init';
        sol = zeros(length(xlist),length(tlist));
        size(sol);
        sol(:,1) = init';
        mesh = size(xlist);
        m = mesh(2);
        disp(m);
        u_new = zeros(size(u));



        for i = 1 : time - 1

            for j = 1 : m   %number of mesh points, richtmyer method

                if j == 1   %left boundary

                    u_1 = (u(j + 1) + u(j))/2 - deltat/(2 * deltax) * (fun(u(j + 1)) - fun(u(j)));  %u_+1/2
                    u_2 = (u(j) + 0)/2 - deltat/(2 * deltax) * (fun(u(j)) - fun(0));  %u_-1/2

                    u_new(j) = u(j) - deltat/deltax * (fun(u_1) - fun(u_2));     %f(u_i+1/2)-f(u_i-1/2)
                    u_new(j) = 0;

                elseif j == m   %right boundary

                    u_1 = (0 + u(j))/2 - deltat/(2 * deltax) * (fun(0) - fun(u(j)));  %u_+1/2
                    u_2 = (u(j) + u(j - 1))/2 - deltat/(2 * deltax) * (fun(u(j)) - fun(u(j - 1)));  %u_-1/2

                    u_new(j) = u(j) - deltat/deltax * (fun(u_1) - fun(u_2));     %f(u_i+1/2)-f(u_i-1/2)
                    u_new(j) = 0;

                else    %in the middle
        
                    u_1 = (u(j + 1) + u(j))/2 - deltat/(2 * deltax) * (fun(u(j + 1)) - fun(u(j)) );  %u_+1/2
                    u_2 = (u(j) + u(j - 1))/2 - deltat/(2 * deltax) * (fun(u(j)) - fun(u(j - 1)) );  %u_-1/2
                    u_new(j) = u(j) - deltat/deltax * (fun(u_1) - fun(u_2));     %f(u_i+1/2)-f(u_i-1/2)
                end

            end
            %disp(u);
            sol(:,i + 1) = u_new;
            u = u_new;
        end
        
        
        
        
        
        
    elseif strcmp(bdcnd,'Neumann') 
        
        u = init;
        sol = zeros(length(xlist),length(tlist));
        size(sol);
        sol(:,1) = init;
        mesh = size(xlist);
        m = mesh(2) - 2;
        disp(m);
        u_new = zeros(size(u));


        for i = 1 : time - 1

            for j = 1 : m   %number of mesh points, richtmyer method

                if j == 1   %left boundary
                    left_ghost = u(j + 1) - 2*deltax*1;                 %u_x(0) = 1
                    
                    u_1 = (u(j + 1) + u(j))/2 - deltat/(2 * deltax) * (fun(u(j + 1)) - fun(u(j)));  %u_+1/2
                    u_2 = (u(j) + left_ghost)/2 - deltat/(2 * deltax) * (fun(u(j)) - fun(left_ghost));  %u_-1/2

                    u_new(j) = u(j) - deltat/deltax * (fun(u_1) - fun(u_2));     %f(u_i+1/2)-f(u_i-1/2)

                elseif j == m   %right boundary

                    u_1 = (0 + u(j))/2 - deltat/(2 * deltax) * (fun(0) - fun(u(j)));  %u_+1/2
                    u_2 = (u(j) + u(j - 1))/2 - deltat/(2 * deltax) * (fun(u(j)) - fun(u(j - 1)));  %u_-1/2

                    u_new(j) = u(j) - deltat/deltax * (fun(u_1) - fun(u_2));     %f(u_i+1/2)-f(u_i-1/2)
                    u_new(j) = 0;

                else    %in the middle
          
                    u_1 = (u(j + 1) + u(j))/2 - deltat/(2 * deltax) * (fun(u(j + 1)) - fun(u(j)));  %u_+1/2
                    u_2 = (u(j) + u(j - 1))/2 - deltat/(2 * deltax) * (fun(u(j)) - fun(u(j - 1)));  %u_-1/2
                    u_new(j) = u(j) - deltat/deltax * (fun(u_1) - fun(u_2));     %f(u_i+1/2)-f(u_i-1/2)
                end

            end
            %disp(u);
            sol(:,i + 1) = u_new;
            u = u_new;
            
        end 
            
            
            
    elseif strcmp(bdcnd,'Robin') 
        
        u = init;
        sol = zeros(length(xlist),length(tlist));
        size(sol);
        sol(:,1) = init;
        mesh = size(xlist);
        m = mesh(2) - 2;
        disp(m);
        u_new = zeros(size(u));


        for i = 1 : time - 1

            for j = 1 : m   %number of mesh points, richtmyer method

                if j == 1   %left boundary

                    u_1 = (u(j + 1) + u(j))/2 - deltat/(2 * deltax) * (fun(u(j + 1)) - fun(u(j)));  %u_+1/2
                    u_2 = (u(j) + 0)/2 - deltat/(2 * deltax) * (fun(u(j)) - fun(0));  %u_-1/2

                    u_new(j) = u(j) - deltat/deltax * (fun(u_1) - fun(u_2));     %f(u_i+1/2)-f(u_i-1/2)
                    u_new(j) = 0;

                elseif j == m   %right 
                    
                    right_ghost = u(j - 1) + 2*deltax*(1 - u(j));                 %u_x(L) + u(L) = 1
                    u_1 = (right_ghost + u(j))/2 - deltat/(2 * deltax) * (fun(right_ghost) - fun(u(j)));  %u_+1/2
                    u_2 = (u(j) + u(j - 1))/2 - deltat/(2 * deltax) * (fun(u(j)) - fun(u(j - 1)));  %u_-1/2

                    u_new(j) = u(j) - deltat/deltax * (fun(u_1) - fun(u_2));     %f(u_i+1/2)-f(u_i-1/2)


                else    %in the middle
         
                    u_1 = (u(j + 1) + u(j))/2 - deltat/(2 * deltax) * (fun(u(j + 1)) - fun(u(j)));  %u_+1/2
                    u_2 = (u(j) + u(j - 1))/2 - deltat/(2 * deltax) * (fun(u(j)) - fun(u(j - 1)));  %u_-1/2
                    u_new(j) = u(j) - deltat/deltax * (fun(u_1) - fun(u_2));     %f(u_i+1/2)-f(u_i-1/2)
                end

            end
            %disp(u);
            sol(:,i + 1) = u_new;
            u = u_new;
        end
    else
        error('No boundary condition matches the input');
    end


    surf(xlist,tlist,sol') 
    title('Numerical solution')
    xlabel('Distance x')
    ylabel('Time t')

end