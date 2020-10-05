function sol = solve_boundary_control(alpha, deltat, deltax, init, time, xlist, tlist, lambda)
    % choose alpha < 0
    u = init;
    sol = zeros(length(xlist),length(tlist));
    size(sol);
    size2 = size(xlist);
    m = size2(2);
    sol(:,1) = init;
    if alpha > 0
        for i = 1 : time
                %disp(size(u));
            u = u - deltat/deltax*alpha*(u - [0, u(1:(length(u) - 1))]);
            u(1) = 0;

            sol(:,i) = u;    
        end
    end

    if alpha < 0
        for i = 1 : time
            control = -sum(u)/(m - 1); %boundary control
            input = zeros(1, m);
            
            for j = 1 : m
                input(j) = (j - 1)/10 * u(1);
            end
            
            u = u - deltat/deltax*alpha*([u(2:length(u)), control] - u) + deltat * lambda * input;  %update range of lambda here
  
            sol(:,i) = u;    
        end
    end
        
    
            


    figure;
    surf(xlist,tlist,sol') 
    title('Solution with Boundary Control')
    xlabel('Distance x')
    ylabel('Time t')
end