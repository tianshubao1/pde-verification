function reach_3d_boundary_control(alpha, deltat, deltax, init_min, init_max, time, xlist, tlist, lambda_1, lambda_2)

    
    u = init_min;
    sol_min = zeros(length(xlist),length(tlist));
    size(sol_min);
    size2 = size(xlist);
    m = size2(2);
    sol(:,1) = init_min;
    
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
            
            u_1 = u - deltat/deltax*alpha*([u(2:length(u)), control] - u) + deltat * lambda_1 * input;  %update range of lambda here
            u_2 = u - deltat/deltax*alpha*([u(2:length(u)), control] - u) + deltat * lambda_2 * input;  %update range of lambda here
            
            for j = 1 : m
                u(j) = min(u_1(j), u_2(j));
            end            
            
            sol_min(:,i) = u;    
        end
    end

    
    
    u = init_max;
    sol_max = zeros(length(xlist),length(tlist));
    size(init_max);
    size2 = size(xlist);
    m = size2(2);
    sol_max(:,1) = init_max;
    
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
            
            u_1 = u - deltat/deltax*alpha*([u(2:length(u)), control] - u) + deltat * lambda_1 * input;  %update range of lambda here
            u_2 = u - deltat/deltax*alpha*([u(2:length(u)), control] - u) + deltat * lambda_2 * input;  %update range of lambda here
            
            for j = 1 : m
                u(j) = max(u_1(j), u_2(j));
            end            
            
            sol_max(:,i) = u;    
        end
    end    
    
    plot_3dbox(sol_min', sol_max', deltax, deltat, xlist, tlist);


end