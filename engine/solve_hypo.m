function sol = solve_hypo(alpha, deltat, deltax, init, time, xlist, tlist, bdcnd)
%one side scheme
%     if bdcnd == 'Dirichlet' 
    if strcmp(bdcnd,'Dirichlet')
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
                u(m) = 0;
                sol(:,i) = u;    
            end
        end

        if alpha < 0
            for i = 1 : time
                u = u - deltat/deltax*alpha*([u(2:length(u)),0] - u);
                u(1) = 0;
                u(m) = 0;            
                sol(:,i) = u;    
            end
        end
        
        
%     elseif bdcnd == 'Neumann'

    elseif strcmp(bdcnd,'Neumann')        
        u = init;
        sol = zeros(length(xlist),length(tlist));
        size(sol);
        size2 = size(xlist);
        m = size2(2);
        sol(:,1) = init;
        if alpha > 0
            for i = 1 : time
                %disp(size(u));
                u(2:length(u)) = u(2:length(u)) - deltat/deltax*alpha*(u(2:length(u)) - u(1:(length(u) - 1)) );     %ignore left boundary point and replace it by derivative, u_x = 1
                u(1) = u(1) - deltat*alpha; 
                u(m) = 0;
                sol(:,i) = u;    
            end
        end

        if alpha < 0
            for i = 1 : time
                u(1:(length(u) - 1)) = u(1:(length(u) - 1)) - deltat/deltax*alpha*(u(2:length(u)) - u(1:(length(u) - 1)));   %ignore right boundary point and replace it by derivative, u_x = 1
                u(1) = 0;
                u(m) = u(m) - deltat*alpha;           
                sol(:,i) = u;    
            end
        end
        
        
%     elseif bdcnd == 'Robin'
    elseif strcmp(bdcnd,'Robin')          
        u = init;
        sol = zeros(length(xlist),length(tlist));
        size(sol);
        size2 = size(xlist);
        m = size2(2);
        sol(:,1) = init;
        if alpha > 0
            for i = 1 : time
                %disp(size(u));
                u(2:length(u)) = u(2:length(u)) - deltat/deltax*alpha*(u(2:length(u)) - u(1:(length(u) - 1)) );     %ignore left boundary point and replace it by derivative, u_x = 1 - u(m)
                u(1) = u(1) - deltat*alpha *(0 - u(1)); 
                u(m) = 0;
                sol(:,i) = u;    
            end
        end

        if alpha < 0
            for i = 1 : time
                u(1:(length(u) - 1)) = u(1:(length(u) - 1)) - deltat/deltax*alpha*(u(2:length(u)) - u(1:(length(u) - 1)));   %ignore right boundary point and replace it by derivative, u_x = 1 - u(m)
                u(1) = 0;
                u(m) = u(m) - deltat*alpha *(0 - u(m));           
                sol(:,i) = u;    
            end
        end
        
    else
        error('No boundary condition matches the input');
    end

    
    
    
    
    
%     figure;
%     surf(xlist,tlist,sol') 
%     title('Numerical solution')
%     xlabel('Distance x')
%     ylabel('Time t')


end
    