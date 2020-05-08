function [sol_min, sol_max] = reach_nonlnhypo(deltat, deltax, init_min, init_max, time, xlist, tlist, bdcnd)

u_min = init_min;
u_max = init_max;
sol_min = zeros(length(xlist),length(tlist));
sol_max = zeros(length(xlist),length(tlist));

size(sol_min);
sol_min(:,1) = init_min;
sol_max(:,1) = init_max;

mesh = size(xlist);
m = mesh(2);    %interior mesh points
% disp(m);
u_min_new = zeros(size(init_min));
u_max_new = zeros(size(init_min));


    function f = fun(x) %flux function
        f = x * x;
    end


sol_min = solve_nonlnhypo(deltat, deltax, init_min, time, xlist, tlist, bdcnd);
sol_max = solve_nonlnhypo(deltat, deltax, init_max, time, xlist, tlist, bdcnd);


for i = 1 : time - 1    

    u_min_new = min_nonlnhypo(u_min, u_max, deltat, deltax, xlist);
    
    if strcmp(bdcnd,'Dirichlet')    
        u_min_new(1) = 0;
        u_min_new(m) = 0;
    end    
    u_max_new = max_nonlnhypo(u_min, u_max, deltat, deltax, xlist);
    
    if strcmp(bdcnd,'Dirichlet')    
        u_max_new(1) = 0;
        u_max_new(m) = 0;
    end
    
    for j = 2 : m - 1       %for the purpose of soundness, adding max and min value at each step
        if u_min_new(j) > sol_min(j, i + 1)
            u_min_new(j) = sol_min(j, i + 1);
        end
        
        if u_max_new(j) < sol_max(j, i + 1)
            u_max_new(j) = sol_max(j, i + 1);
        end
    end
    
    sol_min(:,i + 1) = u_min_new;
    sol_max(:,i + 1) = u_max_new;
    
    u_min = u_min_new;
    u_max = u_max_new;
    
 
end
% sol_min = sol_min';
% sol_max = sol_max';

sol_sample1 = solve_nonlnhypo(deltat, deltax, (init_min + init_max)/2, time, xlist, tlist, bdcnd);
sol_sample2 = solve_nonlnhypo(deltat, deltax, (init_max - init_min)/3 + init_min, time, xlist, tlist, bdcnd);
sol_sample3 = solve_nonlnhypo(deltat, deltax, (init_max - init_min)/3 * 2 + init_min, time, xlist, tlist, bdcnd);


% plot_2dbox(sol_min(:,41)', sol_max(:,41)', sol_sample1(:,41), sol_sample2(:,41), sol_sample3(:,41), deltax, xlist);


% %--------------plot safe region for nonlin eq--------------------------------% 
%     hold on;
%     t = text(0, 1.15,'Unsafe region');
%     t.Color = [1 0 0];
%     safe_line = [1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1]; 
%     plot(linspace(-2, 12, 15), safe_line, '--r'); 
% %--------annotation for nonlinear eq------------------%  
%     x = [0.6,0.57];
%     y = [0.9,0.87];
%     annotation('textarrow',x,y,'String','Uncertain')
% 
%     p = [0.35,0.42];
%     q = [0.4,0.4];
%     annotation('textarrow',p,q,'String','Single traces')
% 
%     dim = [.47 .83 .09 .06];
%     annotation('ellipse',dim)
%     
% % --------------plot fixed location reachable set ------------------------%
%     plot_fixlocbox(sol_min(5, :)', sol_max(5, :)', deltat, tlist);
    
    
figure;
surf(xlist,tlist,sol_min') 
title('Numerical solution for min value')
xlabel('Distance x')
ylabel('Time t')

figure;
surf(xlist,tlist,sol_max') 
title('Numerical solution for max value')
xlabel('Distance x')
ylabel('Time t')
end