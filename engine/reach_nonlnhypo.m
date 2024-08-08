function [sol_min, sol_max] = reach_nonlnhypo(deltat, deltax, init_min, init_max, time, xlist, tlist, bdcnd, h_min, h_max)

u_min = init_min;
u_max = init_max;
sol_min = zeros(length(xlist),length(tlist));   %[11, 51] [x, t]
sol_max = zeros(length(xlist),length(tlist));

size(sol_min);
sol_min(:,1) = init_min; %[x, t]
sol_max(:,1) = init_max;

mesh = size(xlist);
m = mesh(2);    %interior mesh points

h_min = h_min';
h_max = h_max';

    function f = fun(x) %flux function
        f = 1/2 * x * x;
    end

for i = 1 : time - 1    
    
    h_mincurr = min(h_min(:, i), h_max(:, i));  %[x, t]
    h_maxcurr = max(h_min(:, i), h_max(:, i));
    u_min1= u_min;
    u_max1 = u_max;
    u_min = min(u_min1, u_max1); 
    u_max = max(u_min1, u_max1); 
    
%     u_min_new = min_nonlnhypo_local(u_min, u_max, deltat, deltax, xlist, h_mincurr, h_maxcurr);   % local optmization   
%     u_max_new = max_nonlnhypo_local(u_min, u_max, deltat, deltax, xlist, h_mincurr, h_maxcurr);
    u_min_new = min_nonlnhypo_global(u_min, u_max, deltat, deltax, xlist, h_mincurr, h_maxcurr);   % whole layer optmization 
    u_max_new = max_nonlnhypo_global(u_min, u_max, deltat, deltax, xlist, h_mincurr, h_maxcurr);
    
    if strcmp(bdcnd,'Dirichlet')    
        u_max_new(1) = 0;
        u_max_new(m) = 0;
    end
    if strcmp(bdcnd,'Dirichlet')    
        u_min_new(1) = 0;
        u_min_new(m) = 0;
    end 
    
    sol_min(:,i + 1) = u_min_new;
    sol_max(:,i + 1) = u_max_new;
    
    u_min = u_min_new;
    u_max = u_max_new;
    
 
end

    
%     for j = 2 : m - 1       %for the purpose of soundness, adding max and min value at each step
%         if u_min_new(j) > sol_min(j, i + 1)
%             u_min_new(j) = sol_min(j, i + 1);
%         end
%         
%         if u_max_new(j) < sol_max(j, i + 1)
%             u_max_new(j) = sol_max(j, i + 1);
%         end
%     end


% sol_min = sol_min';
% sol_max = sol_max';

% sol_sample1 = solve_nonlnhypo(deltat, deltax, (init_min + init_max)/2, time, xlist, tlist, bdcnd);
% sol_sample2 = solve_nonlnhypo(deltat, deltax, (init_max - init_min)/3 + init_min, time, xlist, tlist, bdcnd);
% sol_sample3 = solve_nonlnhypo(deltat, deltax, (init_max - init_min)/3 * 2 + init_min, time, xlist, tlist, bdcnd);


% plot_2dboxbald(sol_min(:,4/deltat + 1)', sol_max(:,4/deltat + 1)', deltax, xlist, 0); % plot t = 4s

% hold on;
% rectangle('Position',[2   0.5  2.3  0.5], 'FaceColor',[0 0 1], 'EdgeColor',[0 0 1]); 
% hold;
% T = [0,0; 2,0.5; 2,1];
% triplot(T);
% hold;
    
%--------------plot safe region for nonlin eq--------------------------------% 
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
% %     p = [0.35,0.42];
% %     q = [0.4,0.4];
% %     annotation('textarrow',p,q,'String','Single traces')
% 
%     dim = [.47 .83 .09 .06];
%     annotation('ellipse',dim)
%     
% % --------------plot fixed location reachable set ------------------------%
%     plot_fixlocbox(sol_min(5, :)', sol_max(5, :)', deltat, tlist);
%     
    
% figure;
% surf(xlist,tlist,sol_min') 
% title('Numerical solution for min value')
% xlabel('Distance x')
% ylabel('Time t')
% 
% figure;
% surf(xlist,tlist,sol_max') 
% title('Numerical solution for max value')
% xlabel('Distance x')
% ylabel('Time t')
end