function [sol_min, sol_max] = reach_sys_nonlnhypo(gamma, deltat, deltax, init_min, init_max, time, xlist, tlist)
    tic;
    u_min = init_min;
    u_max = init_max;
    mesh = size(xlist);
    m = mesh(2);    %mesh points

%     sol_min = solve_sys_nonlnhypo(gamma, deltat, deltax, init_min, time, xlist, tlist);
%     sol_max = solve_sys_nonlnhypo(gamma, deltat, deltax, init_max, time, xlist, tlist);

    sol_min = zeros(3, m, time);
    sol_max = zeros(3, m, time);
    
    for i = 1 : time - 1    

        u_min_new = min_sys_nonlnhypo(u_min, u_max, deltat, deltax, xlist, gamma);
        u_max_new = max_sys_nonlnhypo(u_min, u_max, deltat, deltax, xlist, gamma); 
% 
%         for j = 1 : m       %for the purpose of soundness, adding max and min value at each step
%             for k = 1 : 3
% 
%                 if u_min_new(k, j) > sol_min(k, j, i + 1)
%                     u_min_new(k, j) = sol_min(k, j, i + 1);
%                 end
% 
%                 if u_max_new(k, j) < sol_max(k, j, i + 1)
%                     u_max_new(k, j) = sol_max(k, j, i + 1);
%                 end
% 
%             end
%         end
% 
        sol_min(:, :, i + 1) = u_min_new;
        sol_max(:, :, i + 1) = u_max_new;

        u_min = u_min_new;
        u_max = u_max_new;

        disp(u_min);
        disp(u_max);
    end

    toc;
    
    
%     plot_2dboxbald(sol_min(1, :,41), sol_max(1, :,41), deltax, xlist, 0);      % t = 4      
%     plot_fixlocbox(sol_min(3, 5, :), sol_max(3, 5, :), deltat, tlist);      % x = 4

%---------------rho plotting--------------------------------
    sol_min_rho = squeeze(sol_min(1, :, :))';
    sol_max_rho = squeeze(sol_max(1, :, :))';
    plot_3dbox(sol_min_rho, sol_max_rho, deltax, deltat, xlist, tlist);
%     title('3-Dimension Reachable Sets for Density')
    
%---------------energy plotting--------------------------------    
    sol_min_engy = squeeze(sol_min(3, :, :))';
    sol_max_engy = squeeze(sol_max(3, :, :))';
    plot_3dbox(sol_min_engy, sol_max_engy, deltax, deltat, xlist, tlist);
%     title('3-Dimension Reachable Sets for Energy')
    
end