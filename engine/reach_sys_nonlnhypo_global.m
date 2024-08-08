function [sol_min, sol_max] = reach_sys_nonlnhypo_global(gamma, deltat, deltax, init_min, init_max, time, xlist, tlist, h_min, h_max)
    tic;
    u_min = init_min;
    u_max = init_max;
    mesh = size(xlist);
    m = mesh(2);    %mesh points

    sol_min = zeros(3, m, time);        %[3, x, t]
    sol_max = zeros(3, m, time);
    
    for i = 1 : time - 1    

%         u_min_new = min_sys_nonlnhypo(u_min, u_max, deltat, deltax, xlist, gamma, h_min(:, i), h_max(:, i));
%         u_max_new = max_sys_nonlnhypo(u_min, u_max, deltat, deltax, xlist, gamma, h_min(:, i), h_max(:, i)); 
        u_min_new = min_sys_nonlnhypo_global(u_min, u_max, deltat, deltax, xlist, gamma, h_min(:, i), h_max(:, i));
        u_max_new = max_sys_nonlnhypo_global(u_min, u_max, deltat, deltax, xlist, gamma, h_min(:, i), h_max(:, i)); 
 
        u_min_new(1, 1) = 0.5;      %left boundary
        u_min_new(2, 1) = 0.5;
        u_min_new(3, 1) = 0.5;       

        u_min_new(1, end) = 0.5;
        u_min_new(2, end) = 0.5;
        u_min_new(3, end) = 0.5;  
        
        u_max_new(1, 1) = 0.8;      %right boundary
        u_max_new(2, 1) = 0.8;
        u_max_new(3, 1) = 0.8;
        u_max_new(1, end) = 0.8;
        u_max_new(2, end) = 0.8;
        u_max_new(3, end) = 0.8;
        
        
        disp(i);
%         disp(u_min_new);
%         disp(u_max_new);

        sol_min(:, :, i + 1) = u_min_new;
        sol_max(:, :, i + 1) = u_max_new;
 
%         sol_max(1, 1, i + 1) = 0.01;    %left boundary max value
%         sol_max(2, 1, i + 1) = 0.01;
%         sol_max(3, 1, i + 1) = 0.01;
%         
%         sol_max(1, end, i + 1) = 0.01;    %right boundary max value
%         sol_max(2, end, i + 1) = 0.01;
%         sol_max(3, end, i + 1) = 0.01;
%         
%         sol_max(1, 1, i + 1) = 0.01;    %left boundary max value
%         sol_max(2, 1, i + 1) = 0.01;
%         sol_max(3, 1, i + 1) = 0.01;
%         
%         sol_max(1, end, i + 1) = 0.01;    %right boundary max value
%         sol_max(2, end, i + 1) = 0.01;
%         sol_max(3, end, i + 1) = 0.01;

        u_min = sol_min(:, :, i + 1);
        u_max =  sol_max(:, :, i + 1);

    end

%     toc;
    
    
%     plot_2dboxbald(sol_min(1, :,41), sol_max(1, :,41), deltax, xlist, 0);      % t = 4      
%     plot_fixlocbox(sol_min(3, 5, :), sol_max(3, 5, :), deltat, tlist);      % x = 4
% 
% %---------------rho plotting--------------------------------
%     sol_min_rho = squeeze(sol_min(1, :, :))';
%     sol_max_rho = squeeze(sol_max(1, :, :))';
%     plot_3dbox(sol_min_rho, sol_max_rho, deltax, deltat, xlist, tlist);
% %     title('3-Dimension Reachable Sets for Density')
%     
% %---------------energy plotting--------------------------------    
%     sol_min_engy = squeeze(sol_min(3, :, :))';
%     sol_max_engy = squeeze(sol_max(3, :, :))';
%     plot_3dbox(sol_min_engy, sol_max_engy, deltax, deltat, xlist, tlist);
% %     title('3-Dimension Reachable Sets for Energy')
%     
end