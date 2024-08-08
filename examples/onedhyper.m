

%     deltat = 0.1;
%     deltax = 1;    
%     xrange = 10;
%     alpha = 0.5;
%     time = 51;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 5, time);     
%     
%     [xgrid, tgrid] = meshgrid(xlist, tlist);
%     hfun = sqrt(2)*exp(-xgrid-tgrid);
%     
%     init = [0, 0, 0, 0, 1.5, 1.5, 1.5, 0, 0, 0, 0];
%     bdcnd = 'Dirichlet';
%     sol = solve_hypo(alpha, deltat, deltax, init, time, xlist, tlist, bdcnd, hfun);
    
    %-------------------2d lin eq one side method reachable set---------------%

    deltat = 0.1;
    deltax = 0.4;
    alpha = 0.5;
    xrange = 10;
    time = 5/deltat + 1;
    numofmesh = xrange/deltax + 1;   
    xlist = linspace(0, xrange, numofmesh);
    tlist = linspace(0, 5, time);
    bdcnd = 'Dirichlet';
    h_coef1 = 0.9;
    h_coef2 = 1.1;
    
    init_min = zeros(1, numofmesh);
    init_min(1, 1: numofmesh/5) = 0.5;    
      
    init_max = zeros(1, numofmesh);
    init_max(1, 1: numofmesh/5 ) = 1;  
    
    [xgrid, tgrid] = meshgrid(xlist, tlist);
    
    h_min = h_coef1*sqrt(2)*exp(-xgrid-tgrid);
    h_max = h_coef2*sqrt(2)*exp(-xgrid-tgrid);
    
    
    [sol_min, sol_max] = reach_linhypo(alpha, deltat, deltax, init_min, init_max, time, xlist, tlist, bdcnd, h_min, h_max);
    %sol_min shape = [x, t]
    

    
    %-------------------error for above ---------------%

    deltat = 0.1;
    deltax = 0.4;
    alpha = 0.5;
    xrange = 10;
    trange = 5;
    time = trange/deltat + 1;
    numofmesh = xrange/deltax + 1;   
    xlist = linspace(0, xrange, numofmesh);
    tlist = linspace(0, trange, time);
    bdcnd = 'Dirichlet';
    h_coef_min = 0.9;
    h_coef_max = 1.1;
    
    init_min = zeros(1, numofmesh);      
    init_max = zeros(1, numofmesh)  + 0.000001;
    
    slope_t_min = (sol_min(:,2:end) - sol_min(:,1:end - 1))/deltat;     %first element is index 1 not 0
    slope_t_min = [zeros(numofmesh, 1), slope_t_min];
    slope_x_min = (sol_min(2:end, :) - sol_min(1:end - 1, :))/deltax;
    slope_x_min = [zeros(1, time); slope_x_min];
    
    slope_t_max = (sol_max(:,2:end) - sol_max(:,1:end - 1))/deltat;
    slope_t_max = [zeros(numofmesh, 1), slope_t_max];
    slope_x_max = (sol_max(2:end, :) - sol_max(1:end - 1, :))/deltax;
    slope_x_max = [zeros(1, time); slope_x_max];
    
    [xgrid, tgrid] = meshgrid(xlist, tlist);
    
    h_min = h_coef_min*sqrt(2)*exp(-xgrid-tgrid) - (slope_t_min + alpha*slope_x_min)';
    h_max = h_coef_max*sqrt(2)*exp(-xgrid-tgrid) - (slope_t_max + alpha*slope_x_max)';
    
    
    [err_sol_min, err_sol_max] = reach_linhypo(alpha, deltat, deltax, init_min, init_max, time, xlist, tlist, bdcnd, h_min, h_max);

    plot_2dboxerr(err_sol_min(:,4/deltat + 1), err_sol_max(:,4/deltat + 1), deltax, xlist, 0, [1 0 0]);   
    hold on;
    
    plot_fixlocboxerr(err_sol_min(5, :), err_sol_max(5, :), deltat, tlist, [1 0 0]);
    hold on;
    
    bloated_min = sol_min + err_sol_min;
    bloated_max = sol_max + err_sol_max;
    
    plot_fixlocbox(bloated_min(5, :), bloated_max(5, :), deltat, tlist, '');
    hold on;  
    
    plot_2dboxcolor(bloated_min(:,4/deltat + 1), bloated_max(:,4/deltat + 1), deltax, xlist, 0, '');   
    hold on;
    
% %----------------------------one side method 3dbox reachable set ------------------------------------
% 
%  %this is correct 
%     deltat = 0.1;
%     deltax = 1;
%     alpha = 0.5;
%     xrange = 10;
%     time = 81;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 8, time);
%     bdcnd = 'Dirichlet';
%     init_min = [0, 0, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0];
%     init_max = [0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0];
% %     tic
%     reach_3dlinhypo(alpha, deltat, deltax, xlist, tlist, init_min, init_max, time, bdcnd);
%     toc


