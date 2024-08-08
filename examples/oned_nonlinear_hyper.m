%-------------------nonlin eq lax-wdf method min_function-------------------------------------------------    
    deltat = 0.4;
    deltax = 1;    
    xrange = 10;
    trange = 4;
    numofmesh = xrange/deltax + 1; 
    time = trange/deltat + 1;
    h_coef1 = 0.9;
    h_coef2 = 1.1;
    
    xlist = linspace(0, xrange, numofmesh);
    tlist = linspace(0, trange, time);    
    
    init_min = zeros(1, numofmesh);
    init_min(1, 1: numofmesh/2) = 0.5;  
    
    init_max = zeros(1, numofmesh) + 0.01;     %[t, x]
    init_max(1, 1: numofmesh/2 ) = 1;  
    
    [xgrid, tgrid] = meshgrid(xlist, tlist);    
    h_min = h_coef1*exp(-xgrid-tgrid);  %[t, x]
    h_max = h_coef2*exp(-xgrid-tgrid); 
    tic;
    [sol_min, sol_max] = reach_nonlnhypo(deltat, deltax, init_min, init_max, time, xlist, tlist, 'Dirichlet', h_min, h_max);
    toc;

    plot_2dboxcolor(sol_min(:,4/deltat + 1), sol_max(:,4/deltat + 1), deltax, xlist, 0, '');   
    hold on;
    
    
%-------------------error for above ---------------%
%     deltat = 0.5;
%     deltax = 1;
%     xrange = 20;
%     trange = 4;
    time = trange/deltat + 1;
    numofmesh = xrange/deltax + 1;   
    xlist = linspace(0, xrange, numofmesh);
    tlist = linspace(0, trange, time);
    bdcnd = 'Dirichlet';
    h_coef_min = 0.9;
    h_coef_max = 1.1;
    
    init_min = zeros(1, numofmesh);      
    init_max = zeros(1, numofmesh) + 0.01;
    
    slope_t_min = (sol_min(:,2:end) - sol_min(:,1:end - 1))/(deltat);     %first element is index 1 not 0
    slope_t_min = [zeros(numofmesh, 1), slope_t_min];
    slope_x_min = (sol_min(2:end, :) - sol_min(1:end - 1, :))/(deltax);
    slope_x_min = [zeros(1, time); slope_x_min];
    
    slope_t_max = (sol_max(:,2:end) - sol_max(:,1:end - 1))/deltat;
    slope_t_max = [zeros(numofmesh, 1), slope_t_max];
    slope_x_max = (sol_max(2:end, :) - sol_max(1:end - 1, :))/deltax;
    slope_x_max = [zeros(1, time); slope_x_max];
    
    [xgrid, tgrid] = meshgrid(xlist, tlist);
    
    
%     h_min = h_coef_min*exp(-xgrid-tgrid);% - (slope_t_min + sol_min*slope_x_min/2)';     %[t, x]
%     h_max = h_coef_max*exp(-xgrid-tgrid);% - (slope_t_max + sol_max*slope_x_max/2)';
    h_min = h_coef_min*exp(-xgrid-tgrid) - (slope_t_min + sol_min.*slope_x_min/2)';     %[t, x]
    h_max = h_coef_max*exp(-xgrid-tgrid) - (slope_t_max + sol_max.*slope_x_max/2)';
    
    [err_sol_min, err_sol_max] = reach_nonlnhypo(deltat, deltax, init_min, init_max, time, xlist, tlist, 'Dirichlet', h_min, h_max);
    
    plot_2dboxerr(err_sol_min(:,4/deltat + 1), err_sol_max(:,4/deltat + 1), deltax, xlist, 0, '');   
    hold on;    
    
    
    bloated_min = sol_min + err_sol_min;
    bloated_max = sol_max + err_sol_max;
    
    plot_2dboxcolor(bloated_min(:,4/deltat + 1), bloated_max(:,4/deltat + 1), deltax, xlist, 0, '');   
    hold on;
    
    plot_fixlocbox_berger(bloated_min(6, :), bloated_max(6, :), deltat, tlist, '');    % x = 5
    hold on;
    