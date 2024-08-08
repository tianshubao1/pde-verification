%-------------------------1d euler equations reachable set-----------------------------------

    deltat = 0.5;
    deltax = 2;    
    xrange = 20;
    trange = 4;
    gamma = 1.4;
    time = trange/deltat + 1; 
    numofmesh = xrange/deltax + 1;
    h_coef1 = 0.9;
    h_coef2 = 1.1;
    
    
    xlist = linspace(0, xrange, numofmesh);
    tlist = linspace(0, trange, time);  
    
    
    init_min1 = zeros(1, numofmesh) + 15;
    init_min1(1, 1: numofmesh/2) = 14;          
    init_max1 = zeros(1, numofmesh) + 18;
    init_max1(1, 1: numofmesh/2 ) = 20;

    init_min2 = zeros(1, numofmesh) ;
    init_min2(1, 1: numofmesh/2) = 2;         
    init_max2 = zeros(1, numofmesh) +1;
    init_max2(1, 1: numofmesh/2 ) = 4;

    init_min3 = zeros(1, numofmesh) ;
    init_min3(1, 1: numofmesh/2) = 2;          
    init_max3 = zeros(1, numofmesh) + 1;
    init_max3(1, 1: numofmesh/2 ) = 4;

     
    init_min = [init_min1; init_min2; init_min3];       %[3, x, t]
    init_max = [init_max1; init_max2; init_max3];
    
    [xgrid, tgrid] = meshgrid(xlist, tlist);    
    h_min = h_coef1*exp(-xgrid-tgrid);  %[t, x]
    h_max = h_coef2*exp(-xgrid-tgrid); 
    
    
    [sol_min, sol_max] = reach_sys_nonlnhypo(gamma, deltat, deltax, init_min, init_max, time, xlist, tlist, h_min', h_max');
    
    plot_2dboxcolor(sol_min(3,:,4/deltat + 1), sol_max(3, :,4/deltat + 1), deltax, xlist, 0, '');   %energy
    hold on;    
    
    
%-------------------error for above ---------------%

    deltat = 0.5;
    deltax = 2;    
    xrange = 20;
    trange = 4;
    gamma = 1.4;
    time = trange/deltat + 1; 
    numofmesh = xrange/deltax + 1;
    h_coef_min = 0.9;
    h_coef_max = 1.1;
    
    
    xlist = linspace(0, xrange, numofmesh);
    tlist = linspace(0, trange, time);  
    
    
    init_min1 = zeros(1, numofmesh);      
    init_max1 = zeros(1, numofmesh) + 0.02;

    init_min2 = zeros(1, numofmesh) ;       
    init_max2 = zeros(1, numofmesh) +0.02;

    init_min3 = zeros(1, numofmesh) ;        
    init_max3 = zeros(1, numofmesh) + 0.02;

     
    init_min = [init_min1; init_min2; init_min3];       %[3, x, t]
    init_max = [init_max1; init_max2; init_max3];

    slope_t_min = (sol_min(1,:,3:end) - sol_min(1,:,1:end - 2))/(2*deltat);     %first element is index 1 not 0, h = [e(-x-t); 0; 0]
    slope_t_min = squeeze(slope_t_min);
    slope_t_min = [zeros(numofmesh, 1), slope_t_min, zeros(numofmesh, 1)];
    
    slope_x_min = (sol_min(2, 3:end, :) - sol_min(2,1:end - 2, :))/(2*deltax);  %[2, x, t]
    slope_x_min = squeeze(slope_x_min);
    slope_x_min = [zeros(1, time); slope_x_min; zeros(1, time)];
    
    
    slope_t_max = (sol_max(1,:,3:end) - sol_max(1,:,1:end - 2))/(2*deltat);     % only in rho
    slope_t_max = squeeze(slope_t_max);
    slope_t_max = [zeros( numofmesh, 1), slope_t_max, zeros( numofmesh, 1)];
    
    slope_x_max = (sol_max(2,3:end, :) - sol_max(2,1:end - 2, :))/deltax;
    slope_x_max = squeeze(slope_x_max);
    slope_x_max = [zeros(1, time); slope_x_max; zeros(1, time)];
    
    [xgrid, tgrid] = meshgrid(xlist, tlist);
    
     
    h_min1 = h_coef_min*exp(-xgrid-tgrid);    %[t, x]
    h_max1 = h_coef_max*exp(-xgrid-tgrid); 
    
%     h_min1 = h_coef_min*exp(-xgrid-tgrid) - (slope_t_min + slope_x_min)';     %[t, x]
%     h_max1 = h_coef_max*exp(-xgrid-tgrid) - (slope_t_max + slope_x_max)';
    
    h_min = min(h_min1, h_max1);
    h_max = max(h_min1, h_max1);
    
    
    [err_sol_min, err_sol_max] = reach_sys_nonlnhypo(gamma, deltat, deltax, init_min, init_max, time, xlist, tlist, h_min', h_max');
    
    plot_2dboxerr(err_sol_min(3,:,4/deltat + 1), err_sol_max(3, :,4/deltat + 1), deltax, xlist, 0, '');   %energy plotting
    hold on;       
  
    
    
    bloated_min = sol_min + err_sol_min;
    bloated_max = sol_max + err_sol_max;
    save('nonsys_bloat_min.mat','bloated_min');
    save('nonsys_bloat_max.mat','bloated_max');
    
    plot_2dboxcolor(bloated_min(3,:,4/deltat + 1), bloated_max(3, :,4/deltat + 1), deltax, xlist, 0, '');   
    hold on;
        
    bloated_min_engy = squeeze(bloated_min(3, :, :))';
    bloated_max_engy = squeeze(bloated_max(3, :, :))';
    plot_3dbox(bloated_min_engy, bloated_max_engy, deltax, deltat, xlist, tlist);
    title('Bloated Continuous Reachable Sets for Energy')   
    
    
    bloated_min_mmt = squeeze(bloated_min(2, :, :))';
    bloated_max_mmt = squeeze(bloated_max(2, :, :))';
    plot_fixlocbox(bloated_min_mmt(4, :), bloated_max_mmt(4, :), deltat, tlist, '');
    hold on;    
    