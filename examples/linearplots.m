function linearplots
    deltat = 0.1;
    deltax = 1;
    alpha = 0.5;
    xrange = 20;
    time = 5/deltat + 1;
    numofmesh = xrange/deltax + 1;   
    xlist = linspace(0, xrange, numofmesh);
    tlist = linspace(0, 5, time);
    bdcnd = 'Dirichlet';
    
    
    init_min = zeros(1, numofmesh);
    init_min(1 , 1: numofmesh/5) = 0.5;    
      
    init_max = zeros(1, numofmesh);
    init_max(1, 1: numofmesh/5 ) = 1;  
    
    tic;
    sol_min = solve_hypo(alpha, deltat, deltax, init_min, time, xlist, tlist, bdcnd);
    sol_max = solve_hypo(alpha, deltat, deltax, init_max, time, xlist, tlist, bdcnd);
    toc;

 %--------2nd plot------------------%     
 
 
    deltat = 0.1;
    deltax2 = 0.1;
    alpha = 0.5;
    xrange = 20;
    time = 5/deltat + 1;
    numofmesh = xrange/deltax2 + 1;   
    xlist2 = linspace(0, xrange, numofmesh);
    tlist = linspace(0, 5, time);
    bdcnd = 'Dirichlet';
    
    
   
    init_min2 = zeros(1, numofmesh);
    init_min2(1 , 1: numofmesh/5) = 0.5;    
      
    init_max2 = zeros(1, numofmesh);
    init_max2(1, 1: numofmesh/5 ) = 1;  
    
    sol_min2 = solve_hypo(alpha, deltat, deltax2, init_min2, time, xlist2, tlist, bdcnd);
    sol_max2 = solve_hypo(alpha, deltat, deltax2, init_max2, time, xlist2, tlist, bdcnd);
    
 %--------annotation for linear eq at t = 4s------------------%  
    figure;
    plot_2dboxcolor(sol_min(:,4/deltat + 1), sol_max(:,4/deltat + 1), deltax, xlist, 0, [0 0.5 0.5]);   
    hold on;
    rectangle('Position',[2   0.5  4  0.5]); 
%     t = text(4, 0.73,'Unsafe region');    
    axis([-2 15 0 1.05])
    hold off;
    
%--------annotation for linear eq at t = 4s------------------%  
    figure;
    plot_2dboxcolor(sol_min2(:,4/deltat + 1), sol_max2(:,4/deltat + 1), deltax2, xlist2, 0, [0 0.5 0.5]);   
    hold on;
    rectangle('Position',[2   0.5  4  0.5]); 
%     t = text(4, 0.73,'Unsafe region');    
    axis([-2 15 0 1.05])
    hold off;
end