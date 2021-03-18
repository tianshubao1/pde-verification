function nonlinearplots
    deltat = 0.5;       % dont change this, use large time step
    deltax1 = 5;    % delta x = 1, 2, 4
    xrange = 100;
    time = 5/deltat + 1;    % simluate from 0 to 5 seconds, too many optimizations will lead to over approximation
    xlist1 = linspace(0, xrange, xrange/deltax1 + 1);
    tlist = linspace(0, 5, time);    
    bdcnd = 'Dirichlet';
    
    m = xrange/deltax1 + 1;
    init_min = zeros(1, m);
    init_min(1, 1: m * 0.3) = 0.5;
    init_max = zeros(1, m);
    init_max(1, 1: m * 0.3) = 1;  
    
    tic
%     [sol_min1, sol_max1] = reach_nonlnhypo(deltat, deltax1, init_min', init_max', time, xlist1, tlist, bdcnd);
    toc
    
    
    deltat = 0.5;       % dont change this, use large time step
    deltax2 = 2;    % delta x = 1, 2, 4
    xrange = 100;
    time = 5/deltat + 1;    % simluate from 0 to 5 seconds, too many optimizations will lead to over approximation
    xlist2 = linspace(0, xrange, xrange/deltax2 + 1);
    tlist = linspace(0, 5, time);    
    bdcnd = 'Dirichlet';
    
    m = xrange/deltax2 + 1;
    init_min = zeros(1, m);
    init_min(1, 1: m * 0.3) = 0.5;
    init_max = zeros(1, m);
    init_max(1, 1: m * 0.3) = 1;   
  
%     tic;
%     [sol_min2, sol_max2] = reach_nonlnhypo(deltat, deltax2, init_min', init_max', time, xlist2, tlist, bdcnd);   
%     toc;
    
    deltat = 0.5;       % dont change this, use large time step
    deltax3 = 1;    % delta x = 1, 2, 4
    xrange = 100;
    time = 5/deltat + 1;    % simluate from 0 to 5 seconds, too many optimizations will lead to over approximation
    xlist3 = linspace(0, xrange, xrange/deltax3 + 1);
    tlist = linspace(0, 5, time);    
    bdcnd = 'Dirichlet';
    
    m = xrange/deltax3 + 1;
    init_min = zeros(1, m);
    init_min(1, 1: m * 0.3) = 0.5;
    init_max = zeros(1, m);
    init_max(1, 1: m * 0.3) = 1;   
    
%     tic;
    [sol_min3, sol_max3] = reach_nonlnhypo(deltat, deltax3, init_min', init_max', time, xlist3, tlist, bdcnd);  
%     toc;
    
    
    figure;   
   
    
    plot_2dboxcolor(sol_min3(:,4/deltat + 1), sol_max3(:,4/deltat + 1), deltax3, xlist3, 0, [0 0.5 0.5]);   % plot t = 4s   
    hold on;
%     plot_2dboxcolor(sol_min2(:,4/deltat + 1), sol_max2(:,4/deltat + 1), deltax2, xlist2, 0, [0 0.5 0.5]);   % plot t = 4s
%     hold on;    
%     plot_2dboxbald(sol_min1(:,4/deltat + 1)', sol_max1(:,4/deltat + 1)', deltax1, xlist1, 0); % plot t = 4s
%     hold on;
    x=[1 31 32 2 1];
    y=[0.5 0.5 1 1 0.5];
    plot(x,y, 'k')
    axis([-5 60 0 3])
%     rectangle('Position',[1   0.5  31  0.5]); 
%     t = text(4, 0.73,'Unsafe region');    
    hold off;
    
%     plot_3dbox(sol_min3', sol_max3', deltax, deltat, xlist, tlist);
%     axis([-5 60 0 3])    
end