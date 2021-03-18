function sysofhypo
%     deltat = 0.1;
%     deltax = 1;    
%     xrange = 10;
%     alpha = 0.5;
%     time = 51;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 5, time);     
%     init = [0, 0, 0, 0, 1.5, 1.5, 1.5, 0, 0, 0, 0];
%     bdcnd = 'Neumann';
%     sol = solve_hypo(alpha, deltat, deltax, init, time, xlist, tlist, bdcnd);
    
%-------------------------solve lin sys using lax-frdch method-------------------------------------------     

%     A = [-0.5, 0; 0, -0.5];
%     
%     deltat = 0.1;
%     deltax = 1;    
%     xrange = 10;
%     time = 51;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 5, time);     
%     init1 = [0, 0, 0, 0, 1.5, 1.5, 1.5, 0, 0, 0, 0];
%     init2 = [0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0];
%     init = [init1; init2];
% %     init = init * V; %transformation
%     solve_jd_hypo(A, deltat, deltax, init, time, xlist, tlist);

%-------------------------solve lin eq using lax-frdch method-------------    
%     A = -0.05;
%     
%     deltat = 0.1;
%     deltax = 0.1;    
%     xrange = 1;
%     time = 51;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 5, time);     
%     init = [0, 1, 0, 0, 1.5, 1.5, 1.5, 0, 0, 0, 0];
% %     init = init * V; %transformation
%     solve_jd_hypo(A, deltat, deltax, init, time, xlist, tlist);
%-------------------------solve nonlin eq using lax-wdf method-------------------------------------------

%     deltat = 0.1;
%     deltax = 1;    
%     xrange = 10;
%     time = 101;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 50, time);
%     bdcnd = 'Neumann';
%     init = [0, 1.2, 1.2, 0, 0, 0, 0, 0, 0, 0, 0];  
% 
% 
% 
%     sol = solve_nonlnhypo(deltat, deltax, init, time, xlist, tlist, bdcnd);


%-------------------2d lin eq one side method reachable set---------------%

%     deltat = 0.1;
%     deltax = 1;
%     alpha = 0.5;
%     xrange = 20;
%     time = 5/deltat + 1;
%     numofmesh = xrange/deltax + 1;   
%     xlist = linspace(0, xrange, numofmesh);
%     tlist = linspace(0, 5, time);
%     bdcnd = 'Dirichlet';
%     
%     
%     init_min = zeros(1, numofmesh);
%     init_min(1, 1: numofmesh/5) = 0.5;    
%       
%     init_max = zeros(1, numofmesh);
%     init_max(1, 1: numofmesh/5 ) = 1;  
%     
%     
%     [sol_min, sol_max] = reach_linhypo(alpha, deltat, deltax, init_min, init_max, time, xlist, tlist, bdcnd);


% %----------------------------one side method 3dbox reachable set ------------------------------------
% 
%  this is correct 
%     deltat = 0.1;
%     deltax = 1;
%     alpha = 0.5;
%     xrange = 10;
%     time = 81;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 8, time);
%     bdcnd = 'Robin';
%     init_min = [0, 0, 0.5, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0];
%     init_max = [0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0];
% %     tic
%     reach_3dlinhypo(alpha, deltat, deltax, xlist, tlist, init_min, init_max, time, bdcnd);
%     toc
%------------------------lin eq lax-frdch method reachable set-------------

%     deltat = 0.1;
%     deltax = 0.2;
%     xrange = 10;
%     time = 51;    
%     xlist = linspace(0, xrange, 51);
%     tlist = linspace(0, 5, time);
%     A = [0.5, 1; 0.5, 0];
%     
%     init_min1 = [0, 0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0,0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0,0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0,0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0,0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0];
%     init_min2 = [0, 0.4, 0.4, 0.5, 0.6, 0, 0, 0.2, 0, 0, 0, 0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0, 0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0, 0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0, 0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0];
%     init_min = [init_min1; init_min2];   
%     
%     init_max1 = [0, 0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0, 0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0, 0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0, 0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0, 0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0];
%     init_max2 = [0, 0.4, 0.4, 0.5, 0.6, 0, 0, 0.2, 0, 0, 0, 0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0, 0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0, 0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0, 0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0];
%     init_max = [init_max1; init_max2]; 
% 
%     [sol_min, sol_max] = reach_linhyposys(A, deltat, deltax, init_min, init_max, time, xlist, tlist);

%----------------------------lin eq lax-fredrich method 3dbox plot------------------------------------

%     deltat = 0.1;
%     deltax = 1;
%     A = [1.5, 0; - 0.5, 0];
%     xrange = 10;
%     time = 81;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 8, time); 
% 
%     init_min1 = [0, 0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0];
%     init_min2 = [0, 0.4, 0.4, 0.5, 0.6, 0, 0, 0.2, 0, 0, 0];
%     init_min = [init_min1; init_min2];   
%     
%     init_max1 = [0, 0.2, 0.2, 0.4, 0.5, 0, 0, 0, 0, 0, 0];
%     init_max2 = [0, 0.4, 0.4, 0.5, 0.6, 0, 0, 0.2, 0, 0, 0];
%     init_max = [init_max1; init_max2]; 
%     
%     tic
%     reach_3dlinhyposys(A, deltax, deltat, xlist, tlist, init_min, init_max, time);
%     toc

%----------------------------2dbox plot------------------------------------

%     deltax = 1; 
%     xlist = linspace(0, 10, 11);
%     init_min = [0, 0.2, 0.2, 0, 0, 0, 0.2, 1, 1, 0, 0];    
%     init_max = [0, 0.5, 0.5, 0, 0.7, 0, 0.8, 1.5, 1.5, 0, 0];
% 
%     plot_2dbox(init_min, init_max, deltax, xlist);

%----------------------------3dbox plot------------------------------------

%     deltat = 0.1;
%     deltax = 1;
%     xrange = 10;
%     time = 51;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 5, time);    
%     min_matrix = [0, 0.2, 0.2, 0.4, 0.5, 0, 0, 1, 1, 0.9, 0;
%                 0, 0.2, 0.2, 0.4, 0.5, 0, 0, 1, 1, 0.9, 0;
%                 0, 0.2, 0.2, 0.4, 0.5, 0, 0, 1, 1, 0.9, 0;
%                 0, 0.2, 0.2, 0.4, 0.5, 0, 0, 1, 1, 0.9, 0;
%                 0, 0.2, 0.2, 0.4, 0.5, 0, 0, 1, 1, 0.9, 0];    
%             
%     max_matrix = [0, 0.4, 0.4, 0.5, 0.6, 0.1, 0.1, 1.2, 1.2, 1, 0;
%                 0, 0.4, 0.4, 0.5, 0.6, 0.1, 0.1, 1.2, 1.2, 1, 0;
%                 0, 0.4, 0.4, 0.5, 0.6, 0.1, 0.1, 1.2, 1.2, 1, 0;
%                 0, 0.4, 0.4, 0.5, 0.6, 0.1, 0.1, 1.2, 1.2, 1, 0;
%                 0, 0.4, 0.4, 0.5, 0.6, 0.1, 0.1, 1.2, 1.2, 1, 0];
% 
%     plot_3dbox(min_matrix, max_matrix, deltax, deltat, xlist, tlist);

%-------------------nonlin eq lax-wdf method 2d and 3d reachable set-------------------------------------------------
% 
    deltat = 0.5;       % dont change this, use large time step
    deltax = 1;    % delta x = 1, 2, 4
    xrange = 100;
    time = 5/deltat + 1;    % simluate from 0 to 5 seconds, too many optimizations will lead to over approximation
    xlist = linspace(0, xrange, xrange/deltax + 1);
    tlist = linspace(0, 5, time);    
    bdcnd = 'Dirichlet';
    
    m = xrange/deltax + 1;
    init_min = zeros(1, m);
    init_min(1, 1: m * 0.3) = 0.5;
    init_max = zeros(1, m);
    init_max(1, 1: m * 0.3) = 1;   
    
    figure;
%     init_min = [0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0]; % deltax = 1
%     init_max = [0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0];
    
%     init_min = [0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0]; 
%     init_max = [0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0];   
    [sol_min, sol_max] = reach_nonlnhypo(deltat, deltax, init_min', init_max', time, xlist, tlist, bdcnd);
    plot_3dbox(sol_min', sol_max', deltax, deltat, xlist, tlist);
    axis([-2 80 0 6])
%-------------------nonlin eq lax-wdf method min_function-------------------------------------------------    
%     deltat = 0.1;
%     deltax = 1;    
%     xrange = 10;
%     time = 31;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 3, time);    
%     u_min = [0, 1.8, 0.9, 0.9, 0.9, 0, 0, 1, 1, 1.4, 0];  
%     u_max = [0, 2.2, 1.2, 1.4, 1.5, 0, 0, 1.2, 1.2, 1.9, 0];     
%     [list, umin] = min_nonlnhypo_local(u_min, u_max, deltat, deltax, xlist)

%-------------------nonlin eq lax-wdf method max_function-------------------------------------------------    
%     deltat = 0.1;
%     deltax = 1;    
%     xrange = 10;
%     time = 31;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 3, time);    
%     u_min = [0, 1.8, 0.9, 0.9, 0.9, 0, 0, 1, 1, 1.4, 0];  
%     u_max = [0, 2.2, 1.2, 1.4, 1.5, 0, 0, 1.2, 1.2, 1.9, 0];   
%     [list, umax] = max_nonlnhypo_local(u_min, u_max, deltat, deltax, xlist)

%----------------------------nonlin eq lax-wdf method method 2d and 3dbox plot------------------------------------   
%  




% 
%     deltat = 0.1;
%     deltax = 1;    
%     xrange = 10;
%     time = 5/deltat + 1;
%     numofmesh = xrange/deltax + 1;
%     xlist = linspace(0, xrange, numofmesh);
%     tlist = linspace(0, 5, time);
%     bdcnd = 'Dirichlet';
%     
%     
%     init_min = zeros(1, numofmesh);
%     init_min(1 , 1: numofmesh/4) = 0.5;    
%       
%     init_max = zeros(1, numofmesh);
%     init_max(1, 1: numofmesh/4 ) = 1;   
%     
%     
%     reach_3dnonlinhypo(deltat, deltax, xlist, tlist, init_min, init_max, time, bdcnd);
%     

%-------------------------solve euler equations using lax-wendroff method-----------------------------------
  
%     deltat = 0.1;
%     deltax = 1;    
%     xrange = 10;
%     gamma = 1.4;
%     time = 51;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 5, time);  
% 
%     init1 = [1, 1, 1, 1, 1, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7];
%     init2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
%     init3 = [2.5, 2.5, 2.5, 2.5, 2.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5];
%     
%     init = [init1; init2; init3];
%     
%     
%     solve_sys_nonlnhypo(gamma, deltat, deltax, init, time, xlist, tlist);


%-------------------------1d euler equations reachable set-----------------------------------
% 
%     deltat = 0.5;
%     deltax = 5;    
%     xrange = 50;
%     gamma = 1.4;
%     time = 5/deltat + 1; 
%     numofmesh = xrange/deltax + 1;
% 
%     xlist = linspace(0, xrange, numofmesh);
%     tlist = linspace(0, 5, time);  
% 
% 
%     init_min1 = zeros(1, numofmesh);
%     init_min1(1 , 1: numofmesh/4) = 0.7;    
%       
%     init_max1 = zeros(1, numofmesh);
%     init_max1(1, 1: numofmesh/4 ) = 1;
% 
%     init_min2 = zeros(1, numofmesh);
%     init_min2(1 , 1: numofmesh/4) = 0;    
%       
%     init_max2 = zeros(1, numofmesh);
%     init_max2(1, 1: numofmesh/4 ) = 0.5;
% 
%     init_min3 = zeros(1, numofmesh);
%     init_min3(1 , 1: numofmesh/4) = 1.5;    
%       
%     init_max3 = zeros(1, numofmesh);
%     init_max3(1, 1: numofmesh/4 ) = 2.5;
% 
% 
% %     init_min1 = [1, 1, 1, 1, 1, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7];
% %     init_min2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
% %     init_min3 = [2.5, 2.5, 2.5, 2.5, 2.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5];
% %     
%     init_min = [init_min1; init_min2; init_min3];
% %     
% %     init_max1 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3];
% %     init_max2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
% %     init_max3 = [1.5, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
% %     
%     init_max = [init_max1; init_max2; init_max3];
% %     
%     reach_sys_nonlnhypo(gamma, deltat, deltax, init_min, init_max, time, xlist, tlist);
   


%-------------------euler eqs min, max function-------------------------------------------------    
%     deltat = 0.1;
%     deltax = 1;    
%     xrange = 10;
%     xlist = linspace(0, xrange, 11);
%     gamma = 1.4;
%     
%     init_min1 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3];
%     init_min2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
%     init_min3 = [1.5, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];  
%     
%     init_min = [init_min1; init_min2; init_min3];   
%     
%     init_max1 = [1, 1, 1, 1, 1, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7];
%     init_max2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
%     init_max3 = [2.5, 2.5, 2.5, 2.5, 2.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5];
%     
%     init_max = [init_max1; init_max2; init_max3]; 
%     
%     list = min_sys_nonlnhypo(init_min, init_max, deltat, deltax, xlist, gamma);
%     disp('min is');
%     disp(list);
%     
%     list = max_sys_nonlnhypo(init_min, init_max, deltat, deltax, xlist, gamma);
%     disp('max is');    
%     disp(list);

%-------------------euler eqs reachable sets-------------------------------------------------    
%     deltat = 1;
%     deltax = 1;    
%     xrange = 10;
%     time = 5/deltat + 1;
%     xlist = linspace(0, xrange, xrange/deltax + 1);
%     gamma = 1.4;
%     tlist = linspace(0, 5, 5/deltat + 1);  
%     numofmesh = xrange/deltax + 1;
%     
%     
% %     init_min1 = zeros(1, numofmesh);
% %     init_min1(1 , 1: int8(numofmesh/2)) = 0.7;    
% %       
% %     init_max1 = zeros(1, numofmesh);
% %     init_max1(1, 1: int8(numofmesh/2) ) = 1;
% % 
% %     init_min2 = zeros(1, numofmesh);
% %     init_min2(1 , 1: int8(numofmesh/2)) = 0;    
% %       
% %     init_max2 = zeros(1, numofmesh);
% %     init_max2(1, 1: int8(numofmesh/2) ) = 0.5;
% % 
% %     init_min3 = zeros(1, numofmesh);
% %     init_min3(1 , 1: int8(numofmesh/2)) = 1.5;    
% %       
% %     init_max3 = zeros(1, numofmesh);
% %     init_max3(1, 1: int8(numofmesh/2) ) = 2.5;
% 
%     init_min1 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3];
%     init_min2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
%     init_min3 = [1.5, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];  
% 
% %     init_min1 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.3 ];
% %     init_min2 = [0, 0, 0, 0, 0, 0, ];
% %     init_min3 = [1.5, 1.5, 1.5, 1.5, 1.5, 0.5];      
%     init_min = [init_min1; init_min2; init_min3];   
%     
%     init_max1 = [1, 1, 1, 1, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7];
%     init_max2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
%     init_max3 = [2.5, 2.5, 2.5, 2.5, 2.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5];
% %     
% %     init_max1 = [1, 1, 1, 1, 0.7, 0.7];
% %     init_max2 = [0, 0, 0, 0, 0, 0];
% %     init_max3 = [2.5, 2.5, 2.5, 2.5, 2.5, 1.5];    
%     init_max = [init_max1; init_max2; init_max3]; 
%     disp(init_min);
%     disp(init_max);
%     reach_sys_nonlnhypo(gamma, deltat, deltax, init_min, init_max, time, xlist, tlist);


%-------------------------solve boundary control problem-------------------------------------------  
%     deltat = 0.1;
%     deltax = 1;    
%     xrange = 10;
%     alpha = -0.5;
%     time = 101;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 10, time);     
%     init = [0, 0, 0, 0, 1.5, 1.5, 1.5, 0, 0, 0, 0];
%     lambda = 0.5;
%     sol = solve_boundary_control(alpha, deltat, deltax, init, time, xlist, tlist, lambda);
    
%-------------------------boundary control problem reachable sets-------------------------------------------      
%     deltat = 0.1;
%     deltax = 1;
%     alpha = -0.5;
%     xrange = 10;
%     time = 81;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 8, time);
%     lambda_1 = 0.4;
%     lambda_2 = 0.6;
%     init_min = [0, 0, 0, 0, 1.0, 1.0, 1.0, 0, 0, 0, 0];    
%     init_max = [0, 0, 0, 0, 1.3, 1.2, 1.2, 0, 0, 0, 0];
%     tic
%     reach_3d_boundary_control(alpha, deltat, deltax, init_min, init_max, time, xlist, tlist, lambda_1, lambda_2);
%     toc


%-------------------------solve heat eq using tradiagonal method------------------------------------------- 

%     deltat = 0.05;
%     deltax = 0.01;    
%     xrange = 1;
%     trange = 5;    
%     time = 101;    
%     xlist = linspace(0, xrange, 101);
%     tlist = linspace(0, trange, time);     
%     init = zeros(101,1);
%     init(20:60) = 1;
% 
%     sol = solve_implicit_tradiagonal(deltat, deltax, init, time, xlist, tlist);
    
%-------------------------reach 1d heat eq using gminres and tradiagonal method------------------------------------------- 

%     deltat = 0.1;
%     deltax = 0.01/8;    
%     
% %     xrange = 10;
%     xrange = 9;
%     
%     trange = 10;
%     time = trange/deltat + 1;
%     space = xrange/deltax + 1;
%     
%     xlist = linspace(0, xrange, xrange/deltax + 1);
%     tlist = linspace(0, trange, trange/deltat + 1);     
%     init_min = zeros(space,1);
%         
%     init_min(20:60) = 0.5;
%     lambda_min = 0.1;
%     lambda_max = 0.2;
%     
%     init_max = zeros(space,1);
%     
%     init_max(20:60) = 1;  
%     
%     reach_implicit_1d(deltat, deltax, init_min, init_max, lambda_min, lambda_max, time, xlist, tlist);    

%-------------------------solve 2d heat eq using CN and gminres ------------------------------------------- 

%     deltat = 0.1;
%     deltax = 0.25;  
%     deltay = 0.25;
%     xrange = 10;
%     yrange = 10;    
%     trange = 10;
%     alpha = 0.5;
%     nummesh = xrange/deltax + 1;    
%     time = trange/deltat;    
%     xlist = linspace(0, xrange, nummesh);
%     ylist = linspace(0, yrange, nummesh);    
%     tlist = linspace(0, trange, time);   
%     init = zeros(nummesh*nummesh, 1);
%         
%     for i = 1 : nummesh
%         for j = 1 : nummesh
%             
%             if(i >= nummesh/4 && i <= nummesh/4*3 && ...
%                 j >= nummesh/4 && j <= nummesh/4*3)
%             
%                 init((i - 1) * nummesh + j) = 1;
%             else
%                 init((i - 1) * nummesh + j) = 0;
%             end
%         end
%     end
%     
%     lambda = 0.2;    
%     sol = solve_CN(deltat, deltax, deltay, lambda, alpha, init, time, xlist, ylist, tlist);   
    
%-------------------------reach 2d heat eq using CN and gminres ------------------------------------------- 

%     deltat = 0.1;
%     deltax = 0.2;  
%     deltay = 0.2;
%     xrange = 10;
%     yrange = 10;    
%     trange = 10;
%     alpha = 0.5;
%     time = trange/deltat; 
%     nummesh = xrange/deltax + 1;
%     xlist = linspace(0, xrange, nummesh);
%     ylist = linspace(0, yrange, nummesh);    
%     tlist = linspace(0, trange, time);   
%     
%     
%     
%     init_min = zeros(nummesh * nummesh, 1);
%         
%     for i = 1 : nummesh
%         for j = 1 : nummesh
%             
%             if(i >= nummesh/4 && i <= nummesh/4*3 && ...
%                 j >= nummesh/4 && j <= nummesh/4*3)
%             
%                 init_min((i - 1) * nummesh + j) = 1;
%             else
%                 init_min((i - 1) * nummesh + j) = 0;
%             end
%         end
%     end
% 
%     
%     
%     
%     init_max = zeros(nummesh * nummesh, 1);
%         
%     for i = 1 : nummesh
%         for j = 1 : nummesh
%             
%             if(i >= nummesh/4 && i <= nummesh/4*3 && ...
%                 j >= nummesh/4 && j <= nummesh/4*3)
%             
%                 init_max((i - 1) * nummesh + j) = 1.5;
%             else
%                 init_max((i - 1) * nummesh + j) = 0.5;
%             end
%         end
%     end
%     
%     lambda_min = 0.1;  
%     lambda_max = 0.2; 
%     
%     reach_3d_CN_2dheat(deltat, deltax, deltay, lambda_min, lambda_max, alpha, init_min, init_max, time, xlist, ylist, tlist);

%-------------------------solve 2d heat eq using ADI, tridiagonal ------------------------------------------- 

%     deltat = 0.1;
%     deltax = 0.2;  
%     deltay = 0.2;
%     xrange = 10;
%     yrange = 10;    
%     trange = 10;
%     alpha = 0.5;
%     time = trange/deltat; 
%     nummesh = xrange/deltax + 1;
%     xlist = linspace(0, xrange, nummesh);
%     ylist = linspace(0, yrange, nummesh);    
%     tlist = linspace(0, trange, time);   
%     
%     
%     
%     init = zeros(nummesh * nummesh, 1);
%         
%     for i = 1 : nummesh
%         for j = 1 : nummesh
%             
%             if(i >= nummesh/4 && i <= nummesh/4*3 && ...
%                 j >= nummesh/4 && j <= nummesh/4*3)
%             
%                 init((i - 1) * nummesh + j) = 1.5;
%             else
%                 init((i - 1) * nummesh + j) = 0.5;
%             end
%         end
%     end
% 
%     
%     
%     
%     lambda = 0.15;
%     
%     solve_peaceman(alpha, deltat, deltax, deltay, lambda, init, time, xlist, ylist, tlist);

%-------------------------test block tridiagonal -------------------------------------------

%     b = 2*ones(1,9);
%     b(1) = 0;
%     b(2) = 0;
%     b(3) = 0;    
%     a = ones(1,9);
%     c = 3*ones(1,9);
%     c(7) = 0;
%     c(8) = 0;
%     c(9) = 0;     
%     f = 10*rand(9, 1);
%     y = block_tridiag(a, b, c, 3, f);
%     
%     
%     
%     
%     A = diag(ones(1, 9)) + diag(3 * ones(1, 6), 3) + diag(2*ones(1, 6), -3); 
%     disp(A);
%     z = A\f;
%     disp(y);
%     disp(z);    
%     disp(y - z);
%     

%-------------------------reach 2d heat eq using ADI ------------------------------------------- 

%     deltat = 0.1;
%     deltax = 0.2;  
%     deltay = 0.2;
%     xrange = 10;
%     yrange = 10;    
%     trange = 10;
%     alpha = 0.5;
%     time = trange/deltat; 
%     nummesh = xrange/deltax + 1;
%     xlist = linspace(0, xrange, nummesh);
%     ylist = linspace(0, yrange, nummesh);    
%     tlist = linspace(0, trange, time);   
%     
%     
%     
%     init_min = zeros(nummesh * nummesh, 1);
%         
%     for i = 1 : nummesh
%         for j = 1 : nummesh
%             
%             if(i >= nummesh/4 && i <= nummesh/4*3 && ...
%                 j >= nummesh/4 && j <= nummesh/4*3)
%             
%                 init_min((i - 1) * nummesh + j) = 1;
%             else
%                 init_min((i - 1) * nummesh + j) = 0;
%             end
%         end
%     end
% 
%     
%     
%     
%     init_max = zeros(nummesh * nummesh, 1);
%         
%     for i = 1 : nummesh
%         for j = 1 : nummesh
%             
%             if(i >= nummesh/4 && i <= nummesh/4*3 && ...
%                 j >= nummesh/4 && j <= nummesh/4*3)
%             
%                 init_max((i - 1) * nummesh + j) = 1.5;
%             else
%                 init_max((i - 1) * nummesh + j) = 0.5;
%             end
%         end
%     end
%     
%     lambda_min = 0.1;  
%     lambda_max = 0.2; 
%     
%     reach_3d_ADI_2dheat(deltat, deltax, deltay, lambda_min, lambda_max, alpha, init_min, init_max, time, xlist, ylist, tlist);
    
%-------------------------solve 1d heat eq using ASC-N------------------------------------------- 

%     deltat = 0.1;
%     deltax = 0.1;   
%     alpha = 0.5;
%     xrange = 10;
%     trange = 10;    
%     time = trange/deltat + 1;    
%     xlist = linspace(0, xrange, xrange/deltax);
%     tlist = linspace(0, trange, time);     
%     init = zeros(xrange/deltax, 1);
%     lambda = 0.15;
%     init(35:70) = 1;
% 
%     sol = solve_ASCN(deltat, deltax, alpha, init, time, xlist, tlist, lambda);    
    
%-------------------------reach 1d heat eq using ASC-N-------------------------------------------

%     deltat = 0.1;
%     deltax = 0.01/8;   
%     alpha = 0.5;
%     
% %     xrange = 10;
% 
%     xrange = 9; %900 points
%     trange = 10;    
%     time = trange/deltat + 1;    
%     xlist = linspace(deltax, xrange, xrange/deltax);        %100 points in stead of 101
%     tlist = linspace(0, trange, time);  
%     lambda_min = 0.1;
%     lambda_max = 0.15;   
%     
%     init_min = zeros(xrange/deltax, 1);
%     init_max = zeros(xrange/deltax, 1);
%     
%     init_min(35:70) = 0.5;
%     init_max(35:70) = 1;
%     
%     
%     reach_3d_ASCN_1dheat(deltat, deltax, lambda_min, lambda_max, alpha, init_min, init_max, time, xlist, tlist)



end
    
    
    
    
    
    
    
    
    
    
    