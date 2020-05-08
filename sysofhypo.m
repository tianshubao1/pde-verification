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
%     xrange = 10;
%     time = 51;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 5, time);
%     bdcnd = 'Dirichlet';
%     init_min = [0, 0.2, 0.2, 0.4, 0.5, 0, 0, 1, 1, 0.9, 0];    
%     init_max = [0, 0.4, 0.4, 0.5, 0.6, 0.1, 0.1, 1.2, 1.2, 1, 0];
%     [sol_min, sol_max] = reach_linhypo(alpha, deltat, deltax, init_min, init_max, time, xlist, tlist, bdcnd);

% %----------------------------one side method 3dbox reachable set ------------------------------------
% 
%   this is correct 
%     deltat = 0.1;
%     deltax = 1;
%     alpha = 0.5;
%     xrange = 10;
%     time = 81;    
%     xlist = linspace(0, xrange, 41);
%     tlist = linspace(0, 8, time);
%     bdcnd = 'Robin';
%     init_min = [0, 0.2, 0.2, 0.4, 0.5, 0, 0, 1, 1, 0.9];    
%     init_max = [0, 0.4, 0.4, 0.5, 0.6, 0.1, 0.1, 1.2, 1.2, 1];
%     tic
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

%-------------------nonlin eq lax-wdf method reachable set-------------------------------------------------
%     deltat = 0.1;
%     deltax = 1;    
%     xrange = 10;
%     time = 61;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 6, time);    
%     bdcnd = 'Dirichlet';
%     init_min = [0, 0.2, 0.2, 0.4, 0.5, 0, 0, 1, 1, 0.9, 0];  
%     init_max = [0, 0.2, 1.2, 0.4, 1.5, 0, 0, 1.2, 1.2, 1.9, 0]; 
%     [sol_min, sol_max] = reach_nonlnhypo(deltat, deltax, init_min, init_max, time, xlist, tlist, bdcnd);
    
%-------------------nonlin eq lax-wdf method min_function-------------------------------------------------    
%     deltat = 0.1;
%     deltax = 1;    
%     xrange = 10;
%     time = 31;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 3, time);    
%     u_min = [0, 1.8, 0.9, 0.9, 0.9, 0, 0, 1, 1, 1.4, 0];  
%     u_max = [0, 2.2, 1.2, 1.4, 1.5, 0, 0, 1.2, 1.2, 1.9, 0];     
%     [list, umin] = min_nonlnhypo(u_min, u_max, deltat, deltax, xlist)

%-------------------nonlin eq lax-wdf method max_function-------------------------------------------------    
%     deltat = 0.1;
%     deltax = 1;    
%     xrange = 10;
%     time = 31;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 3, time);    
%     u_min = [0, 1.8, 0.9, 0.9, 0.9, 0, 0, 1, 1, 1.4, 0];  
%     u_max = [0, 2.2, 1.2, 1.4, 1.5, 0, 0, 1.2, 1.2, 1.9, 0];   
%     [list, umax] = max_nonlnhypo(u_min, u_max, deltat, deltax, xlist)

%----------------------------nonlin eq lax-wdf method method 2d and 3dbox plot------------------------------------   
 
%     deltat = 0.1;
%     deltax = 1;    
%     xrange = 10;
%     time = 51;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 5, time);
%     bdcnd = 'Neumann';
%     init_min = [0, 0.8, 0.8, 0.4, 0, 0, 0, 0, 0, 0, 0];    
%     init_max = [0, 1, 1, 0.5, 0, 0, 0, 0, 0, 0, 0];
%     reach_3dnonlinhypo(deltat, deltax, xlist, tlist, init_min, init_max, time, bdcnd);
    

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

%     ???
%     deltat = 0.1;
%     deltax = 1;    
%     xrange = 10;
%     gamma = 1.4;
%     time = 51;    
%     xlist = linspace(0, xrange, 11);
%     tlist = linspace(0, 5, time);  
% 
%     init_min1 = [1, 1, 1, 1, 1, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7];
%     init_min2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
%     init_min3 = [2.5, 2.5, 2.5, 2.5, 2.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5];
%     
%     init_min = [init_min1; init_min2; init_min3];
%     
%     init_max1 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3];
%     init_max2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
%     init_max3 = [1.5, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];
%     
%     init_max = [init_max1; init_max2; init_max3];
%     
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
%     deltat = 0.1;
%     deltax = 1;    
%     xrange = 10;
%     time = 61;
%     xlist = linspace(0, xrange, 11);
%     gamma = 1.4;
%     tlist = linspace(0, 6, time);  
%     
%     
%     init_min1 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3];
%     init_min2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
%     init_min3 = [1.5, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];  
%     
%     init_min = [init_min1; init_min2; init_min3];   
%     
%     init_max1 = [1, 1, 1, 1, 1, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7];
%     init_max2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
%     init_max3 = [2.5, 2.5, 2.5, 2.5, 2.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5];
%     
%     init_max = [init_max1; init_max2; init_max3]; 
%     tic
%     reach_sys_nonlnhypo(gamma, deltat, deltax, init_min, init_max, time, xlist, tlist);
%     toc

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


%-------------------------solve heat eq using gminres subspace------------------------------------------- 

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
%     sol = solve_implicit(deltat, deltax, init, time, xlist, tlist);
    
%-------------------------reach heat eq using gminres subspace------------------------------------------- 

    deltat = 0.05;
    deltax = 0.01;    
    xrange = 1;
    trange = 5;
    time = 101;    
    xlist = linspace(0, xrange, 101);
    tlist = linspace(0, trange, time);     
    init_min = zeros(101,1);
        
    init_min(20:60) = 0.8;
    lambda_min = 0.4;
    lambda_max = 1.2;
    
    init_max = zeros(101,1);
    
    init_max(20:60) = 1;  
    
    sol = reach_implicit_krylov(deltat, deltax, init_min, init_max, lambda_min, lambda_max, time, xlist, tlist);    
end
    
    
    
    
    
    
    
    
    
    
    