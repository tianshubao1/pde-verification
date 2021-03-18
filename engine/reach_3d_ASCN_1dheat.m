function reach_3d_ASCN_1dheat(deltat, deltax, lambda_min, lambda_max, alpha, init_min, init_max, time, xlist, tlist)
    

    sol_min = solve_ASCN(deltat, deltax, alpha, init_min, time, xlist, tlist, lambda_min);               
             
    sol_max = solve_ASCN(deltat, deltax, alpha, init_max, time, xlist, tlist, lambda_max);

%     sol_min = solve_ASCN5(deltat, deltax, alpha, init_min, time, xlist, tlist, lambda_min);               
%              
%     sol_max = solve_ASCN5(deltat, deltax, alpha, init_max, time, xlist, tlist, lambda_max);
    plot_2dbox(sol_min(:,4/deltat + 1), sol_max(:,4/deltat + 1), deltax, xlist);   
    axis([-1 11 0 inf])
    
    plot_3dbox(sol_min', sol_max', deltax, deltat, xlist, tlist)
end