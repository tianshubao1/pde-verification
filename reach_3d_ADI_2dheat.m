function reach_3d_ADI_2dheat(deltat, deltax, deltay, lambda_min, lambda_max, alpha, init_min, init_max, time, xlist, ylist, tlist)
    
    tic
    sol_min = solve_peaceman(alpha, deltat, deltax, deltay, lambda_min, init_min, time, xlist, ylist, tlist);

    sol_max = solve_peaceman(alpha, deltat, deltax, deltay, lambda_max, init_max, time, xlist, ylist, tlist);
    toc
    
    plot_3dbox_xy(sol_min, sol_max, deltax, deltay, xlist, ylist);
end