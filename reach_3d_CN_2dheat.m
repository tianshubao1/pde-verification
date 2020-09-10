function reach_3d_CN_2dheat(deltat, deltax, deltay, lambda_min, lambda_max, alpha, init_min, init_max, time, xlist, ylist, tlist)
    tic
    sol_min = solve_CN(deltat, deltax, deltay, lambda_min, alpha, init_min, time, xlist, ylist, tlist, "min");

    sol_max = solve_CN(deltat, deltax, deltay, lambda_max, alpha, init_max, time, xlist, ylist, tlist, "max");
    toc
    plot_3dbox_xy(sol_min, sol_max, deltax, deltay, xlist, ylist);
end