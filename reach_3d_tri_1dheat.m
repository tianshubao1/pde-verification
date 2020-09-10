function [sol1,sol2] = reach_implicit_1d(deltat, deltax, init_min, init_max, lambda_min, lambda_max, time, xlist, tlist)

    solve_implicit_tradiagonal(deltat, deltax, init, time, xlist, tlist)
    


    plot_3dbox(sol_min', sol_max', deltax, deltat, xlist, tlist);
end