function reach_3dlinhypo(alpha, deltat, deltax, xlist, tlist, init_min, init_max, time, bdcnd)

    sol_min = solve_hypo(alpha, deltat, deltax, init_min, time, xlist, tlist, bdcnd);
    sol_max = solve_hypo(alpha, deltat, deltax, init_max, time, xlist, tlist, bdcnd);
    
    
    plot_3dbox(sol_min', sol_max', deltax, deltat, xlist, tlist);


end