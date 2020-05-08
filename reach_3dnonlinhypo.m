function reach_3dnonlinhypo(deltat, deltax, xlist, tlist, init_min, init_max, time, bdcnd)

    [sol_min, sol_max] = reach_nonlnhypo(deltat, deltax, init_min, init_max, time, xlist, tlist, bdcnd);
    
    plot_3dbox(sol_min', sol_max', deltax, deltat, xlist, tlist);


end