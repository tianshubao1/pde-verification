function reach_3dlinhyposys(A, deltax, deltat, xlist, tlist, init_min, init_max, time)
    
    sol_min = solve_jd_hypo(A, deltat, deltax, init_min, time, xlist, tlist);    
    sol_max = solve_jd_hypo(A, deltat, deltax, init_max, time, xlist, tlist);   
    
    sol_min = squeeze(sol_min(2, :, :))';
    sol_max = squeeze(sol_max(2, :, :))';
    plot_3dbox(sol_min, sol_max, deltax, deltat, xlist, tlist);

end