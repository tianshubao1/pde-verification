function [sol_min, sol_max] = reach_linhyposys(A, deltat, deltax, init_min, init_max, time, xlist, tlist)

    sol_min = solve_jd_hypo(A, deltat, deltax, init_min, time, xlist, tlist);    
    sol_max = solve_jd_hypo(A, deltat, deltax, init_max, time, xlist, tlist);

%     plot_2dboxbald(sol_min(1, :,41), sol_max(1, :,41), deltax, xlist);     %variable 1, at t = 4s
    plot_2dboxbald(sol_min(2, :,41), sol_max(2, :,4/deltat + 1), deltax, xlist);       %variable 2, at t = 4s
    
    plot_fixlocbox(sol_min(2, 5, :), sol_max(2, 5, :), deltat, tlist);    
end