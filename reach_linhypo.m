function [sol_min, sol_max] = reach_linhypo(alpha, deltat, deltax, init_min, init_max, time, xlist, tlist, bdcnd)


    sol_min = solve_hypo(alpha, deltat, deltax, init_min, time, xlist, tlist, bdcnd);
    sol_max = solve_hypo(alpha, deltat, deltax, init_max, time, xlist, tlist, bdcnd);
    sol_sample1 = solve_hypo(alpha, deltat, deltax, (init_min + init_max)/2, time, xlist, tlist, bdcnd);
    sol_sample2 = solve_hypo(alpha, deltat, deltax, (init_max - init_min)/3 + init_min, time, xlist, tlist, bdcnd);
    sol_sample3 = solve_hypo(alpha, deltat, deltax, (init_max - init_min)/3 * 2 + init_min, time, xlist, tlist, bdcnd); 
    

    plot_2dbox(sol_min(:,40), sol_max(:,40), sol_sample1(:,40), sol_sample2(:,40), sol_sample3(:,40), deltax, xlist);     
%--------annotation for linear eq------------------%    
    x = [0.5,0.5];
    y = [0.46,0.39];
    a = annotation('textarrow',x,y,'String','Single traces');
% --------------plot safe region for lin eq--------------------------------% 
    safe_line = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9];
    hold on;
    ylim([0 1.1])
    t = text(4, 1.0,'Unsafe region');
    t.Color = [1 0 0];
    plot(linspace(-2, 12, 15), safe_line, '--r'); 
    
% --------------plot fixed location reachable set ------------------------%     
    plot_fixlocbox(sol_min(5, :), sol_max(5, :), deltat, tlist);
    hold on;
    rectangle('Position',[4   0.45  0.5  0.15], 'FaceColor',[1 0 0], 'EdgeColor',[1 0 0]); 
    t = text(4, 0.43,'Unsafe region');
    hold;
end