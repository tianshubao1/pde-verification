function [sol_min, sol_max] = reach_linhypo(alpha, deltat, deltax, init_min, init_max, time, xlist, tlist, bdcnd)


    sol_min = solve_hypo(alpha, deltat, deltax, init_min, time, xlist, tlist, bdcnd);
    sol_max = solve_hypo(alpha, deltat, deltax, init_max, time, xlist, tlist, bdcnd);
    true_value = [0, 0, 0, 0, 0, 0, 0.6, 0.6, 0.6, 0.6, 0]; 
    
    
%--------annotation for linear eq at t = 4s------------------%  
    plot_2dboxbald(sol_min(:,4/deltat + 1), sol_max(:,4/deltat + 1), deltax, xlist, 0.5);    
    
    
%--------annotation for linear eq------------------%    
%     x = [0.5,0.5];
%     y = [0.46,0.39];
%     a = annotation('textarrow',x,y,'String','Single traces');
% --------------plot safe region for lin eq--------------------------------% 
    safe_line = 1.0*ones(15,1);
    hold on;
    ylim([0 1.1])
    t = text(4, 1.05,'Unsafe region');
    t.Color = [1 0 0];
    plot(linspace(-2, 12, 15), safe_line, '--r'); 
    
% --------------plot unsafe region ------------------------%     
    plot_fixlocbox(sol_min(5, :), sol_max(5, :), deltat, tlist);
    hold on;
    rectangle('Position',[4   0.75  0.5  0.2], 'FaceColor',[1 0 0], 'EdgeColor',[1 0 0]); 
    t = text(4, 0.73,'Unsafe region');
    hold;
end