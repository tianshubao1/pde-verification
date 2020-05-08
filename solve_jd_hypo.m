function sol = solve_jd_hypo(A, deltat, deltax, init, time, xlist, tlist)
%-------------------- Lax-Friedrich scheme-------------------------%
    u = init;
    len = size(u);
    n = len(1);
    m = len(2);
    sol = zeros(n, length(xlist), length(tlist));
    sol(:, :, 1) = init; %x is horizontal, u is vertical, t is z
    
    for i = 1 : time - 1
            u_new = 1/2*([u(:, 2:m), zeros(n, 1)] + [zeros(n, 1), u(:, 1:(m - 1))]) - A * deltat/(2*deltax) * ([u(:, 2:m), zeros(n, 1)] - [zeros(n, 1), u(:, 1:(m - 1))]);
            sol(:, :, i) = u_new;
            u = u_new;
    end
     
    figure;
    sol2 = squeeze(sol(2,:,:))';
    size(xlist);
    size(tlist);
    size(sol2);
    
    surf(xlist,tlist,sol2) 
    title('Numerical solution for density')
    xlabel('Distance x')
    ylabel('Time t')  
    
end
    