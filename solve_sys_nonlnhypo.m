function sol = solve_sys_nonlnhypo(gamma, deltat, deltax, init, time, xlist, tlist)         %solve 1D euler equations

u = init;
len = size(u);
m = len(2);
% disp(m);
sol = zeros(3, length(xlist), length(tlist));
sol(:, :, 1) = init; %x is horizontal, 3 is vertical, t is z
u_new = zeros(3, m);

    function result = fun(rho, mmt, engy, gamma) %flux function
        f1 = mmt;
        f2 = 1/2*(3 - gamma)*mmt^2/rho + (gamma - 1)*engy;
        f3 = gamma*mmt*engy/rho - 1/2*(gamma - 1)*mmt^2/(rho^2);
        result = [f1; f2; f3];
    end

for i = 1 : time - 1
    
    for j = 1 : m   %number of mesh points, richtmyer method
      
        if j == 1   %left boundary
            
            u_1 = (u(:, j + 1) + u(:, j))/2 - deltat/(2 * deltax) * ...
                (fun(u(1, j + 1), u(2, j + 1), u(3, j + 1), gamma) - fun(u(1, j), u(2, j), u(3, j), gamma));  %u_+1/2
            
            u_2 = (u(:, j) + [1; 0; 2])/2 - deltat/(2 * deltax) * ...
                (fun(u(1, j), u(2, j), u(3, j), gamma) - fun(1, 0, 2, gamma));  %u_-1/2, left bndy equals 0.7

            u_new(:, j) = u(:, j) - deltat/deltax * ...
                (fun(u_1(1), u_1(2), u_1(3), gamma) - fun(u_2(1), u_2(2), u_2(3), gamma));     %f(u_i+1/2)-f(u_i-1/2)

            u_new(:, j) = [0.5; 0; 2];
            
        elseif j == m   %right boundary
            
            u_1 = ([0.7; 0; 1.0] + u(:, j))/2 - deltat/(2 * deltax) * ...
                (fun(0.7, 0, 1.0, gamma) - fun(u(1, j), u(2, j), u(3, j), gamma));  %u_+1/2
            
            u_2 = (u(:, j) + u(:, j - 1))/2 - deltat/(2 * deltax) * ...
                (fun(u(1, j), u(2, j), u(3, j), gamma) - fun(u(1, j - 1), u(2, j - 1), u(3, j - 1), gamma));  %u_-1/2
            
            u_new(:, j) = u(:, j) - deltat/deltax * ...
                (fun(u_1(1), u_1(2), u_1(3), gamma) - fun(u_2(1), u_2(2), u_2(3), gamma));     %f(u_i+1/2)-f(u_i-1/2)

            u_new(:, j) = [0.5; 0; 1.0];  
            
        else    %in the middle
                      
            u_1 = (u(:, j + 1) + u(:, j))/2 - deltat/(2 * deltax) * ...
                (fun(u(1, j + 1), u(2, j + 1), u(3, j + 1), gamma) - fun(u(1, j), u(2, j), u(3, j), gamma));  %u_+1/2
            u_2 = (u(:, j) + u(:, j - 1))/2 - deltat/(2 * deltax) * ...
                (fun(u(1, j), u(2, j), u(3, j), gamma) - fun(u(1, j - 1), u(2, j - 1), u(3, j - 1), gamma));  %u_-1/2
            u_new(:, j) = u(:, j) - deltat/deltax * (fun(u_1(1), u_1(2), u_1(3), gamma) - fun(u_2(1), u_2(2), u_2(3), gamma));     %f(u_i+1/2)-f(u_i-1/2)
        end

    end
%     disp(u);
    sol(:, :, i + 1) = u_new;
    u = u_new;
end

sol1 = squeeze(sol(1,:,:))';
figure;
surf(xlist,tlist,sol1) 
title('Rho Numerical solution')
xlabel('Distance x')
ylabel('Time t')

end