function switchedtrafficflow


%deltat and h must be equal
h = 0.1;
deltat = 0.1;
m = 10/h + 1; %number of mesh points, 51 mesh points including 2 edge points
t = 5/deltat + 1; %number of time steps
v1 = 3;
v2 = -1;
rho_c = 1;

u = zeros(m - 2, t);
u0 = zeros(m - 2, 1);
u0(10) = 1.2;
u0(11) = 1.2;
u0(12) = 1.2;
u0(20) = 0.8;
u0(21) = 0.8;
u0(35) = 1.2;
u0(36) = 1.2;
u0(37) = 1.2;
u0(39) = 0.8;
u0(40) = 0.8;
u(:,1) = u0;


for i = 1 : t - 1  %time point
    for j = 1 : m - 2   %mesh point    
        if u(j, i) >= rho_c && j + v2 >= 1  %v2 = -1, backward wave
            u(j + v2, i + 1) = u(j, i) + u(j + v2, i + 1);    
        end
        
        if u(j, i) < rho_c && j + v1 <= m - 2
            u(j + v1, i + 1) = u(j, i) + u(j + v1,i + 1);   %v1 = 4, forward wave
        end
        
    end
    %size(u)
    %disp(u(:,1));
end


xlist = linspace(0,10,m);
tlist = linspace(0,5,t);
u = [zeros(1, t); u; zeros(1, t)];
u = u';


size(xlist);
size(tlist);
size(u);



surf(xlist,tlist,u) 
title('Numerical solution computed with 100 mesh points.')
xlabel('Distance x')
ylabel('Time t')

figure
plot(xlist,u(11,:))
title('Solution at t = 1')
xlabel('Distance x')
ylabel('u')


figure
plot(xlist,u(31,:))
title('Solution at t = 3')
xlabel('Distance x')
ylabel('u')

end