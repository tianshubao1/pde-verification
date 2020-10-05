function plot_3dbox(min_matrix, max_matrix, deltax, deltat, xlist, tlist)
    m = size(xlist);
    m = m(2);
    n = size(tlist);
    n = n(2);
  
    figure;
    title('3-Dimension Reachable Sets')
    xlabel('Distance x')
    ylabel('Time t')
    zlabel('u(x, t)')
                
    for j = 1 : n     %t direction
        for i = 1 : m   %x direction
            if j == 1 && i == 1 %left below corner
                
                v1 = (0 + min_matrix(j, i))/2;
                v2 = (min_matrix(j + 1, i) + min_matrix(j, i))/2;
                v3 = (min_matrix(j, i) + 0)/2;
                v4 = (min_matrix(j, i) + min_matrix(j, i + 1))/2;
                
                v6 = (0 + max_matrix(j, i))/2;
                v7 = (max_matrix(j + 1, i) + max_matrix(j, i))/2;
                v8 = (max_matrix(j, i) + 0)/2;
                v9 = (max_matrix(j, i) + max_matrix(j, i + 1))/2;
                
            elseif j == 1 && i == m     %right below corner
              
                v1 = (0 + min_matrix(j, i))/2;
                v2 = (min_matrix(j + 1, i) + min_matrix(j, i))/2;
                v3 = (min_matrix(j, i) + min_matrix(j, i - 1))/2;
                v4 = (min_matrix(j, i) + 0)/2;
                
                v6 = (0 + max_matrix(j, i))/2;
                v7 = (max_matrix(j + 1, i) + max_matrix(j, i))/2;
                v8 = (max_matrix(j, i) + max_matrix(j, i - 1))/2;
                v9 = (max_matrix(j, i) + 0)/2;   
                
            elseif j == n && i == 1     %left above corner

                v1 = (min_matrix(j - 1, i) + min_matrix(j, i))/2;
                v2 = (0 + min_matrix(j, i))/2;
                v3 = (min_matrix(j, i) + 0)/2;
                v4 = (min_matrix(j, i) + min_matrix(j, i + 1))/2;
                
                v6 = (max_matrix(j - 1, i) + max_matrix(j, i))/2;
                v7 = (0 + max_matrix(j, i))/2;
                v8 = (max_matrix(j, i) + 0)/2;
                v9 = (max_matrix(j, i) + max_matrix(j, i + 1))/2; 
                
            elseif j == n && i == m
               
                v1 = (min_matrix(j - 1, i) + min_matrix(j, i))/2;
                v2 = (0 + min_matrix(j, i))/2;
                v3 = (min_matrix(j, i) + min_matrix(j, i - 1))/2;
                v4 = (min_matrix(j, i) + 0)/2;
                
                v6 = (max_matrix(j - 1, i) + max_matrix(j, i))/2;
                v7 = (0 + max_matrix(j, i))/2;
                v8 = (max_matrix(j, i) + max_matrix(j, i - 1))/2;
                v9 = (max_matrix(j, i) + 0)/2;
                
            elseif j == 1  
                v1 = (0 + min_matrix(j, i))/2;
                v2 = (min_matrix(j + 1, i) + min_matrix(j, i))/2;
                v3 = (min_matrix(j, i) + 0)/2;
                v4 = (min_matrix(j, i) + min_matrix(j, i + 1))/2;
                
                v6 = (0 + max_matrix(j, i))/2;
                v7 = (max_matrix(j + 1, i) + max_matrix(j, i))/2;
                v8 = (max_matrix(j, i) + 0)/2;
                v9 = (max_matrix(j, i) + max_matrix(j, i + 1))/2;
                              
            elseif j == n
                
                v1 = (min_matrix(j - 1, i) + min_matrix(j, i))/2;
                v2 = (0 + min_matrix(j, i))/2;              
                v3 = (min_matrix(j, i) + min_matrix(j, i - 1))/2;
                v4 = (min_matrix(j, i) + min_matrix(j, i + 1))/2;
                
                v6 = (max_matrix(j - 1, i) + max_matrix(j, i))/2;
                v7 = (0 + max_matrix(j, i))/2;
                v8 = (max_matrix(j, i) + max_matrix(j, i - 1))/2;
                v9 = (max_matrix(j, i) + max_matrix(j, i + 1))/2;            
            
            elseif i == 1   %left boundary
                v1 = (min_matrix(j - 1, i) + min_matrix(j, i))/2;
                v2 = (min_matrix(j + 1, i) + min_matrix(j, i))/2;
                v3 = (min_matrix(j, i) + 0)/2;
                v4 = (min_matrix(j, i) + min_matrix(j, i + 1))/2;
                
                v6 = (max_matrix(j - 1, i) + max_matrix(j, i))/2;
                v7 = (max_matrix(j + 1, i) + max_matrix(j, i))/2;
                v8 = (max_matrix(j, i) + 0)/2;
                v9 = (max_matrix(j, i) + max_matrix(j, i + 1))/2;              

            elseif i == m   %right boundary
                v1 = (min_matrix(j - 1, i) + min_matrix(j, i))/2;
                v2 = (min_matrix(j + 1, i) + min_matrix(j, i))/2;
                v3 = (min_matrix(j, i) + min_matrix(j, i - 1))/2;
                v4 = (min_matrix(j, i) + 0)/2;
                
                v6 = (max_matrix(j - 1, i) + max_matrix(j, i))/2;
                v7 = (max_matrix(j + 1, i) + max_matrix(j, i))/2;
                v8 = (max_matrix(j, i) + max_matrix(j, i - 1))/2;
                v9 = (max_matrix(j, i) + 0)/2;               

            else
                v1 = (min_matrix(j - 1, i) + min_matrix(j, i))/2;
                v2 = (min_matrix(j + 1, i) + min_matrix(j, i))/2;
                v3 = (min_matrix(j, i) + min_matrix(j, i - 1))/2;
                v4 = (min_matrix(j, i) + min_matrix(j, i + 1))/2;
                
                v6 = (max_matrix(j - 1, i) + max_matrix(j, i))/2;
                v7 = (max_matrix(j + 1, i) + max_matrix(j, i))/2;
                v8 = (max_matrix(j, i) + max_matrix(j, i - 1))/2;
                v9 = (max_matrix(j, i) + max_matrix(j, i + 1))/2;
                
            end
            
            v5 = min_matrix(j, i);
            v10 = max_matrix(j, i);
            arr = [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10];
            lb = min(arr);
            ub = max(arr);                
            origin = [xlist(i)-deltax/2, tlist(j)-deltat/2, lb];
            edge = [deltax, deltat, ub];
            plotcube(edge, origin, 1, [0 0.5 0.5]);
            hold;                        
        end
   end
end