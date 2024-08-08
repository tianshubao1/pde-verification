function plot_fixlocboxerr(min_list, max_list, deltat, tlist, color)

    m = size(tlist);
    m = m(2);
    figure;
    title('Discrete Reachable Sets of Error at x = 4')
    xlabel('Time t')
%     ylabel('U_4^n')
    
    if isempty(color)
        color = [0 0 1];
    end
    
    for i = 1 : m 
        if i == 1
            v1 = min_list(i);
            v2 = max_list(i);
            v3 = min_list(i);
            v4 = max_list(i);
            v5 = (min_list(i) + min_list(i + 1))/2;
            v6 = (max_list(i) + max_list(i + 1))/2;
            arr = [v1,v2,v3,v4,v5,v6];
            lb = min(arr);
            ub = max(arr);
            rectangle('Position',[tlist(i)-deltat/3   lb  deltat*2/3  ub-lb], 'FaceColor',color, 'EdgeColor',color);      
            hold;
            
        elseif i == m
            v1 = (min_list(i - 1) + min_list(i))/2;
            v2 = (max_list(i - 1) + max_list(i))/2;
            v3 = min_list(i);
            v4 = max_list(i);
            v5 = min_list(i);
            v6 = max_list(i);
            arr = [v1,v2,v3,v4,v5,v6];
            lb = min(arr);
            ub = max(arr);
            rectangle('Position',[tlist(i)-deltat/3   lb  deltat*2/3  ub-lb], 'FaceColor',color, 'EdgeColor',color);       
            hold;
            
        else
            v1 = (min_list(i - 1) + min_list(i))/2;
            v2 = (max_list(i - 1) + max_list(i))/2;
            v3 = min_list(i);
            v4 = max_list(i);
            v5 = (min_list(i) + min_list(i + 1))/2;
            v6 = (max_list(i) + max_list(i + 1))/2;
            arr = [v1,v2,v3,v4,v5,v6];
            lb = min(arr);
            ub = max(arr);
            rectangle('Position',[tlist(i)-deltat/3   lb  deltat*2/3  ub-lb], 'FaceColor',color, 'EdgeColor',color);
            hold;
        end
    end
    
%     xlim([-0.1 0.1])
    ylim([-0.1 0.2])
end