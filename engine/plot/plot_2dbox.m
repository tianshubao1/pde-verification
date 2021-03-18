% function plot_2dbox(min_list, max_list, sample_list1, deltax, xlist)
function plot_2dbox(min_list, max_list, deltax, xlist)
    m = size(xlist);
    m = m(2);
    figure;
    title('2-Dimension Reachable Sets at t = 4s')
    xlabel('Distance x')
    ylabel('u(x,4)')
    

    for i = 1 : m 
        if i == 1
            v1 = (0 + min_list(i))/2;
            v2 = (0 + max_list(i))/2;
            v3 = min_list(i);
            v4 = max_list(i);
            v5 = (min_list(i) + min_list(i + 1))/2;
            v6 = (max_list(i) + max_list(i + 1))/2;
            arr = [v1,v2,v3,v4,v5,v6];
            lb = min(arr);
            ub = max(arr);
            rectangle('Position',[xlist(i)-deltax/2   lb  deltax  ub-lb], 'FaceColor',[0 0.5 0.5], 'EdgeColor',[0 0.5 0.5]); 
            hold;
            
        elseif i == m
            v1 = (min_list(i - 1) + min_list(i))/2;
            v2 = (max_list(i - 1) + max_list(i))/2;
            v3 = min_list(i);
            v4 = max_list(i);
            v5 = (min_list(i) + 0)/2;
            v6 = (max_list(i) + 0)/2;
            arr = [v1,v2,v3,v4,v5,v6];
            lb = min(arr);
            ub = max(arr);
            rectangle('Position',[xlist(i)-deltax/2   lb  deltax  ub-lb], 'FaceColor',[0 0.5 0.5], 'EdgeColor',[0 0.5 0.5]);     
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
            rectangle('Position',[xlist(i)-deltax/2   lb  deltax  ub-lb], 'FaceColor',[0 0.5 0.5], 'EdgeColor',[0 0.5 0.5]);
            hold;
        end
    end
    
%------plot interior trajectories---------------------%    
%     hold on; 
%     plot(xlist, sample_list1, 'black');
%     hold on;    
%     plot(xlist, sample_list2, 'black');
%     hold on;
%     plot(xlist, sample_list3, 'black'); 
    

end