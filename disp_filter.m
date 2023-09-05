function [b, pvalues] = disp_filter(A, alpha)
% this function returns the list of p-values computed and the edge list 
% of the backbone validated by the disparity filter for a directed network

    k_in  = full(sum(A > 0));   % IN Degree sequence
    k_out = full(sum(A' > 0));  % OUT Degree sequence
    s_in  = full(sum(A));       % IN Strength sequence
    s_out  = full(sum(A'));     % OUT Strength sequence

    % Finding indices of non-zero entries in A (i.e., links)
    [ind1,ind2] = find(A > 0); 

    b = []; % Empty array to store links in backbone
    pvalues = []; % Empty array to store p-values

    % Loop on links
    for i = 1:length(ind1) 
        
        w = A(ind1(i),ind2(i)); % Weight on current link
    
        % P-value of current link relative to the destination node
        p_in = (1 - w/s_in(ind2(i)))^k_in(ind2(i)); 

        % P-value of current link relative to the starting node
        p_out = (1 - w/s_out(ind1(i)))^k_out(ind1(i)); 
        
        % chek on the p-values to find the lower one 
        if k_in(ind2(i)) == 1
            P = p_out;
        elseif k_out(ind1(i)) == 1
            P = p_in;
        else 
            P = min(p_in, p_out); %keep the minor one 
        end

        pvalues = [pvalues; P];

        % If the p-value falls below the significance level in input, the
        % corresponding link is stored in the backbone
        if P < alpha
           b = [b; ind1(i) ind2(i)];
        end
end





function b = disp_filter_old(A,alpha)
% this function returns the edge list of the backbone validated by the
% disparity filter for an undirected network

    k = full(sum(A > 0)); % Degree sequence
    s = full(sum(A)); % Strength sequence

    [ind1,ind2] = find(A > 0); % Finding indices of non-zero entries in A (i.e., links)

    b = []; % Empty array to store links in backbone

    for i = 1:length(ind1) % Loop on links
        
        w = A(ind1(i),ind2(i)); % Weight on current link
    
        p = (1 - w/s(ind1(i)))^k(ind1(i)); % P-value of current link
        
        % If the p-value falls below the significance level in input, the
        % corresponding link is stored in the backbone
        if p < alpha
           b = [b; ind1(i) ind2(i)];
        end

end