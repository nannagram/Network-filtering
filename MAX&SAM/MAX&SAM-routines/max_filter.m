%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function is returns the p-value list and the backbone estracted by
%%% the max entropy filter for directed, weighted network according to the 
%%% weighted configuration model. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [b, p] = max_filter(A, solxy, alpha)
    % A is the weighted adiaciency matrix of the network
    % solxy is a .mat file containing the 2N lagrangian multiplier of the
    % solved system of equations 
    % alpha is the significance level

    sol = load(solxy);
    xy = sol.sol;

    % number of nodes
    N = length(A(:,1));
    % number of links
    L = nnz(A);

    % find indeces of links
    [ind1, ind2] = find(A>0);

    % initialize backbone
    b = [];

    % compute p-values
    p = zeros(L,1);
    for i=1:length(ind1)
        p(i) = (xy(ind1(i), 1).* xy(ind2(i), 2)).^A(ind1(i),ind2(i));

        % If the p-value falls below the significance level in input, the
        % corresponding link is stored in the backbone
        if p(i) < alpha
           b = [b; ind1(i) ind2(i)];  
        end        
    end
end