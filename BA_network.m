function A = BA_network(N,m)    

        % Initialising adjacency matrix of initial core (fully connected)
        A = ones(m) - eye(m);
        A = sparse(A);
        
        % Loop on number of nodes to be added 
        for i = 1:N-m
        
            k = full(sum(A)); % Degree sequence
            p = k/sum(k); % Probability of linking to a node
        
            %%% Our strategy to simulate preferential attachement will be
            %%% to convert the vector p of probabilities into a vector of
            %%% cumulative probabilities (i.e., a vector that sums up to 1,
            %%% which is akin to the partition of the unit interval into k
            %%% sub-intervals, whose length is proportional to the
            %%% probability of connecting to nodes).
            
            p = cumsum(p); 
            r = rand(m,1); % m random numbers in [0,1]
        
            %%% Adding new nodes and links
            ind = [];
        
            %%% Loop on new links formed by new node: the links are formed
            %%% with the preexisting nodes whose corresponding
            %%% sub-interval in the vector p contains the random numbers in
            %%% r
            for j = 1:m
                        
                aux = p - r(j);
                aux = find(aux > 0);
                ind = [ind; aux(1)];
        
            end
        
            ind = unique(ind);
        
            %%% Creating new rows and columns in adjacency matrix
            A = [A; zeros(1,size(A,2))];
            A = [A zeros(size(A,1),1)];
        
            %%% Adding new links
            A(end,ind) = 1;
            A(ind,end) = 1;
        
        end     

end