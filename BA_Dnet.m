function A = BA_Dnet(N, m0, m, p, alfa)
%%% function to generate a DIRECTED BA network 
%%% N total number of nodes at the end
%%% m0 the initial number of nodes
%%% m links for each new node
%%% p the probability that drives the copying attachement 
%%% alfa the probability to generate a link from the new node from the
%%% existing one 

    % Initialising adjacency matrix of initial core (fully connected)
    A = ones(m0) - eye(m0);
    A = sparse(A);

    % Loop on number of nodes to be added 
    for i = 1:N-m0
    
        k_in = full(sum(A)); % Incoming degree sequence (sum over a column)
        p_in = k_in/sum(k_in); % Probability of linking to a node

        k_out = full(sum(A')); % Outgoing degree (sum over a line) 
        p_out = k_out / sum(k_out); % Probability of starting a link
    
        %%% Our strategy to simulate preferential attachement will be
        %%% to convert the vector p of probabilities into a vector of
        %%% cumulative probabilities (i.e., a vector that sums up to 1,
        %%% which is akin to the partition of the unit interval into k
        %%% sub-intervals, whose length is proportional to the
        %%% probability of connecting to nodes).
        
        p_in = cumsum(p_in);  % cumulative probability of linking to a node 
        p_out = cumsum(p_out); % cumulative prob of starting a link 

        r = rand(m,1); % m random numbers in [0,1]
    
        %%% Adding new nodes and links
        ind1 = [];
        ind2 = [];
    
        %%% Loop on new links formed by new node: the links are formed
        %%% with the preexisting nodes whose corresponding
        %%% sub-interval in the vector p_in contains the random numbers in
        %%% r
        %%% we use p_in since we want to introduce links
        for j = 1:m
            % chose the direction of the link to add
            if rand < alfa
                % randomly select a node, weighted on his incoming degree        
                aux1 = p_in - r(j);
                aux1 = find(aux1 > 0);
                ind1 = [ind1; aux1(1)];
            else 
                %randomly select a node, weigthed on his outgoing degree
                aux2 = p_out - r(j);
                aux2 = find(aux2 > 0);
                ind2 = [ind2; aux2(1)];
            end
    
        end
    
        ind1 = unique(ind1);
        ind2 = unique(ind2);
        
        %%% Creating new rows and columns in adjacency matrix
        A = [A; zeros(1,size(A,2))];
        A = [A zeros(size(A,1),1)];

        %%% FOR THE LINK FROM EXISTING TO NEW
        A(ind2, end) = 1;

        %%% FOR THE LINK FROM NEW TO EXISTING 
        if length(ind1) > 0 
            n_rand = rand(length(ind1),1); % generate a rand num for each new link
        
            for k = 1:length(n_rand) % loop sui 3 link da aggiungere 
        
                if n_rand(k) < p % collego il link al nodo scelto con prob p
                    A(end, ind1(k)) = 1;
        
                else
                    % find the nodes tho whom the selected one is connected
                    f = find(A(ind1(k), :) == 1);
                    % se il nodo scelto non ha link in uscita, connetti al
                    % nodo stesso
                    if length(f) == 0
                        A(end, ind1(k)) = 1;
                    else
                         % select a random node among these
                        index = randsample(f,1);
                        % copy the link
                        A(end, index) = 1;
                    end
                end
            end  
        end
    end     
end
