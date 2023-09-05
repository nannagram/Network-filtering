function A = ER_network(N,p)

    %%% Generating an ER network adjacency matrix in 4 steps
    %%% 1. "rand(N) < p" generates an NxN matrix of random numbers in [0,1],
    %%% and sets its entries that are below p to 1 and those that are above
    %%% p to 0.
    %%% 2. "sparse" converts the above matrix to a sparse object, from now
    %%% on everithing will be in the form of a sparse object (list of
    %%% indeces and value)
    %%% 3. triu(...,1) extracts the upper triangular part of the matrix,
    %%% the ",1" means that it excludes the main diagonal 
    %%% 4. A+A' symmetrizes the matrix

    A = triu(sparse(rand(N) < p),1);
    A = A+A';

end
