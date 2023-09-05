function A = ER_net_directed(N,p)
%%% this function generate the adiacency matrix of a directed ER network
    A = sparse(rand(N) < p);
end
