clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file generates and saves a weighted directed Barabasi-Albert network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% 1) GENERATE THE NETWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 4000;   % Number of final nodes
m0 = 6;     % Starting fully connected nodes
m = 4;      % Number of new links added at each step
p = 1;      % probability to connect to the randomly chosen link (1-p is the prob to copy one of the links of the chosen node)
alfa = 0.5; % probability of generating a link outgoing from the new node

A = BA_Dnet(N, m0, m, p, alfa);

disp('network generated')

%%% 2) RANDOMLY ASSIGN THE WEIGHTS (POWER LAW) %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% W = sparse(N,N);  %initialize weights matrix
% 
% a = 2.5;            % power law exponent
% k = 1;           % power law coefficient 
% 
% L = nnz(A);  % count the number of links
% 
% % randomly extract the weights out of a power law
% r_weight = randraw('pareto', [k,a], L); 
% 
% disp('extracted weights list')
% 
% aux = find(A>0);
% 
% W(aux) = r_weight;
% 
% save('/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/saveweightedBA.mat', "W")
% 
% disp('saved network')

%%% 2) RANDOMLY ASSIGN THE WEIGHTS (GAUSSIAN) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Y = sparse(N,N);  %initialize weights matrix
% 
% mu = 100;            % mean 
% sigma = 10;          % variance 
% 
% L = nnz(A);  % count the number of links
% 
% % randomly extract the weights out of a gaussian distribution
% rand = normrnd(mu, sigma, L,1); 
% rand(rand<0) = normrnd(mu, sigma, 1,1);
% 
% disp('extracted second weights list')
% 
% aux = find(A>0);
% 
% disp('weights assigned')
% 
% Y(aux) = rand;
% 
% save('/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/saveweightedBA_gauss.mat', "W")
% 
% disp('second network saved')


%%% 3) PARTIALLY RANDOMLY ASSIGN THE WEIGHTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extract the weight list, then sort it from heaviest to smallest
%%% extract the link to which assign the weights based on the degree of the
%%% two nodes. 
%%% statistically the heavier links will go to higher connected nodes. 

W = sparse(N,N);  %initialize weights matrix

a = 2;         % power law exponent
k = 1;           % power law coefficient 

L = nnz(A);      % count the number of links

r_weight = randraw('pareto', [k,a], L); % extract the weights out of a power law

r_weight = sort(r_weight, 'descend'); % sort the weights from highest to lowest

k_in = sum(A);    % get inward and outward degree of each nodes
k_out= sum(A');

p_in = k_in ./ (sum(k_in)); % get probabilities based on degrees
p_out = k_out ./ (sum(k_out));
 
[row, col] = find(A>0);     % get indeces of links

p_link = p_out(row) .* p_in(col); % get probability of each link
p_link = p_link ./ (sum(p_link));

p = cumsum(p_link);         % get the cumulative for the random extraction 

while isempty(r_weight) == 0 % check there are still weights not assigned 

    r   = rand(1,1);     % extract a random number 
    aux = p - r;         % compute differences with p comulatives 
    aux = find(aux > 0); % find the ones that are 

    % check if the link found already hasn't a weight assigned
    if W(row(aux(1)),col(aux(1))) == 0
        % assign the weight
        W(row(aux(1)),col(aux(1))) = r_weight(1);
        % remove that weight from the list 
        r_weight(1) = [];
    else 
        continue 
    end

end

save('/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/saveweightedBA2.mat', "W")

disp('saved network')


%%% PLOT DEGREE DISTRIBUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if you want to plot the network with the nodes size proportional to the
% strength of the node:
% plot(G, 'Layout','auto', 'EdgeColor',[162/255 158/255 254/255], 'EdgeAlpha', 0.2, 'NodeColor',[105/255 95/255 255/255], 'MarkerSize', s/2.5e5)

% single index of weigths in the adiacency matrix
wind = sub2ind(size(W), row, col);

h = histogram(W(wind), 50, 'Normalization','pdf');

% r_weight = randraw('pareto', [k,a], 100000); 
% h = histogram(r_weight, 100, 'Normalization','pdf');

x = h.BinEdges;
x = x(2:end);
y = h.Values;

% subplot(2,1,1)
loglog(x,y,'ob','MarkerSize',6,'MarkerFaceColor','b')
hold on
% power-law fit
loglog(x, a*k^a*x.^-(a+1), 'b-','LineWidth',1.5)


% xlabel('$x$','Interpreter','latex')
% ylabel('$p(x)$','Interpreter','latex')
% set(gca,'FontSize',22)
% 
% subplot(2,1,2)
% histogram(Y(aux), 'BinCounts',100, 'Normalisation', 'pdf')



