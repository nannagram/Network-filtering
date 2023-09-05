clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file generates an ER network, then assign the weights to the links
%%% and then save the network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% GENERATE THE NETWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 3000; % Number of nodes
p = 0.005; % Connectivity parameter

A = ER_net_directed(N,p); % Generating an ER network's adjacency matrix

disp("end network generation")

%%% RANDOMLY ASSIGN THE WEIGHTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = sparse(N,N);  %initialize weights matrix

% mu = 70;            % power law exponent
% sigma = 10;           % power law coefficient 

L = nnz(A);  % count the number of links

% % randomly extract the weights out of a gaussian distribution 
% rand = normrnd(mu, sigma, L,1); 

% randomly extract the weights out of a power-law distribution
a = 2;            % power law exponent
k = 1;           % power law coefficient 

r_weight = randraw('pareto', [k,a], L); 

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

disp("end weight generation")

save('/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/saveweightedER_powelaw.mat', "W")

disp("end saving network")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k_in = full(sum(A));   % Computing the incoming degree sequence
k_out = full(sum(A')); % Computing the outgoing degree sequence

k = k_in + k_out;

subplot(1,2,1)
histogram(k_in, 'Normalization', 'pdf');

hold on

subplot(1,2,2)
histogram(k_out);


%%% VISUALIZING THE NETWORK

G = digraph(A); %first create the graph object from the adiacience matric #F5D9A4 #AC291E


% plot con link lilla semitrasparenti, nodi viola e piccoli
plot(G, 'Layout','auto', 'EdgeColor',[162/255 158/255 254/255], ...
    'EdgeAlpha', 0.2, 'NodeColor',[105/255 95/255 255/255], 'MarkerSize',.9)

% plot(G, 'Layout','auto', 'EdgeColor','white', 'NodeColor','#F5D9A4') 
% %the layout specifics how to dispose the nodes
% set(gca,'Color','k')




