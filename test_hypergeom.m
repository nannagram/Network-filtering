clear all
close all

addpath("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/PLOT/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/PLOTPLOT/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/MAX&SAM/MAX&SAM-routines/")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file aims to test the dependence of the hypergeometric filter with
% the total strength of the network. 
% to do so, we build an ER network, and then we assign to each link a
% weight extracted from a PARETO distribution with a normalisation
% coefficient that that span several order of magnitude
% The expected behaviour is a saturation of the fraction of validated link
% to 1 when the total stregth became too big.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% GENERATE THE NETWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 3000; % Number of nodes
p = 0.01; % Connectivity parameter

ER = ER_net_directed(N,p); % Generating an ER network's adjacency matrix

disp("end network generation")

%%% PARAMETERS
alpha = 0.05;    % univariate significance level
L = nnz(ER);      % Computing the number of links in A
alpha = alpha/L; % Bonferroni correction

iter = 20;
k = logspace(-1, 2, iter); % vector of values for the k parameter of pareto distribution 
a = 2;    

fr_link = [];
fr_nodes = [];
S = [];

%%% ASSIGN THE WEIGHTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aux = find(ER>0);

          

W = sparse(N,N);  %initialize weights matrix

for i = 1:iter
    disp(['Filtro', num2str(i), 'di', num2str(iter)])
    rand = randraw('pareto', [k(i),a], L); 
    W(aux) = rand;

    [backbone, pvalues] = hypergeom_filter(W,alpha);  
    fr_nodes(i) = length(unique(backbone)) / size(W,1);
    fr_link(i) = length(backbone) / L;
    S(i) = sum(sum(W)); 
end

figure(1)
semilogx(k, fr_link, k, fr_nodes, 'LineWidth',1.5);

figure(2)
semilogx(S, fr_link, S, fr_nodes, 'LineWidth',1.5);
xlim([S(1), S(iter)])
xlabel('Total strength S', 'interpreter', 'latex')
ylabel('Fraction of links/nodes in the backbone', 'interpreter', 'latex')
title('Dipendence of Hypergeometric filter from total strength', 'interpreter', 'latex')
legend('fraction of links', 'fraction of nodes', 'Location', 'NorthEast', 'interpreter', 'latex')
set(gca, 'fontsize', 20, 'ticklabelinterpreter', 'latex')






function A = ER_net_directed(N,p)
%%% this function generate the adiacency matrix of a directed ER network
    A = sparse(rand(N) < p);
end

