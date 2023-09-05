clear all
close all 

addpath("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/PLOT/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/PLOTPLOT/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/MAX&SAM/MAX&SAM-routines/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/3_BA_NET/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/1_ER_NET/")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file aims to plot measurement of the different performance of the
%%% 4 filters implemented in response to the heterogeneity of the networks,
%%% in function of the powerlaw exponent of the weight distribution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% GENERATE THE TOPOLOGY OF THE NETWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% BA NETWORK %%%%%%%%%%%%%%%%%%%%
N = 6000;   % Number of final nodes
m0 = 6;     % Starting fully connected nodes
m = 4;      % Number of new links added at each step
p = 1;      % probability to connect to the randomly chosen link (1-p is the prob to copy one of the links of the chosen node)
alfa = 0.5; % probability of generating a link outgoing from the new node

A = BA_Dnet(N, m0, m, p, alfa);

disp('network generated')


% %%% ER NETWORK %%%%%%%%%%%%%%%%%%%%
% N = 3000; % Number of nodes
% p = 0.005; % Connectivity parameter
% 
% A = ER_net_directed(N,p); % Generating an ER network's adjacency matrix
% disp("end network generation")


L = nnz(A);  % count the number of links
aux = find(A>0);

%%% GENERATE THE WEIGHTED NETWORK AND FILTER IT %%%%%%%%%%%%%%%%%%%%%%%%%%%

fr_node_DF     = [];   % disparity filter
fr_node_HF     = [];   % hypergeometric filter
fr_node_PF_min = [];   % polya filter min 
fr_node_PF_max = [];   % polya filter max 
fr_node_PF_1   = [];   % polya filter a=1
fr_node_PF_ML  = [];   % polya filter a maximum likelihood

fr_link_DF     = [];   % disparity filter
fr_link_HF     = [];   % hypergeometric filter
fr_link_PF_min = [];   % polya filter min 
fr_link_PF_max = [];   % polya filter max
fr_link_PF_1   = [];   % polya filter a=1
fr_link_PF_ML  = [];   % polya filter a maximum 

fr_weig_DF     = [];   % disparity filter
fr_weig_HF     = [];   % hypergeometric filter
fr_weig_PF_min = [];   % polya filter min 
fr_weig_PF_max = [];   % polya filter max
fr_weig_PF_1   = [];   % polya filter a=1
fr_weig_PF_ML  = [];   % polya filter a maximum 



%%% PARAMETERS OF THE WEIGHT DISTRIBUTION 
a = linspace(1, 3, 20); % array of 10 exponents 
k = 1;           % power law coefficient 

%%% PARAMETERS OF THE FILTERS
alpha   = 0.05/L;  % Bonferroni corrected statical significance level
apr_lvl = 1e23;    % approximation level for the polya filter (NO APPROXIMATION)
a1      = 0.01;    % min val for a in polya filter 
a2      = 2;       % max val for a in polya filter 

aml = [];


for i=1:length(a)

    disp(i)

    % initialize matrix
    W = sparse(N,N);
    
    % randomly extract the weights out of a power law
    r_weight = randraw('pareto', [k,a(i)], L); 

    disp('extracted weights list')

    W(aux) = r_weight;

    S = sum(sum(W)) % total strength

    disp('weighted network completed')

    %%% FILTER THE NETWORKS

    [bb_DF, p_DF]            = disp_filter(W, alpha);
    [bb_HF, p_HF]            = hypergeom_filter(W, alpha);
    [bb_PF_min, p_PF_min]    = PF(W, a1, alpha, apr_lvl, 0);
    [bb_PF_max, p_PF_max]    = PF(W, a2, alpha, apr_lvl, 0);
    [bb_PF_1, p_PF_1]        = PF(W, 1, alpha, apr_lvl, 0);
    [bb_PF_ML, p_PF_ML, mle] = PF(W, -1, alpha, apr_lvl, 0);

    aml(i) = mle; 


    %%% COMPUTE THE FRACTION OF LINKS 
    fr_link_DF(i)     = length(bb_DF)/L;
    fr_link_HF(i)     = length(bb_HF)/L;
    fr_link_PF_min(i) = length(bb_PF_min)/L;
    fr_link_PF_max(i) = length(bb_PF_max)/L;
    fr_link_PF_1(i)   = length(bb_PF_1)/L;
    fr_link_PF_ML(i)  = length(bb_PF_ML)/L;


    %%% COMPUTE THE FRACTION OF NODES
    [ind_HF] = sub2ind(size(W), bb_HF(:,1), bb_HF(:,2));
    fr_node_HF(i) = length(unique(bb_HF))/N;

    if isempty(bb_DF) == 0
        [ind_DF] = sub2ind(size(W), bb_DF(:,1), bb_DF(:,2));
        fr_node_DF(i) = length(unique(bb_DF))/N;
    else 
        ind_DF = [];
        fr_nodes_DF(i) = 0;
    end

    if isempty(bb_PF_min) == 0
        [ind_PF_min] = sub2ind(size(W), bb_PF_min(:,1), bb_PF_min(:,2));
        fr_node_PF_min(i) = length(unique(bb_PF_min))/N;
    else
        ind_PF_min = [];
        fr_node_PF_min(i) = 0;
    end


    if isempty(bb_PF_max) == 0
        [ind_PF_max] = sub2ind(size(W), bb_PF_max(:,1), bb_PF_max(:,2));
        fr_node_PF_max(i) = length(unique(bb_PF_max))/N;
    else
        ind_PF_max = [];
        fr_node_PF_max(i) = 0;
    end

    if isempty(bb_PF_1) == 0
        [ind_PF_1] = sub2ind(size(W), bb_PF_1(:,1), bb_PF_1(:,2));
        fr_node_PF_1(i) = length(unique(bb_PF_1))/N;
    else 
        ind_PF_1 = [];
        fr_node_PF_1(i) = 0;
    end
 
    [ind_PF_ML] = sub2ind(size(W), bb_PF_ML(:,1), bb_PF_ML(:,2));
    fr_node_PF_ML(i) = length(unique(bb_PF_ML))/N;


    %%% COMPUTE THE FRACTION TOTAL WEIGHT
    fr_weig_DF(i)     = sum(sum(W(ind_DF))) / S;
    fr_weig_HF(i)     = sum(sum(W(ind_HF))) / S;
    fr_weig_PF_min(i) = sum(sum(W(ind_PF_min))) / S;
    fr_weig_PF_max(i) = sum(sum(W(ind_PF_max))) / S;
    fr_weig_PF_1(i)   = sum(sum(W(ind_PF_1))) / S;
    fr_weig_PF_ML(i)  = sum(sum(W(ind_PF_ML))) / S;

end


%%% PLOT THE FRACTION OF LINKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)

semilogy(a, fr_link_DF, a, fr_link_HF, ...
    a, fr_link_PF_max, a, fr_link_PF_min, ...
    a, fr_link_PF_1, a, fr_link_PF_ML, 'LineWidth',1.5);

title('Responsiveness to heterogeneity')
xlabel("Pareto exponent a")
ylabel("Fraction of links retained")
legend('DF', 'HF', 'PF_{MAX}', 'PF_{MIN}', 'PF (a=1)', ' PF (a MLE)', 'Location','southwest')

set(gca, 'FontName', 'latex')
savefig(gcf, 'PLOT/Het_Fr_links_BA.fig')



%%% PLOT THE FRACTION OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)

plot(a, fr_node_DF, a, fr_node_HF, ...
    a, fr_node_PF_max, a, fr_node_PF_min, ...
    a, fr_node_PF_1, a, fr_node_PF_ML, 'LineWidth',1.5);

title('Responsiveness to heterogeneity')
xlabel("Pareto exponent a")
ylabel("Fraction of nodes retained")
legend('DF', 'HF', 'PF_{MAX}', 'PF_{MIN}', 'PF (a=1)', ' PF (a MLE)', 'Location','southeast')

set(gca, 'FontName', 'latex')
savefig(gcf, 'PLOT/Het_Fr_node_BA.fig')



%%% PLOT THE FRACTION OF WEIGHT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)

plot(a, fr_weig_DF, a, fr_weig_HF, ...
    a, fr_weig_PF_max, a, fr_weig_PF_min, ...
    a, fr_weig_PF_1, a, fr_weig_PF_ML, 'LineWidth',1.5);

title('Responsiveness to heterogeneity')
xlabel("Pareto exponent a")
ylabel("Fraction of weight retained")
legend('DF', 'HF', 'PF_{MAX}', 'PF_{MIN}', 'PF (a=1)', ' PF (a MLE)', 'Location','northeast')

set(gca, 'FontName', 'latex')
savefig(gcf, 'PLOT/Het_Fr_weight_BA.fig')


%%% PLOT THE POLYA PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)

semilogy(a, aml, '.-' ,'color', [162/255 158/255 254/255], ...
    'Linewidth', 1, 'MarkerSize', 20, 'MarkerEdgeColor', [105/255 95/255 255/255])

title('Polya parameter  vs Pareto exponent')
xlabel('Pareto exponent')
ylabel('Polya parameter a_{ML}')

set(gca, 'FontName', 'latex')
savefig(gcf, 'PLOT/Het_PolyaPar.fig')


%%% PLOT THE PARETO DISTRIBUTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)

var = logspace(0, 7, 100); 
loglog(var, a(1).* k^a(1).*var.^(-a(1)-1), 'LineWidth',1)
hold on
for i = 2:length(a)
    loglog(var, a(i).* k^a(i).*var.^(-a(i)-1), 'LineWidth',1)
end

set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 15)
title('Pareto distribution', 'interpreter', 'latex', 'FontSize', 20)
xlabel('x', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$p(x \> | \> \gamma, x_m) $', 'FontSize',20, 'Interpreter', 'latex')
xlim([1, 1e7])
ylim([1e-28, 1e3])

savefig(gcf, 'PLOT/pareto.fig')






