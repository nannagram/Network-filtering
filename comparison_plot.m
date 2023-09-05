clear all 
close all

addpath("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/PLOT/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/PLOTPLOT/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/MAX&SAM/MAX&SAM-routines/")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file aims to plot measurement of the different performance of the
%%% 4 filters implemented and, in function of the multivariate
%%% significance level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% LOAD THE DESIRED NETWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % load the US air traffic network
% A = matfile("saveA.mat"); 
% A = A.A;

% load the WIOT network 
A = readmatrix('WIOT.txt');

% % load the RESIDENCE HALL network 
% A = readmatrix('Residence.txt');
% A = sparse(A(:,1), A(:,2), A(:,3));

% % load the FLORIDA FOODWEB IN THE DRY SEASON network 
% A = readmatrix('FloridaFoodweb.txt');
% A = sparse(A(:,1), A(:,2), A(:,3));

% % load the ER weighted network 3000 nodes
% ER = load('saveweightedER_powelaw.mat');
% A = ER.W;

% % load the BA weighted directed network 400 nodes 
% A = load("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/saveweightedBA400.mat");
% A = A.W;

% % load the BA weighted directed network 3000 nodes
% % A = load("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/saveweightedBA2.mat");
% % A = A.W;

% % load the BA weighted directed network 6000 nodes
% A = load("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/saveweightedBA_coupled.mat");
% A = A.W;





% GET THE NUMBER OF NODES CONNECTED IN THE NETWORK

N = length(A(:,1)); % total number of nodes

% compute the real number of nodes connected in the network
[row,col] = find(A>0);
new = [row col];
N_A = length(unique(new));

% GET THE NUMBER OF LINK IN THE NETWORK
S = sum(sum(A));
L = nnz(A); 
aux = find(A>0);

wsort = A(aux);
wsort = sort(wsort, 'descend');

%%% FRACTION OF NODES RETAINED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize vectors of nodes retained for 
fr_node_DF     = [];   % disparity filter
fr_node_HF     = [];   % hypergeometric filter
fr_node_PF_min = [];   % polya filter min 
fr_node_PF_max = [];   % polya filter max 
fr_node_PF_1   = [];   % polya filter a=1
fr_node_PF_ML   = [];   % polya filter a maximum likelihood
fr_node_ME     = [];   % max entropy filter 

%%% FRACTION OF LINKS RETAINED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize vectors of links retained for 
fr_link_DF     = [];   % disparity filter
fr_link_HF     = [];   % hypergeometric filter
fr_link_PF_min = [];   % polya filter min 
fr_link_PF_max = [];   % polya filter max
fr_link_PF_1   = [];   % polya filter a=1
fr_link_PF_ML  = [];   % polya filter a maximum likelihood 
fr_link_ME     = [];   % max entropy filter 

% initialize vectors of weight retained for 
fr_wei_DF     = [];   % disparity filter
fr_wei_HF     = [];   % hypergeometric filter
fr_wei_PF_min = [];   % polya filter min 
fr_wei_PF_max = [];   % polya filter max
fr_wei_PF_1   = [];   % polya filter a=1
fr_wei_PF_ML  = [];   % polya filter a maximum likelihood 
fr_wei_ME     = [];   % max entropy filter 

%%% JACCARD SIMILARITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize vectors of jaccard similarity for 
jacc_DF     = [];   % disparity filter
jacc_HF     = [];   % hypergeometric filter
jacc_PF_min = [];   % polya filter min 
jacc_PF_max = [];   % polya filter max
jacc_PF_1   = [];   % polya filter a=1
jacc_PF_ML   = [];   % polya filter a maximum likelihood
jacc_ME     = [];   % max entropy filter 

% initialize vectors of jaccard similarity for 
jacc_DF2     = [];   % disparity filter
jacc_HF2     = [];   % hypergeometric filter
jacc_PF_min2 = [];   % polya filter min 
jacc_PF_max2 = [];   % polya filter max
jacc_PF_12   = [];   % polya filter a=1
jacc_PF_ML2   = [];   % polya filter a maximum likelihood
jacc_ME2     = [];   % max entropy filter 

%%% SALIENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize vectors of salience for 
sal_DF     = [];   % disparity filter
sal_HF     = [];   % hypergeometric filter
sal_PF_min = [];   % polya filter min 
sal_PF_max = [];   % polya filter max
sal_PF_1   = [];   % polya filter a=1
sal_PF_ML  = [];   % polya filter a maximum likelihood
sal_ME     = [];   % max entropy filter 

%%% SALIENCE X JACCARD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize vectors of salience for 
O1_DF     = [];   % disparity filter
O1_HF     = [];   % hypergeometric filter
O1_PF_min = [];   % polya filter min 
O1_PF_max = [];   % polya filter max
O1_PF_1   = [];   % polya filter a=1
O1_PF_ML  = [];   % polya filter a maximum likelihood
O1_ME     = [];   % max entropy filter 



% generate num_alpha values of multivariate significance level 
% equally spaced between 10^-8 and 10^-1

num_alpha = 50;
alpha = logspace(-8, -1, num_alpha);


apr_lvl = 1e23; % approximation level for the polya filter (NO APPROXIMATION)
a1 = 0.01;      % min val for a in polya filter 
a2 = 6;         % max val for a in polya filter 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EXTRACT THE P-VALUE AND THE BACKBONE LIST FOR EACH FILTER 

% sol = Lag_mult(A);

% [bb_ME, p_ME]         = max_filter(A, 'solxy.mat', alpha(1));
[bb_DF, p_DF]         = disp_filter(A, alpha(1));
[bb_HF, p_HF]         = hypergeom_filter(A, alpha(1));
[bb_PF_min, p_PF_min] = PF(A, a1, alpha(1),apr_lvl, 0);
[bb_PF_max, p_PF_max] = PF(A, a2, alpha(1),apr_lvl, 0);
[bb_PF_1, p_PF_1]     = PF(A, 1, alpha(1),apr_lvl, 0);
[bb_PF_ML, p_PF_ML]   = PF(A, -1, alpha(1),apr_lvl, 0);



% initialize binary matrix to store the checks on the p-values based on
% each value of the significance level alpha 
% fin_ME     = zeros(L,length(alpha));
fin_DF     = zeros(L,length(alpha));
fin_HF     = zeros(L,length(alpha));
fin_PF_min = zeros(L,length(alpha));
fin_PF_max = zeros(L,length(alpha));
fin_PF_1   = zeros(L,length(alpha));
fin_PF_ML   = zeros(L,length(alpha));

% check whether p is under/above the significance level 
for j = 1:length(alpha)
%     fin_ME(:,j)     = (p_ME < alpha(j));
    fin_DF(:,j)     = (p_DF < alpha(j));
    fin_HF(:,j)     = (p_HF < alpha(j));
    fin_PF_min(:,j) = (p_PF_min < alpha(j));
    fin_PF_max(:,j) = (p_PF_max < alpha(j));
    fin_PF_1(:,j)   = (p_PF_1 < alpha(j));
    fin_PF_ML(:,j)  = (p_PF_ML < alpha(j));
end

%%% COMPUTE THE FRACTION OF LINKS RETAINED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fr_link_ME     = sum(fin_ME,1)/L;
fr_link_DF     = sum(fin_DF,1)/L;
fr_link_HF     = sum(fin_HF,1)/L;
fr_link_PF_min = sum(fin_PF_min,1)/L;
fr_link_PF_max = sum(fin_PF_max,1)/L;
fr_link_PF_1   = sum(fin_PF_1,1)/L;
fr_link_PF_ML  = sum(fin_PF_ML,1)/L;

%%% COMPUTE THE FRACTION OF NODES RETAINED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(alpha)
%     [ind_ME] = find(fin_ME(:,i)>0);
%     [row_ME, col_ME] = ind2sub(size(A), aux(ind_ME));
%     new_ME = [row_ME col_ME];
%     fr_node_ME(i) = length(unique(new_ME))/N_A;
%     wlist = sort(A(aux(ind_ME)), 'descend');
%     jacc_ME2(i) = jacc_weighted(wlist, wsort(1:length(wlist)));


    [ind_HF] = find(fin_HF(:,i)>0);
    [row_HF, col_HF] = ind2sub(size(A), aux(ind_HF));
    new_HF = [row_HF col_HF];
    fr_node_HF(i) = length(unique(new_HF))/N_A;
    wlist = sort(A(aux(ind_HF)), 'descend');
    jacc_HF2(i) = jacc_weighted(wlist, wsort(1:length(wlist)));


    [ind_DF] = find(fin_DF(:,i)>0);
    [row_DF, col_DF] = ind2sub(size(A), aux(ind_DF));
    new_DF = [row_DF col_DF];
    fr_node_DF(i) = length(unique(new_DF))/N_A;
    wlist = sort(A(aux(ind_DF)), 'descend');
    jacc_DF2(i) = jacc_weighted(wlist, wsort(1:length(wlist)));


    [ind_PF_min] = find(fin_PF_min(:,i)>0);
    [row_PF_min, col_PF_min] = ind2sub(size(A), aux(ind_PF_min));
    new_PF_min = [row_PF_min col_PF_min];
    fr_node_PF_min(i) = length(unique(new_PF_min))/N_A;
    wlist = sort(A(aux(ind_PF_min)), 'descend');
    jacc_PF_min2(i) = jacc_weighted(wlist, wsort(1:length(wlist)));


    [ind_PF_max] = find(fin_PF_max(:,i)>0);
    [row_PF_max, col_PF_max] = ind2sub(size(A), aux(ind_PF_max));
    new_PF_max = [row_PF_max col_PF_max];
    fr_node_PF_max(i) = length(unique(new_PF_max))/N_A;
    wlist = sort(A(aux(ind_PF_max)), 'descend');
    jacc_PF_max2(i) = jacc_weighted(wlist, wsort(1:length(wlist)));


    [ind_PF_1] = find(fin_PF_1(:,i)>0);
    [row_PF_1, col_PF_1] = ind2sub(size(A), aux(ind_PF_1));
    new_PF_1 = [row_PF_1 col_PF_1];
    fr_node_PF_1(i) = length(unique(new_PF_1))/N_A;
    wlist = sort(A(aux(ind_PF_1)), 'descend');
    jacc_PF_12(i) = jacc_weighted(wlist, wsort(1:length(wlist)));


    [ind_PF_ML] = find(fin_PF_ML(:,i)>0);
    [row_PF_ML, col_PF_ML] = ind2sub(size(A), aux(ind_PF_ML));
    new_PF_ML = [row_PF_ML col_PF_ML];
    fr_node_PF_ML(i) = length(unique(new_PF_ML))/N_A;
    wlist = sort(A(aux(ind_PF_ML)), 'descend');
    jacc_PF_ML2(i) = jacc_weighted(wlist, wsort(1:length(wlist)));

end



%%% COMPUTE THE BACKBONES AND THE MEASURES ON THEM %%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(alpha)
%     fprintf('Alpha value n. %d.\n', i);
    disp(i)

    tic
    %%% MAX ENTROPY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     disp("Max entropy filter")
%     % extract indeces of the backbone for this value of alpha 
%     % memo: fin is an array of size L (numebr of links
%     link_ind = find( fin_ME(:,i)>0 );
% 
%     % retreive right indeces of link selected 
%     back_row = row(link_ind); 
%     back_col = col(link_ind);
% 
%     % convert indeces of links to a one index 
%     link_ind = sub2ind(size(A), back_row, back_col);
% 
%     % build the backbone 
%     BB = sparse(back_row, back_col, A(link_ind));
% 
%     % pad the matrix with zeros to make it square
%     BB = sparse([full(BB) zeros(size(BB, 1), size(A,1)-size(BB,2))]);
%     BB = sparse([full(BB); zeros(size(A,2)-size(BB,1), size(BB,2))]);
%     
%     fr_wei_ME(i) = sum(sum(BB)) / S;
%     jacc_ME(i) = Jaccard_weighted(A, BB);
%     sal_ME(i)  = salience(BB);
%     O1_ME(i)   = jacc_ME(i)*sal_ME(i);


    %%% DISPARITY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp("Disparity filter")
    % extract indeces of the backbone for this value of alpha 
    % memo: fin is an array of size L (number of links)
    link_ind = find( fin_DF(:,i)>0 );

    % retreive right indeces of link selected 
    back_row = row(link_ind); 
    back_col = col(link_ind);

    % convert indeces of links to a one index 
    link_ind = sub2ind(size(A), back_row, back_col);

    % build the backbone 
    BB = sparse(back_row, back_col, A(link_ind));

    % pad the matrix with zeros to make it square
    BB = sparse([full(BB) zeros(size(BB, 1), size(A,1)-size(BB,2))]);
    BB = sparse([full(BB); zeros(size(A,2)-size(BB,1), size(BB,2))]);

    fr_wei_DF(i) = sum(sum(BB)) / S;
    jacc_DF(i) = Jaccard_weighted(A, BB);
    sal_DF(i)  = salience(BB);
    O1_DF(i)   = jacc_DF2(i)*sal_DF(i);

    

    %%% HYPERGEOMETRIC FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    disp("Hypergeometric filter")
    % extract indeces of the backbone for this value of alpha 
    % memo: fin is an array of size L (numebr of links
    link_ind = find( fin_HF(:,i)>0 );

    % retreive right indeces of link selected 
    back_row = row(link_ind); 
    back_col = col(link_ind);

    % convert indeces of links to a one index 
    link_ind = sub2ind(size(A), back_row, back_col);

    % build the backbone 
    BB = sparse(back_row, back_col, A(link_ind));

    % pad the matrix with zeros to make it square
    BB = sparse([full(BB) zeros(size(BB, 1), size(A,1)-size(BB,2))]);
    BB = sparse([full(BB); zeros(size(A,2)-size(BB,1), size(BB,2))]);

    fr_wei_HF(i) = sum(sum(BB)) / S;
    jacc_HF(i) = Jaccard_weighted(A, BB);
    sal_HF(i)  = salience(BB);
    O1_HF(i)   = jacc_HF2(i)*sal_HF(i);

    

    %%% POLYA FILTER A MIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    disp("Polya filter with a min")
    % extract indeces of the backbone for this value of alpha 
    % memo: fin is an array of size L (numebr of links
    link_ind = find( fin_PF_min(:,i)>0 );

    % retreive right indeces of link selected 
    back_row = row(link_ind); 
    back_col = col(link_ind);

    % convert indeces of links to a one index 
    link_ind = sub2ind(size(A), back_row, back_col);

    % build the backbone 
    BB = sparse(back_row, back_col, A(link_ind));

    % pad the matrix with zeros to make it square
    BB = sparse([full(BB) zeros(size(BB, 1), size(A,1)-size(BB,2))]);
    BB = sparse([full(BB); zeros(size(A,2)-size(BB,1), size(BB,2))]);

    fr_wei_PF_min(i) = sum(sum(BB)) / S;
    jacc_PF_min(i) = Jaccard_weighted(A, BB);
    sal_PF_min(i)  = salience(BB);
    O1_PF_min(i)   = jacc_PF_min2(i)*sal_PF_min(i);

    

    %%% POLYA FILTER A MAX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    disp("Polya filter with a max")
    % extract indeces of the backbone for this value of alpha 
    % memo: fin is an array of size L (numebr of links
    link_ind = find( fin_PF_max(:,i)>0 );

    % retreive right indeces of link selected 
    back_row = row(link_ind); 
    back_col = col(link_ind);

    % convert indeces of links to a one index 
    link_ind = sub2ind(size(A), back_row, back_col);

    % build the backbone 
    BB = sparse(back_row, back_col, A(link_ind));

    % pad the matrix with zeros to make it square
    BB = sparse([full(BB) zeros(size(BB, 1), size(A,1)-size(BB,2))]);
    BB = sparse([full(BB); zeros(size(A,2)-size(BB,1), size(BB,2))]);

    fr_wei_PF_max(i) = sum(sum(BB)) / S;
    jacc_PF_max(i) = Jaccard_weighted(A, BB);
    sal_PF_max(i)  = salience(BB);
    O1_PF_max(i)   = jacc_PF_max2(i)*sal_PF_max(i);

    

    %%% POLYA FILTER A = 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    disp("Polya filter with a = 1")
    % extract indeces of the backbone for this value of alpha 
    % memo: fin is an array of size L (numebr of links
    link_ind = find( fin_PF_1(:,i)>0 );

    % retreive right indeces of link selected 
    back_row = row(link_ind); 
    back_col = col(link_ind);

    % convert indeces of links to a one index 
    link_ind = sub2ind(size(A), back_row, back_col);

    % build the backbone 
    BB = sparse(back_row, back_col, A(link_ind));

    % pad the matrix with zeros to make it square
    BB = sparse([full(BB) zeros(size(BB, 1), size(A,1)-size(BB,2))]);
    BB = sparse([full(BB); zeros(size(A,2)-size(BB,1), size(BB,2))]);

    fr_wei_PF_1(i) = sum(sum(BB)) / S;
    jacc_PF_1(i) = Jaccard_weighted(A, BB);
    sal_PF_1(i)  = salience(BB);
    O1_PF_1(i)   = jacc_PF_12(i)*sal_PF_1(i);



    %%% POLYA FILTER A MAX LIKELIHOOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    disp("Polya filter with a  maximum likelihood")
    % extract indeces of the backbone for this value of alpha 
    % memo: fin is an array of size L (numebr of links
    link_ind = find( fin_PF_ML(:,i)>0 );

    % retreive right indeces of link selected 
    back_row = row(link_ind); 
    back_col = col(link_ind);

    % convert indeces of links to a one index 
    link_ind = sub2ind(size(A), back_row, back_col);

    % build the backbone 
    BB = sparse(back_row, back_col, A(link_ind));

    % pad the matrix with zeros to make it square
    BB = sparse([full(BB) zeros(size(BB, 1), size(A,1)-size(BB,2))]);
    BB = sparse([full(BB); zeros(size(A,2)-size(BB,1), size(BB,2))]);

    fr_wei_PF_ML(i) = sum(sum(BB)) / S;
    jacc_PF_ML(i) = Jaccard_weighted(A, BB);
    sal_PF_ML(i)  = salience(BB);
    O1_PF_ML(i)   = jacc_PF_ML2(i)*sal_PF_ML(i);

    toc
  
end


%%% COMPUTE O2 FOR EACH VALUE OF ALPHA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

O2_HF = sal_HF ./ jacc_HF2;
O2_DF = sal_DF ./ jacc_DF2;
O2_PF_min = sal_PF_min ./ jacc_PF_min2;
O2_PF_max = sal_PF_max ./ jacc_PF_max2;
O2_PF_1 = sal_PF_1 ./ jacc_PF_12;
O2_PF_ML = sal_PF_ML ./ jacc_PF_ML2;
O2_ME = sal_ME ./ jacc_ME2;




%%% PLOT GRAPH WITH FRACTION OF RETAINED NODES VS ALPHA %%%%%%%%%%%%%%%%%%%

loglog(alpha, fr_node_DF, alpha, fr_node_HF, alpha, ...
    fr_node_PF_max, alpha, fr_node_PF_min, ...
    alpha, fr_node_PF_1, alpha, fr_node_PF_ML, 'LineWidth',1.5);

% loglog(alpha, fr_node_DF, alpha, fr_node_HF, alpha, ...
%     fr_node_PF_max, alpha, fr_node_PF_min, ...
%     alpha, fr_node_PF_1, alpha, fr_node_PF_ML, alpha, fr_node_ME, 'LineWidth',1.5);

hold on 

fr_node_PF_max2 = fr_node_PF_max;
if nnz(fr_node_PF_max) < 50  
    fr_node_PF_max2(fr_node_PF_max==0) = min(fr_node_PF_max(fr_node_PF_max>0));
end

% set the 0 to a different value so the patch command works:
patch([alpha fliplr(alpha)], [fr_node_PF_min fliplr(fr_node_PF_max2)], [241/255 218/255 170/255],  'FaceAlpha',.3, 'LineStyle','none')

xline(0.05/L, '--')

set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 15)
title('WIOT network', 'interpreter', 'latex', 'fontsize', 20)
xlabel("Bonferroni-corrected $\alpha$", 'interpreter', 'latex', 'FontSize',20)
ylabel('$\frac{N_N^b}{N_N}$', 'interpreter', 'latex', 'rotation', 0, 'verticalalignment', 'middle', 'fontsize', 30, 'horizontalalignment', 'right')
legend('DF', 'HF', 'PF$_{max}$', 'PF$_{min}$', 'PF$_1$', ...
        ' PF$_{ML}$', '', '', 'Location','southoutside', ...
        'orientation', 'horizontal', 'interpreter', 'latex')
% legend('DF', 'HF', 'PF$_{max}$', 'PF$_{min}$', 'PF$_1$', ...
%         ' PF$_{ML}$', 'ME', '', '', 'Location','southoutside', ...
%         'orientation', 'horizontal', 'interpreter', 'latex')

xlim([1e-8 1e-1])

savefig(gcf, 'PLOT/Fr_nodes_WIOT.fig')

hold off 

%%% PLOT GRAPH WITH FRACTION OF RETAINED LINKS VS ALPHA %%%%%%%%%%%%%%%%%%%

loglog(alpha, fr_link_DF, alpha, fr_link_HF, ...
    alpha, fr_link_PF_max, alpha, fr_link_PF_min, ...
    alpha, fr_link_PF_1, alpha, fr_link_PF_ML, 'LineWidth',1.5);

% loglog(alpha, fr_link_DF, alpha, fr_link_HF, ...
%     alpha, fr_link_PF_max, alpha, fr_link_PF_min, ...
%     alpha, fr_link_PF_1, alpha, fr_link_PF_ML, alpha, fr_link_ME, 'LineWidth',1.5);

hold on 

% set the 0 to a different value so the patch command works:

fr_link_PF_max2 = fr_link_PF_max;
if nnz(fr_link_PF_max) < 50
    fr_link_PF_max2(fr_link_PF_max==0) = min(fr_link_PF_max(fr_link_PF_max>0));
end

fill([alpha fliplr(alpha)], [fr_link_PF_min fliplr(fr_link_PF_max2)], [241/255 218/255 170/255],  'FaceAlpha',.3, 'LineStyle','none')

xline(0.05/L, '--')

set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 15)
title('WIOT network', 'interpreter', 'latex', 'fontsize', 20)
xlabel("Bonferroni-corrected $\alpha$", 'interpreter', 'latex', 'FontSize',20)
ylabel('$\frac{N_L^b}{N_L}$', 'interpreter', 'latex', 'rotation', 0, 'verticalalignment', 'middle', 'fontsize', 30, 'horizontalalignment', 'right')
legend('DF', 'HF', 'PF$_{max}$', 'PF$_{min}$', 'PF$_1$', ...
        ' PF$_{ML}$', '', '', 'Location','southoutside', ...
        'orientation', 'horizontal', 'interpreter', 'latex')
% legend('DF', 'HF', 'PF$_{max}$', 'PF$_{min}$', 'PF$_1$', ...
%         ' PF$_{ML}$', 'ME', '', '', 'Location','southoutside', ...
%         'orientation', 'horizontal', 'interpreter', 'latex')

xlim([1e-8 1e-1])
savefig(gcf, 'PLOT/Fr_links_WIOT.fig')

hold off 



%%% PLOT GRAPH WITH FRACTION OF RETAINED WEIGHT VS ALPHA %%%%%%%%%%%%%%%%%%

loglog(alpha, fr_wei_DF, alpha, fr_wei_HF, ...
    alpha, fr_wei_PF_max, alpha, fr_wei_PF_min, ...
    alpha, fr_wei_PF_1, alpha, fr_wei_PF_ML, 'LineWidth',1.5);

% loglog(alpha, fr_wei_DF, alpha, fr_wei_HF, ...
%     alpha, fr_wei_PF_max, alpha, fr_wei_PF_min, ...
%     alpha, fr_wei_PF_1, alpha, fr_wei_PF_ML,alpha, fr_wei_ME, 'LineWidth',1.5);

hold on 

% set the 0 to a different value so the patch command works:
fr_wei_PF_max2 = fr_wei_PF_max;
if nnz(fr_wei_PF_max) < 50
    fr_wei_PF_max2(fr_wei_PF_max==0) = min(fr_wei_PF_max(fr_wei_PF_max>0));
end

fill([alpha fliplr(alpha)], [fr_wei_PF_min fliplr(fr_wei_PF_max2)], [241/255 218/255 170/255],  'FaceAlpha',.3, 'LineStyle','none')

xline(0.05/L, '--')

set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 15)
title('WIOT network', 'interpreter', 'latex', 'fontsize', 20)
xlabel("Bonferroni-corrected $\alpha$", 'interpreter', 'latex', 'FontSize',20)
ylabel('$\frac{S^b}{S}$', 'interpreter', 'latex', 'rotation', 0, 'verticalalignment', 'middle', 'fontsize', 30, 'horizontalalignment', 'right')
legend('DF', 'HF', 'PF$_{max}$', 'PF$_{min}$', 'PF$_1$', ...
        ' PF$_{ML}$', '', '', 'Location','southoutside', ...
        'orientation', 'horizontal', 'interpreter', 'latex')
% legend('DF', 'HF', 'PF$_{max}$', 'PF$_{min}$', 'PF$_1$', ...
%         ' PF$_{ML}$', 'ME', '', '', 'Location','southoutside', ...
%         'orientation', 'horizontal', 'interpreter', 'latex')

xlim([1e-8 1e-1])
savefig(gcf, 'PLOT/Fr_weight_WIOT.fig')

hold off 


%%% PLOT GRAPH WITH JACCARD SIMILARITY VS ALPHA %%%%%%%%%%%%%%%%%%%%%%%%%%%

semilogx(alpha, jacc_DF2, alpha, jacc_HF2, ...
    alpha, jacc_PF_max2, alpha, jacc_PF_min2, ...
    alpha, jacc_PF_12, alpha, jacc_PF_ML2, 'LineWidth',1.5);

% semilogx(alpha, jacc_DF2, alpha, jacc_HF2, ...
%     alpha, jacc_PF_max2, alpha, jacc_PF_min2, ...
%     alpha, jacc_PF_12, alpha, jacc_PF_ML2, alpha, jacc_ME, 'LineWidth',1.5);

hold on 

% set the 0 to a different value so the patch command works:
jacc_PF_max = jacc_PF_max2;
if nnz(jacc_PF_max2) < 50
    jacc_PF_max(isnan(jacc_PF_max2)) = max(jacc_PF_max2(jacc_PF_max2>0));
end

fill([alpha fliplr(alpha)], [jacc_PF_min2 fliplr(jacc_PF_max)], [241/255 218/255 170/255],  'FaceAlpha',.3, 'LineStyle','none')

xline(0.05/L, '--')

set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 15)
title('US air traffic network', 'interpreter', 'latex', 'fontsize', 20)
xlabel("Bonferroni-corrected $\alpha$", 'interpreter', 'latex', 'FontSize',20)
ylabel('J', 'interpreter', 'latex', 'rotation', 0, 'verticalalignment', 'middle', 'fontsize', 30, 'horizontalalignment', 'right')
legend('DF', 'HF', 'PF$_{max}$', 'PF$_{min}$', 'PF$_1$', ...
        ' PF$_{ML}$', '', '', 'Location','southoutside', ...
        'orientation', 'horizontal', 'interpreter', 'latex')
% legend('DF', 'HF', 'PF$_{max}$', 'PF$_{min}$', 'PF$_1$', ...
%         ' PF$_{ML}$', 'ME', '', '', 'Location','southoutside', ...
%         'orientation', 'horizontal', 'interpreter', 'latex')

xlim([1e-8 1e-1])
savefig(gcf, 'PLOT/Jacc_airport.fig')

hold off 

%%% PLOT GRAPH WITH SALIENCE VS ALPHA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

semilogx(alpha, sal_DF, alpha, sal_HF, ...
    alpha, sal_PF_max, alpha, sal_PF_min, alpha, ...
    sal_PF_1, alpha, sal_PF_ML, 'LineWidth',1.5);

% semilogx(alpha, sal_DF, alpha, sal_HF, ...
%     alpha, sal_PF_max, alpha, sal_PF_min, alpha, ...
%     sal_PF_1, alpha, sal_PF_ML, alpha, sal_ME, 'LineWidth',1.5);

hold on 

% set the 0 to a different value so the patch command works:
sal_PF_max2 = sal_PF_max;
if sum(isnan(sal_PF_max)) > 0
    sal_PF_max2(isnan(sal_PF_max)) = max(sal_PF_max);
end

patch([alpha fliplr(alpha)], [sal_PF_min fliplr(sal_PF_max2)], [241/255 218/255 170/255],  'FaceAlpha',.3, 'LineStyle','none')

xline(0.05/L, '--')

set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 15)
title('WIOT network', 'interpreter', 'latex', 'fontsize', 20)
xlabel("Bonferroni-corrected $\alpha$", 'interpreter', 'latex', 'FontSize',20)
ylabel('$\langle S \rangle$', 'interpreter', 'latex', 'rotation', 0, 'verticalalignment', 'middle', 'fontsize', 30, 'horizontalalignment', 'right')
legend('DF', 'HF', 'PF$_{max}$', 'PF$_{min}$', 'PF$_1$', ...
        ' PF$_{ML}$', '', '', 'Location','southoutside', ...
        'orientation', 'horizontal', 'interpreter', 'latex')
% legend('DF', 'HF', 'PF$_{max}$', 'PF$_{min}$', 'PF$_1$', ...
%         ' PF$_{ML}$', 'ME', '', '', 'Location','southoutside', ...
%         'orientation', 'horizontal', 'interpreter', 'latex')

xlim([1e-8 1e-1])
savefig(gcf, 'PLOT/sal_WIOT.fig')

hold off 

%%% PLOT GRAPH WITH O1 MEASURE VS ALPHA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% O1_PF_max(isnan(O1_PF_max)) = 1e-5;

semilogx(alpha, O1_DF, alpha, O1_HF, ...
    alpha, O1_PF_max, alpha, O1_PF_min, ...
    alpha, O1_PF_1, alpha, O1_PF_ML, 'LineWidth',1.5);

% loglog(alpha, O1_DF, alpha, O1_HF, ...
%     alpha, O1_PF_max, alpha, O1_PF_min, ...
%     alpha, O1_PF_1, alpha, O1_PF_ML, alpha, O1_ME, 'LineWidth',1.5);


hold on 

% set the 0 to a different value so the patch command works:
O1_PF_max2 = O1_PF_max;
if sum(isnan(O1_PF_max)) > 0
    O1_PF_max2(isnan(O1_PF_max)) = max(O1_PF_max(O1_PF_max>0));
end

patch([alpha fliplr(alpha)], [O1_PF_min fliplr(O1_PF_max2)], [241/255 218/255 170/255],  'FaceAlpha',.3, 'LineStyle','none')

xline(0.05/L, '--')

set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 15)
title('WIOT network', 'interpreter', 'latex', 'fontsize', 20)
xlabel("Bonferroni-corrected $\alpha$", 'interpreter', 'latex', 'FontSize',20)
ylabel('O$_1$', 'interpreter', 'latex', 'rotation', 0, 'verticalalignment', 'middle', 'fontsize', 30, 'horizontalalignment', 'right')
legend('DF', 'HF', 'PF$_{max}$', 'PF$_{min}$', 'PF$_1$', ...
        ' PF$_{ML}$', '', '', 'Location','southoutside', ...
        'orientation', 'horizontal', 'interpreter', 'latex')
% legend('DF', 'HF', 'PF$_{max}$', 'PF$_{min}$', 'PF$_1$', ...
%         ' PF$_{ML}$', 'ME', '', '', 'Location','southoutside', ...
%         'orientation', 'horizontal', 'interpreter', 'latex')

xlim([1e-8 1e-1])
savefig(gcf, 'PLOT/O1_WIOT.fig')

hold off 


%%% PLOT GRAPH WITH O2 MEASURE VS ALPHA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

semilogx(alpha, O2_DF, alpha, O2_HF, ...
    alpha, O2_PF_max, alpha, O2_PF_min, ...
    alpha, O2_PF_1, alpha, O2_PF_ML, 'LineWidth',1.5);

hold on 

% set the 0 to a different value so the patch command works:
if sum(isnan(O2_PF_max)) > 0
    O2_PF_max(isnan(O2_PF_max)) = min(O2_PF_max(O2_PF_max>0));
end

patch([alpha fliplr(alpha)], [O2_PF_min fliplr(O2_PF_max)], [241/255 218/255 170/255],  'FaceAlpha',.3, 'LineStyle','none')

xline(0.05/L, '--')

set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 15)
title('WIOT network', 'interpreter', 'latex', 'fontsize', 20)
xlabel("Bonferroni-corrected $\alpha$", 'interpreter', 'latex', 'FontSize',20)
ylabel('O$_2$', 'interpreter', 'latex', 'rotation', 0, 'verticalalignment', 'middle', 'fontsize', 30, 'horizontalalignment', 'right')
legend('DF', 'HF', 'PF$_{max}$', 'PF$_{min}$', 'PF$_1$', ...
        ' PF$_{ML}$', '', '', 'Location','southoutside', ...
        'orientation', 'horizontal', 'interpreter', 'latex')
xlim([1e-8 1e-1])
savefig(gcf, 'PLOT/O2_WIOT.fig')

hold off




%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [j] = jacc_weighted(a,b)

    j = sum(min(a,b)) / sum(max(a, b));
end

function [J] = Jaccard_weighted(W,B)
% this function compute the weighted Jaccard similarity between two network 
% whose adiacency matrices are W and B (ideally W is the original network 
% and B is the backbone filtered. 
% For each node, it compares the weights of outgoing links of the node in W
% and in B, and then averages over all the nodes that has at least one 
% outgoing link. 
% W and B are sparse objects

    N = length(W); % total number of nodes

    % compute the real number of nodes connected in the network
    [row,col] = find(W>0);
    %new = [row col];

    % list of indeces of nodes which has at least an outgoing link
    ind = unique(row);  
    N_A = length(ind);

    % initialize vector of jaccard similarities for each node 
    jacc = [];

    % loop over connected nodes
    for i = 1:N_A

        jacc(i) = sum(min( W(ind(i),:), B(ind(i),:) )) / sum( max(W(ind(i),:), B(ind(i),:)) );
    end
    % compute the average
    J = mean(jacc);

end