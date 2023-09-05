clear all
close all

addpath("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/PLOT/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/PLOTPLOT/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/MAX&SAM/MAX&SAM-routines/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/3_BA_NET/")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FAKE LINKS GROUND TRUTH TEST 1
%%% This file aims to perform a ground truth analysis on EMPIRICAL networks
%%% Given a statistical significance level it filters an untouched network,
%%% then it choses a random fraction of nodes that aren't connected to each
%%% other, create a link between them with weight gamma * median, gamma is
%%% the scaling factor, and then filter again the network. 
%%% at the end it calculate:
%%% whether or not if the new links have been validated. 
%%% and the proportion with a control group 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% 1) LOAD THE NETWORK AND FILTER IT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % load the original network (airports)
% A = matfile("saveA.mat"); 
% W = A.A;
 
% load the WIOT network 
W = readmatrix('/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/WIOT.txt');


l_ind   = find(W>0); % list of all links
w_list  = W(l_ind);  % list of all weights in the network
L       = nnz(W);    % count the number of links
N       = length(W); % count number of nodes
alpha   = 0.05 / L;  % multivariate significance level 

ap1     = 0.1;       % minimum Polya parameter 
ap2     = 4;         % maximum Polya parameter 
apr_lvl = 1e23;      % approximation level of polya filter 

% % extract backbone and p-values list 
% % [bb_ME, p_ME]          = max_filter(A, 'solxy.mat', alpha);
% [bb_DF, p_DF]            = disp_filter(W, alpha);
% [bb_HF, p_HF]            = hypergeom_filter(W, alpha);
% [bb_PF_min, p_PF_min]    = PF(W, ap1, alpha, apr_lvl, 0);
% [bb_PF_max, p_PF_max]    = PF(W, ap2, alpha, apr_lvl, 0);
% [bb_PF_1, p_PF_1]        = PF(W, 1, alpha, apr_lvl, 0);
[bb_PF_ML, p_PF_ML, aml] = PF(W, -1, alpha, apr_lvl, 0);


%%% 2) START THE BOOTSTRAP ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nbs   = 1;                     % number opf resampling in bootstrap analysis
frac  = 0.01;                  % fraction of link manipulated
nlink = round(frac*L);         % number of link to be added
gamma = [1 5 10 50 100, 500];  % scaling factors for manipulated weights


% initialize vectors of True Positive (link addend and validated)
TP_DF     = []; TP_HF   = []; TP_PF_min = [];
TP_PF_max = []; TP_PF_1 = []; TP_PF_ML  = [];

% initialize vectore of False Negative (link added and NOT validated)
FN_DF     = []; FN_HF   = []; FN_PF_min = [];
FN_PF_max = []; FN_PF_1 = []; FN_PF_ML  = [];

% initialize vectore of links in the control group validated 
CON_DF     = []; CON_HF   = []; CON_PF_min = [];
CON_PF_max = []; CON_PF_1 = []; CON_PF_ML  = [];

% initialize vectore of goodnesses
goodness_DF     = []; goodness_HF   = []; goodness_PF_min = [];
goodness_PF_max = []; goodness_PF_1 = []; goodness_PF_ML  = [];

% initialize vector of length of backbones extracted
lb_DF     = []; lb_HF   = []; lb_PF_min = [];
lb_PF_max = []; lb_PF_1 = []; lb_PF_ML  = [];


%%% 3) SELECT A FRACTION OF NODES AND ADD A LINK BETWEEN THEM %%%%%%%%%%%%%


for i=1:Nbs
    % restore the weights matrix
    W = readmatrix("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/WIOT.txt");
%     W = A.A;

   
    %%% GENERATE THE SET OF NEW LINKS

    % find nodes from whom no links are starting
    s   = find(all(W == 0, 2));
    % find nodes to whom no link are arriving 
    t = find(all(W==0));

    % create a matrix of the same size of W
    entries = ones(N,N);

    % set to 0 every row and col without pre-existing connections 
    entries(s, :) = 0;
    entries(:, t) = 0; 

    % find indeces of link in the net 
    [r,c] = ind2sub([N, N], l_ind);
    % set to 0 their corrispondent in the new matrix
    for j = 1:length(r)
        entries(r(j),c(j)) = 0;
    end


    % extract nlink random entries in order to have node with other connection 
    % but without the one that we're extracting to create new links   
    ind = round(randsample(find(entries==1), nlink));

    % extract nlink weights from the weight list
%     r_weight = randsample(w_list, nlink);
    r_weight = median(w_list) * ones(nlink,1);



    for j = 1:length(gamma)

        disp(j)

        tic
        % assign the weights 
        W(ind) = r_weight * gamma(j);


        % GENERATE THE CONTROL GROUP (LINKS WITH A WEIGHT ON AVERAGE 
        % SIMILAR TO GAMMA * MEDIAN

        sigma = gamma(j) * median(w_list) * 2/3;
        med   = gamma(j) * median(w_list);

        c_ind1 = find(w_list < med + sigma);
        c_ind2 = find(w_list > med - sigma);
        c_ind = intersect(c_ind1, c_ind2);
        %c_ind = c_ind(1:nlink);
        c_ind = randsample(c_ind, nlink);

        % weight matrix indeces of control group 
        control = l_ind(c_ind);


        %%% 4) FILTER AGAIN THE NETWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % [bb_ME2, p_ME2]         = max_filter(W, 'solxy.mat', alpha);
        [bb_DF2, p_DF2]         = disp_filter(W, alpha);
        lb_DF(i,j) = length(bb_DF2); 
        disp('end disparity filter')

        [bb_HF2, p_HF2]         = hypergeom_filter(W, alpha);
        disp('end hypergeometric filter')

        [bb_PF_min2, p_PF_min2] = PF(W, ap1, alpha, apr_lvl, 0);
        lb_PF_min(i,j) = length(bb_PF_min2);
        disp('end polya filter a=0.1')

        [bb_PF_max2, p_PF_max2] = PF(W, ap2, alpha, apr_lvl, 0);
        lb_PF_max(i,j) = length(bb_PF_max2);
        disp('end polya filter a=4')

        [bb_PF_12, p_PF_12]     = PF(W, 1, alpha, apr_lvl, 0);
        lb_PF_1(i,j) = length(bb_PF_12);
        disp('end polya filter a=1')
        %%% ATTENZIONE - ADESSO C'E' A = ALL'A TROVATO DAL ML DELLA RETE NON
        %%% MODIFICATA - PER FAR RICERCARE A METTERE - 1 AL SUO POSTO
        [bb_PF_ML2, p_PF_ML2]   = PF(W, aml, alpha, apr_lvl, 0);
        lb_PF_ML(i,j) = length(bb_PF_ML2);
        disp('end polya filter a=4.5')

        %%% PROTECT THE EMPTY BACKBONES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(bb_DF2) == 0
            b_DF2     = sub2ind(size(W), bb_DF2(:,1), bb_DF2(:,2));
        else
            b_DF2 = 0;
        end
        if isempty(bb_HF2) == 0
            b_HF2     = sub2ind(size(W), bb_HF2(:,1), bb_HF2(:,2));
        else
            b_HF2 = 0;
       end
        if isempty(bb_PF_min2) == 0
            b_PF_min2 = sub2ind(size(W), bb_PF_min2(:,1), bb_PF_min2(:,2));
        else
            b_PF_min2 = 0;
        end
        if isempty(bb_PF_max2) == 0
            b_PF_max2 = sub2ind(size(W), bb_PF_max2(:,1), bb_PF_max2(:,2));
        else
            b_PF_max2 = 0;
        end
        if isempty(bb_PF_12) == 0
            b_PF_12   = sub2ind(size(W), bb_PF_12(:,1), bb_PF_12(:,2));
        else
            b_PF_12 = 0;
        end
        if isempty(bb_PF_ML2) == 0
            b_PF_ML2  = sub2ind(size(W), bb_PF_ML2(:,1), bb_PF_ML2(:,2));
        end

        %%% 5) GET THE PARATEMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % get the fraction of true positive (manipulated and validated) 
        TP_DF(i,j)     = sum(ismember(b_DF2, ind)) / nlink;
        TP_HF(i,j)     = sum(ismember(b_HF2, ind)) / nlink;
        TP_PF_min(i,j) = sum(ismember(b_PF_min2, ind)) / nlink;
        TP_PF_max(i,j) = sum(ismember(b_PF_max2, ind)) / nlink;
        TP_PF_1(i,j)   = sum(ismember(b_PF_12, ind)) / nlink;
        TP_PF_ML(i,j)  = sum(ismember(b_PF_ML2, ind)) / nlink;
    
        % get the fraction of false negative (manipulated NON validated)
        FN_DF(i,j)     = 1 - TP_DF(i,j); 
        FN_HF(i,j)    = 1 - TP_HF(i,j); 
        FN_PF_min(i,j) = 1 - TP_PF_min(i,j); 
        FN_PF_max(i,j) = 1 - TP_PF_max(i,j); 
        FN_PF_1(i,j)   = 1 - TP_PF_1(i,j); 
        FN_PF_ML(i,j)   = 1 - TP_PF_ML(i,j); 

        % get the fraction of links in the control group validated 
        CON_DF(i,j)     = sum(ismember(b_DF2, c_ind)) / nlink;
        CON_HF(i,j)     = sum(ismember(b_HF2, c_ind)) / nlink;
        CON_PF_min(i,j) = sum(ismember(b_PF_min2, c_ind)) / nlink;
        CON_PF_max(i,j) = sum(ismember(b_PF_max2, c_ind)) / nlink;
        CON_PF_1(i,j)   = sum(ismember(b_PF_12, c_ind)) / nlink;
        CON_PF_ML(i,j)  = sum(ismember(b_PF_ML2, c_ind)) / nlink;
    
    
        % goodness measure: distance from "perfect point" (TP=1, FP=0)
        % più è piccola e più il filtro è vicino al filtro perfetto 
        
        goodness_DF(i,j)     = (1-TP_DF(i,j))^2     + FN_DF(i,j)^2;    
        goodness_HF(i,j)     = (1-TP_HF(i,j))^2     + FN_HF(i,j)^2;  
        goodness_PF_min(i,j) = (1-TP_PF_min(i,j))^2 + FN_PF_min(i,j)^2;
        goodness_PF_max(i,j) = (1-TP_PF_max(i,j))^2 + FN_PF_max(i,j)^2;
        goodness_PF_1(i,j)   = (1-TP_PF_1(i,j))^2   + FN_PF_1(i,j)^2;
        goodness_PF_ML(i,j)  = (1-TP_PF_ML(i,j))^2  + FN_PF_ML(i,j)^2;
    
        toc
    end

end
% disp(['DF goodness = ', num2str(mean(goodness_DF))])
% disp(['HF goodness =',num2str(mean(goodness_HF))])
% disp(['PF min goodness = ', num2str(mean(goodness_PF_min))])
% disp(['PF max goodness = ', num2str(mean(goodness_PF_max))])
% disp(['PF 1 goodness = ', num2str(mean(goodness_PF_1))])
% disp(['PF ML goodness = ', num2str(mean(goodness_PF_ML))])

av_good = [mean(goodness_DF), mean(goodness_HF), mean(goodness_PF_min), mean(goodness_PF_max), mean(goodness_PF_1), mean(goodness_PF_ML)];

m = 1/length(b_DF2) + 1/length(b_PF_min2) + 1/length(b_PF_max2) + 1/length(b_PF_12)+ 1/length(b_PF_ML2);


%%% 6) PLOT THE ROC CURVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
semilogx(gamma, mean(TP_DF,1)./mean(CON_DF,1), '.-', 'MarkerSize', 15, 'linewidth', 1.3)
hold on
semilogx(gamma, mean(TP_HF,1)./mean(CON_HF,1), '.-', 'MarkerSize', 15, 'LineWidth',1.3)
semilogx(gamma, mean(TP_PF_min,1)./mean(CON_PF_min,1), '.-', 'MarkerSize', 15, 'LineWidth', 1.3)
semilogx(gamma, mean(TP_PF_max,1)./mean(CON_PF_max,1), '.-', 'MarkerSize', 15, 'LineWidth',1.3)
semilogx(gamma, mean(TP_PF_1,1)./mean(CON_PF_1,1), '.-', 'MarkerSize', 15, 'LineWidth',1.3)
semilogx(gamma, mean(TP_PF_ML,1)./mean(CON_PF_ML,1), '.-', 'MarkerSize', 15, 'LineWidth',1.3)

title('\textbf{WIOT network}', 'interpreter', 'latex', 'FontSize',12)
xlabel('$\mathbf{\gamma:}$ \textbf{scale factor of fake weigths}', 'Interpreter', 'latex', 'Fontsize', 12)
ylabel('$\mathbf{\frac{fake\> link\> validated}{control \> link \> validated}}$', 'Interpreter', 'latex', 'fontsize', 12)
set(gca, 'FontName', 'latex', 'TickLabelInterpreter', 'latex')
legend('DF', 'HF', 'PF $a_{MIN}$', 'PF $a_{MAX}$', 'PF $a=1$', 'PF $a_{ML}$', ...
    'interpreter', 'latex', 'location', 'southoutside', 'orientation', 'horizontal')

savefig(gcf, 'PLOT/TPcontrol-vs-gamma-WIOT.fig')


figure(2);

semilogx(gamma, mean(TP_DF,1), '.-', 'MarkerSize', 15, 'LineWidth',1.3)
hold on
semilogx(gamma, mean(TP_HF,1), '.-', 'MarkerSize', 15, 'LineWidth',1.3)
semilogx(gamma, mean(TP_PF_min,1), '.-', 'MarkerSize', 15, 'LineWidth',1.3)
semilogx(gamma, mean(TP_PF_max,1), '.-', 'MarkerSize', 15, 'LineWidth',1.3)
semilogx(gamma, mean(TP_PF_1,1), '.-', 'MarkerSize', 15, 'LineWidth',1.3)
semilogx(gamma, mean(TP_PF_ML,1), '.-', 'MarkerSize', 15, 'LineWidth',1.3)

title('\textbf{WIOT network}', 'interpreter', 'latex', 'fontsize', 12)
xlabel('$\mathbf{\gamma:}$ \textbf{scale factor of fake weigths}', 'Interpreter', 'latex', 'Fontsize', 12)
ylabel('\textbf{Fraction of fake links validated}', 'Interpreter', 'latex', 'fontsize', 12)
set(gca, 'FontName', 'latex', 'TickLabelInterpreter', 'latex', 'FontSize', 20)
legend('DF', 'HF', 'PF $a_{MIN}$', 'PF $a_{MAX}$', 'PF $a=1$', 'PF $a_{ML}$', ...
    'interpreter', 'latex', 'location', 'southoutside', 'orientation', 'horizontal')

savefig(gcf, 'PLOT/TP-vs-gamma-WIOT.fig')

% scatter(TP_DF./length(b_DF2)./m , FN_DF./length(b_DF2)./m, 'filled')
% hold on 
% scatter(TP_HF , FN_HF, 'filled')
% hold on 
% scatter(TP_PF_min./length(b_PF_min2)./m, FN_PF_min./length(b_PF_min2)./m, 'filled')
% hold on
% scatter(TP_PF_max./length(b_PF_max2)./m, FN_PF_max./length(b_PF_max2)./m, 'filled')
% hold on
% scatter(TP_PF_1./length(b_PF_12)./m,     FN_PF_1./length(b_PF_12)./m, 'filled')
% hold on
% scatter(TP_PF_ML./length(b_PF_ML2)./m,   FN_PF_ML./length(b_PF_ML2)./m, 'filled')
% hold on

% title('Airport network with a fraction of fake connections')
% xlabel('TP: fraction of fake links that are validated')
% ylabel('FP: fraction of fake links that are NOT validated')
% legend('DF', 'HF', 'PF MIN', 'PF MAX', 'PF 1', 'PF ML', 'location', 'northeast')
% savefig(gcf, 'PLOT/ROC_airport_weight.fig')

% figure(3)
% num = linspace(1, 6, 6);
% 
% semilogy(num, av_good, '*', 'MarkerSize', 20, 'color', '#6D69DF', 'LineWidth',2)
% 
% xticks(num)
% xticklabels({'DF', 'HF', 'PF MIN', 'PF MAX', 'PF 1', 'PF ML'})
% 
% xlabel('Filter')
% ylabel('Average goodness measure')
% title('Average goodness in detecting manipulated links in the airport network')
% 
% axis([0 7 1e-4 10]) % set axis limits - first 2 for the x second 2 for the y 
% savefig(gcf, 'PLOT/goodness_airport.fig')

