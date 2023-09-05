clear all
close all

addpath("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/PLOT/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/PLOTPLOT/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/MAX&SAM/MAX&SAM-routines/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/3_BA_NET/")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MANIPULATED LINKS GROUND TRUTH TEST 2
%%% This file aims to perform a ground truth analysis on synthetic or 
%%% EMPIRICAL networks
%%% Given a statistical significance level it filters an untouched network,
%%% then it choses a random fraction of links among the unfiltered ones, it
%%% manipulates their weights and then filter again the network. 
%%% at the end it calculate:
%%% the fraction of link manipulated retained in the backbone (True Positive), 
%%% the fraction of link NON manipulated NON retained (True Negative),
%%% the fraction of link manipulated NON retained (False Negative),
%%% the fraction of link NON manipulated retained (False Positive),

%%% The manipulation is made scaling the weight by a factor (1 + gamma)
%%% where gamma spans several order of magnitude

%%% it's already implememted the possibility to do a bootstrap analysis
%%% extracting multiple control sets. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% 1) LOAD THE NETWORK AND FILTER IT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % load the ER weighted network with power-law weight distribution
% ER = load('saveweightedER_powelaw.mat');
% W = ER.W;

% % load the BA weighted directed network 
% % (weights extracted from power law distribution correlated in first order to degree)
% A = load("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/saveweightedBA_coupled.mat");
% W = A.W;

% % load the BA weighted directed network 
% % (weights extracted from power law distribution)
% A = load("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/saveweightedBA.mat");
% W = A.W;

% % load the original network (airports)
% A = matfile("saveA.mat"); 
% W = A.A;
 
% load the WIOT network 
W = readmatrix('/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/WIOT.txt');



l_ind   = find(W>0); % list of all links
L       = nnz(W);    % count the number of links
alpha   = 0.05 / L;  % multivariate significance level 

ap1     = 0.1;       % minimum Polya parameter 
ap2     = 4;         % maximum Polya parameter 
apr_lvl = 1e23;      % approximation level of polya filter 

% extract backbone and p-values list 
% [bb_ME, p_ME]          = max_filter(A, 'solxy.mat', alpha);
[bb_DF, p_DF]            = disp_filter(W, alpha);
[bb_HF, p_HF]            = hypergeom_filter(W, alpha);
[bb_PF_min, p_PF_min]    = PF(W, ap1, alpha, apr_lvl, 0);
[bb_PF_max, p_PF_max]    = PF(W, ap2, alpha, apr_lvl, 0);
[bb_PF_1, p_PF_1]        = PF(W, 1, alpha, apr_lvl, 0);
[bb_PF_ML, p_PF_ML, aml] = PF(W, -1, alpha, apr_lvl, 0);

% convert edge list in linear index if isempty(bb_DF) == 0
    if isempty(bb_DF) == 0
        b_DF     = sub2ind(size(W), bb_DF(:,1), bb_DF(:,2));
    else
        b_DF = 0;
    end
    if isempty(bb_HF) == 0
        b_HF     = sub2ind(size(W), bb_HF(:,1), bb_HF(:,2));
    else
        b_HF = 0;
    end
    if isempty(bb_PF_min) == 0
        b_PF_min = sub2ind(size(W), bb_PF_min(:,1), bb_PF_min(:,2));
    else
        b_PF_min = 0;
    end
    if isempty(bb_PF_max) == 0
        b_PF_max = sub2ind(size(W), bb_PF_max(:,1), bb_PF_max(:,2));
    else
        b_PF_max = 0;
    end
    if isempty(bb_PF_1) == 0
        b_PF_1   = sub2ind(size(W), bb_PF_1(:,1), bb_PF_1(:,2));
    else
        b_PF_1 = 0;
    end
    if isempty(bb_PF_ML) == 0
        b_PF_ML  = sub2ind(size(W), bb_PF_ML(:,1), bb_PF_ML(:,2));
    else
        b_PF_ML  = 0;
    end



%%% 3) START THE BOOTSTRAP ANALYSIS OF TP, TN, FP, FN %%%%%%%%%%%%%%%%%%%%%
Nbs = 1;                % number opf resampling in bootstrap analysis
frac = 0.05;             % fraction of link manipulated
nlink = round(frac*L);   % number of link manipulated 

% initialize all vectors needed 
TP_DF     = []; TP_HF     = []; TP_PF_min = [];
TP_PF_max = []; TP_PF_1   = []; TP_PF_ML  = []; 

TN_DF     = []; TN_HF     = []; TN_PF_min = [];
TN_PF_max = []; TN_PF_1   = []; TN_PF_ML  = []; 

FP_DF     = []; FP_HF     = []; FP_PF_min = [];
FP_PF_max = []; FP_PF_1   = []; FP_PF_ML  = []; 

FN_DF     = []; FN_HF     = []; FN_PF_min = [];
FN_PF_max = []; FN_PF_1   = []; FN_PF_ML  = []; 

goodness_DF     = []; goodness_HF   = []; goodness_PF_min = [];
goodness_PF_max = []; goodness_PF_1 = []; goodness_PF_ML = [];



%%% 3.1) SELECT A FRACTION OF THEN AND MANIPULATE THEIR WEIGHTS %%%%%%%%%%%
% array of 
gamma = [ 0.1 0.5 1 5 10 50 ];


for j = 1:length(gamma)
    disp(j)

    % restore the weights matrix
    W = readmatrix("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/WIOT.txt");
    % W = A.A;

    %%% DISPARITY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % take off the validated links from the list
    antib = l_ind;
    antib(ismember(antib, b_DF)) = [];

    % extract nlink random links out of the list of non validated ones 
    ind = randsample(antib, nlink);

    % take off also manipulated links from the list 
    antib2 = antib;
    antib2(ismember(antib, ind)) =[];
    
    % assign the weights
    W(ind) = W(ind) * (1 + gamma(j));    


    for i=1:Nbs

        % extract the control list, of non validated and non modified   
        ind2 = randsample(antib2, nlink);
    
        % filter the network
        [bb_DF2, p_DF2] = disp_filter(W, alpha);

        if isempty(bb_DF2) == 0
            b_DF2     = sub2ind(size(W), bb_DF2(:,1), bb_DF2(:,2));
        else
            b_DF2 = 0;
        end

        % get the fraction of TP - FN - FP - TN 
        TP_DF(j,i)     = sum(ismember(b_DF2, ind)) / nlink;
        FN_DF(j,i)     = 1 - TP_DF(i); 
        FP_DF(j,i)     = sum(ismember(b_DF2, ind2)) / nlink;
        TN_DF(j,i)     = 1 - FP_DF(i); 

    end

    %%% HYPERGEOMETRIC FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % restore the weights matrix
    W = readmatrix("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/WIOT.txt");
    % W = A.A;

    % take off the validated links from the list
    antib = l_ind;
    antib(ismember(antib, b_HF)) = [];

    % extract nlink random links out of the list of non validated ones 
    ind = randsample(antib, nlink);

    % take off also manipulated links from the list 
    antib2 = antib;
    antib2(ismember(antib, ind)) =[];
    
    % assign the weights
    W(ind) = W(ind) * (1 + gamma(j));    

    for i=1:Nbs
    
       % extract the control list, of non validated and non modified 
       ind2 = randsample(antib2, nlink);
    
       % filter the network
       [bb_HF2, p_HF2] = hypergeom_filter(W, alpha);
    
       if isempty(bb_HF2) == 0
           b_HF2     = sub2ind(size(W), bb_HF2(:,1), bb_HF2(:,2));
       else
           b_HF2 = 0;
       end
    
       % get the fraction of TP - FN - FP - TN 
       TP_HF(j,i)     = sum(ismember(b_HF2, ind)) / nlink;
       FN_HF(j,i)     = 1 - TP_HF(i); 
       FP_HF(j,i)     = sum(ismember(b_HF2, ind2)) / nlink;
       TN_HF(j,i)     = 1 - FP_HF(i); 
    
    end

    %%% POLYA FILTER A MIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % restore the weights matrix
    W = readmatrix("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/WIOT.txt");
    % W = A.A;

    % take off the validated links from the list
    antib = l_ind;
    antib(ismember(antib, b_PF_min)) = [];

    % extract nlink random links out of the list of non validated ones 
    ind = randsample(antib, nlink);

    % take off also manipulated links from the list 
    antib2 = antib;
    antib2(ismember(antib, ind)) =[];
    
    % assign the weights
    W(ind) = W(ind) * (1 + gamma(j));    

     for i=1:Nbs

        % extract the control list, of non validated and non modified 
        ind2 = randsample(antib2, nlink);
    
        % filter the network
        [bb_PF_min2, p_PF_min2] = PF(W, ap1, alpha, apr_lvl, 0);

        b_PF_min2     = sub2ind(size(W), bb_PF_min2(:,1), bb_PF_min2(:,2));

        % get the fraction of TP - FN - FP - TN 
        TP_PF_min(j,i)     = sum(ismember(b_PF_min2, ind)) / nlink;
        FN_PF_min(j,i)     = 1 - TP_PF_min(i); 
        FP_PF_min(j,i)     = sum(ismember(b_PF_min2, ind2)) / nlink;
        TN_PF_min(j,i)     = 1 - FN_PF_min(i); 

     end

     %%% POLYA FILTER A MAX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % restore the weights matrix
    W = readmatrix("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/WIOT.txt");
    % W = A.A;

    % take off the validated links from the list
    antib = l_ind;
    antib(ismember(antib, b_PF_max)) = [];

    % extract nlink random links out of the list of non validated ones 
    ind = randsample(antib, nlink);

    % take off also manipulated links from the list 
    antib2 = antib;
    antib2(ismember(antib, ind)) =[];
    
    % assign the weights
    W(ind) = W(ind) * (1 + gamma(j));    

     for i=1:Nbs

        % extract the control list, of non validated and non modified 
        ind2 = randsample(antib2, nlink);
    
        % filter the network
        [bb_PF_max2, p_PF_max2] = PF(W, ap2, alpha, apr_lvl, 0);

        b_PF_max2     = sub2ind(size(W), bb_PF_max2(:,1), bb_PF_max2(:,2));

        % get the fraction of TP - FN - FP - TN 
        TP_PF_max(j,i)     = sum(ismember(b_PF_max2, ind)) / nlink;
        FN_PF_max(j,i)     = 1 - TP_PF_max(i); 
        FP_PF_max(j,i)     = sum(ismember(b_PF_max2, ind2)) / nlink;
        TN_PF_max(j,i)     = 1 - FN_PF_max(i); 

     end

     %%% POLYA FILTER A = 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % restore the weights matrix
    W = readmatrix("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/WIOT.txt");
    % W = A.A;

    % take off the validated links from the list
    antib = l_ind;
    antib(ismember(antib, b_PF_1)) = [];

    % extract nlink random links out of the list of non validated ones 
    ind = randsample(antib, nlink);

    % take off also manipulated links from the list 
    antib2 = antib;
    antib2(ismember(antib, ind)) =[];
    
    % assign the weights
    W(ind) = W(ind) * (1 + gamma(j));    

     for i=1:Nbs

        % extract the control list, of non validated and non modified 
        ind2 = randsample(antib2, nlink);
    
        % filter the network
        [bb_PF_12, p_PF_12] = PF(W, 1, alpha, apr_lvl, 0);

        b_PF_12     = sub2ind(size(W), bb_PF_12(:,1), bb_PF_12(:,2));

        % get the fraction of TP - FN - FP - TN 
        TP_PF_1(j,i)     = sum(ismember(b_PF_12, ind)) / nlink;
        FN_PF_1(j,i)     = 1 - TP_PF_1(i); 
        FP_PF_1(j,i)     = sum(ismember(b_PF_12, ind2)) / nlink;
        TN_PF_1(j,i)     = 1 - FN_PF_1(i); 

     end

     %%% POLYA FILTER A MLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % restore the weights matrix
    W = readmatrix("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/WIOT.txt");
    % W = A.A;

    % take off the validated links from the list
    antib = l_ind;
    antib(ismember(antib, b_PF_ML)) = [];

    % extract nlink random links out of the list of non validated ones 
    ind = randsample(antib, nlink);

    % take off also manipulated links from the list 
    antib2 = antib;
    antib2(ismember(antib, ind)) =[];
    
    % assign the weights
    W(ind) = W(ind) * (1 + gamma(j));    

     for i=1:Nbs

        % extract the control list, of non validated and non modified 
        ind2 = randsample(antib2, nlink);
    
        % filter the network
        [bb_PF_ML2, p_PF_ML2] = PF(W, aml, alpha, apr_lvl, 0);

        b_PF_ML2     = sub2ind(size(W), bb_PF_ML2(:,1), bb_PF_ML2(:,2));

        % get the fraction of TP - FN - FP - TN 
        TP_PF_ML(j,i)     = sum(ismember(b_PF_ML2, ind)) / nlink;
        FN_PF_ML(j,i)     = 1 - TP_PF_ML(i); 
        FP_PF_ML(j,i)     = sum(ismember(b_PF_ML2, ind2)) / nlink;
        TN_PF_ML(j,i)     = 1 - FN_PF_ML(i); 

    end

        % goodness measure: distance from "perfect point" (TP=1, FP=0)
        % più è piccola e più il filtro è vicino al filtro perfetto 
        
%         goodness_DF(j,i)     = (1-TP_DF(j,i))^2     + FP_DF(j,i)^2     + (1-TN_DF(j,i))^2     + FN_DF(j,i)^2 ;
%         goodness_HF(j,i)     = (1-TP_HF(j,i))^2     + FP_HF(j,i)^2     + (1-TN_HF(j,i))^2     + FN_HF(j,i)^2;
%         goodness_PF_min(j,i) = (1-TP_PF_min(j,i))^2 + FP_PF_min(j,i)^2 + (1-TN_PF_min(j,i))^2 + FN_PF_min(j,i)^2 ;
%         goodness_PF_max(j,i) = (1-TP_PF_max(j,i))^2 + FP_PF_max(j,i)^2 + (1-TN_PF_max(j,i))^2 + FN_PF_max(j,i)^2;
%         goodness_PF_1(j,i)   = (1-TP_PF_1(j,i))^2   + FP_PF_1(j,i)^2   + (1-TN_PF_1(j,i))^2   + FN_PF_1(j,i)^2;
%         goodness_PF_ML(j,i)  = (1-TP_PF_ML(j,i))^2  + FP_PF_ML(j,i)^2  + (1-TN_PF_ML(j,i))^2  + FN_PF_ML(j,i)^2;
    
end

% disp(['DF goodness = ', num2str(mean(goodness_DF,2 ))])
% disp(['HF goodness =',num2str(mean(goodness_HF,2))])
% disp(['PF min goodness = ', num2str(mean(goodness_PF_min,2))])
% disp(['PF max goodness = ', num2str(mean(goodness_PF_max,2))])
% disp(['PF 1 goodness = ', num2str(mean(goodness_PF_1,2))])
% disp(['PF ML goodness = ', num2str(mean(goodness_PF_ML,2))])
% 
% av_good = [mean(goodness_DF,2), mean(goodness_HF,2), mean(goodness_PF_min,2), mean(goodness_PF_max,2), mean(goodness_PF_1,2), mean(goodness_PF_ML,2)];

%%% 5) PLOT THE ROC CURVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);

plot(mean(TP_DF,2) , mean(FP_DF,2), '.-', 'MarkerSize',10, 'LineWidth',2)
plot(mean(TP_HF,2) , mean(FP_HF,2), '.-', 'MarkerSize',10, 'LineWidth',1.8)
plot(mean(TP_PF_min,2) , mean(FP_PF_min,2), '.-', 'MarkerSize',10, 'LineWidth',1.6)
plot(mean(TP_PF_max,2) , mean(FP_PF_max,2), '.-', 'MarkerSize',10, 'LineWidth',1.4)
plot(mean(TP_PF_1,2) , mean(FP_PF_1,2), '.-', 'MarkerSize',10, 'LineWidth',1.2)
plot(mean(TP_PF_ML,2) , mean(FP_PF_ML,2), '.-', 'MarkerSize',10, 'LineWidth',1)

title('\textbf{US air traffic network}', 'interpreter', 'latex', 'FontSize',12)
xlabel('\textbf{TP: manipulated and validated links}', 'Interpreter', 'latex', 'Fontsize', 12)
ylabel('\textbf{FP: NON manipulated and validated links}', 'Interpreter', 'latex', 'fontsize', 12)
set(gca, 'FontName', 'latex', 'TickLabelInterpreter', 'latex')
legend('DF', 'HF', 'PF $a_{MIN}$', 'PF $a_{MAX}$', 'PF $a=1$', 'PF $a_{ML}$', ...
    'interpreter', 'latex', 'location', 'southoutside', 'orientation', 'horizontal')

savefig(gcf, 'PLOT/ROC-WIOT.fig')


figure(2);

loglog(gamma, TP_DF,    '.-', 'MarkerSize', 15, 'LineWidth',1.3)
hold on
loglog(gamma, TP_HF,    '.-', 'MarkerSize', 15, 'LineWidth',1.3)
loglog(gamma, TP_PF_min,'.-', 'MarkerSize', 15, 'LineWidth',1.3)
loglog(gamma, TP_PF_max,'.-', 'MarkerSize', 15, 'LineWidth',1.3)
loglog(gamma, TP_PF_1,  '.-', 'MarkerSize', 15, 'LineWidth',1.3)
loglog(gamma, TP_PF_ML, '.-', 'MarkerSize', 15, 'LineWidth',1.3)

title('\textbf{WIOT network with a fraction of fake connections}', 'interpreter', 'latex', 'fontsize', 12)
xlabel('$\mathbf{\gamma:}$ \textbf{scale factor of fake weigths}', 'Interpreter', 'latex', 'Fontsize', 12)
ylabel('\textbf{TP: Fraction of manipulated links validated}', 'Interpreter', 'latex', 'fontsize', 12)
set(gca, 'FontName', 'latex', 'TickLabelInterpreter', 'latex')
legend('DF', 'HF', 'PF $a_{MIN}$', 'PF $a_{MAX}$', 'PF $a=1$', 'PF $a_{ML}$', ...
    'interpreter', 'latex', 'location', 'southoutside', 'orientation', 'horizontal')

savefig(gcf, 'PLOT/TP-vs-gamma-WIOT2.fig')


figure(3);

loglog(gamma, TP_DF./length(b_DF).*L,        '.-', 'MarkerSize', 15, 'LineWidth',1.3)
hold on
loglog(gamma, TP_HF./length(b_HF).*L,        '.-', 'MarkerSize', 15, 'LineWidth',1.3)
loglog(gamma, TP_PF_min./length(b_PF_min).*L,'.-', 'MarkerSize', 15, 'LineWidth',1.3)
loglog(gamma, TP_PF_max./length(b_PF_max).*L,'.-', 'MarkerSize', 15, 'LineWidth',1.3)
loglog(gamma, TP_PF_1./length(b_PF_1).*L,    '.-', 'MarkerSize', 15, 'LineWidth',1.3)
loglog(gamma, TP_PF_ML./length(b_PF_ML).*L,  '.-', 'MarkerSize', 15, 'LineWidth',1.3)

title('\textbf{WIOT network with a fraction of fake connections}', 'interpreter', 'latex', 'fontsize', 12)
xlabel('$\mathbf{\gamma:}$ \textbf{scale factor of fake weigths}', 'Interpreter', 'latex', 'Fontsize', 12)
ylabel('$\mathbf{\frac{validation\> rate\> manipulate \> links}{total \> validation \> rate}}$', 'Interpreter', 'latex', 'fontsize', 12)
set(gca, 'FontName', 'latex', 'TickLabelInterpreter', 'latex')
legend('DF', 'HF', 'PF $a_{MIN}$', 'PF $a_{MAX}$', 'PF $a=1$', 'PF $a_{ML}$', ...
    'interpreter', 'latex', 'location', 'southoutside', 'orientation', 'horizontal')

savefig(gcf, 'PLOT/TPweigh-vs-gamma-WIOT2.fig')




figure(3)
num = linspace(1,length(gamma), length(gamma));
num = [num; num; num; num; num; num];
num = num';
semilogy(num, av_good, '*', 'MarkerSize', 20, 'color', '#6D69DF', 'LineWidth',2)

xticks(linspace(1,6,6))
xticklabels({'DF', 'HF', 'PF MIN', 'PF MAX', 'PF 1', 'PF ML'})

xlabel('Filter')
ylabel('Average goodness measure')
title('Average goodness in detecting manipulated links in a BA network')

axis([0 7 1e-4 10]) % set axis limits - first 2 for the x second 2 for the y 
savefig(gcf, 'PLOT/goodness_BA_pl-gamma.fig')

