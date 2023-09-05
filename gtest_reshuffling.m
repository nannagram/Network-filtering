clear all
close all

addpath("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/PLOT/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/PLOTPLOT/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/MAX&SAM/MAX&SAM-routines/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/3_BA_NET/")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RESHUFFLING TEST
%%% This file aims to perform a test on EMPIRICAL network. 
%%% It filters the untouched network and then reshuffle the weigth assigned
%%% to each link, subsequentetly it filters again the network and see if
%%% the rate of validation changes. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% 1) LOAD THE NETWORK AND FILTER IT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % load the original network (airports)
% A = matfile("saveA.mat"); 
% W = A.A;
 
% load the WIOT network 
W = readmatrix('WIOT.txt');


L       = nnz(W);    % count the number of links
N       = length(W); % count number of nodes
l_ind   = find(W>0); % list of all links

w_list  = W(l_ind);  % list of all weights in the network

alpha   = 0.05 / L;  % multivariate significance level 

ap1     = 0.1;       % minimum Polya parameter 
ap2     = 4;         % maximum Polya parameter 
apr_lvl = 1e23;      % approximation level of polya filter 


k_in  = full(sum(W>0)); % Degree sequence (in-degree) (sum over column)
s_in  = full(sum(W)); % Stregth degree (in-stregth) 
k_out = full(sum(W'>0)); % Degree sequence (in-degree) (sum over column)
s_out = full(sum(W')); % Stregth degree (in-stregth) 



%%% 1) FILTER THE UNTOUCHED NETWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [bb_ME, p_ME]          = max_filter(A, 'solxy.mat', alpha);
[bb_DF, p_DF]            = disp_filter(W, alpha);
disp('end Disparity Filter')
[bb_HF, p_HF]            = hypergeom_filter(W, alpha);
disp('end Hypergeometric Filter')
[bb_PF_min, p_PF_min]    = PF(W, ap1, alpha, apr_lvl, 0);
disp('end Polya Filter a = 0.1')
[bb_PF_max, p_PF_max]    = PF(W, ap2, alpha, apr_lvl, 0);
disp('end Polya Filter a = 4')
[bb_PF_1, p_PF_1]        = PF(W, 1, alpha, apr_lvl, 0);
disp('end Polya Filter a = 1')
[bb_PF_ML, p_PF_ML, aml] = PF(W, -1, alpha, apr_lvl, 0);
disp('end Polya Filter a ML')

fr_link = [length(bb_DF)/L, length(bb_HF)/L, length(bb_PF_min)/L, length(bb_PF_max)/L, length(bb_PF_1)/L, length(bb_PF_ML)/L];

%%% 2) RESHUFFLE THE NETWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsims = 1;
fr_link2 = [];

for i= 1:Nsims
    disp(['loop number ' i])
    W_shuf = W;
    w_shuf  = w_list(randperm(L)); % list of reshuffled weigths
    W_shuf(l_ind) = w_shuf;        % assign the reshuffled weight 

    k_in2  = full(sum(W_shuf>0)); % Degree sequence (in-degree) (sum over column)
    s_in2  = full(sum(W_shuf)); % Stregth degree (in-stregth) 
    k_out2 = full(sum(W_shuf'>0)); % Degree sequence (in-degree) (sum over column)
    s_out2 = full(sum(W_shuf')); % Stregth degree (in-stregth) 


%%% 2) FILTER THE RESHUFFLED NETWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % [bb_ME2, p_ME2]         = max_filter(W, 'solxy.mat', alpha);
    [bb_DF2, p_DF2]         = disp_filter(W_shuf, alpha);
    fr_link2(i,1) = length(bb_DF2)/L;
    disp('end disparity filter')
    
    [bb_HF2, p_HF2]         = hypergeom_filter(W_shuf, alpha);
    fr_link2(i,2) = length(bb_HF2)/L;
    disp('end hypergeometric filter')
    
    [bb_PF_min2, p_PF_min2] = PF(W_shuf, ap1, alpha, apr_lvl, 0);
    fr_link2(i,3) = length(bb_PF_min2)/L;
    disp('end polya filter a=0.1')
    
    [bb_PF_max2, p_PF_max2] = PF(W_shuf, ap2, alpha, apr_lvl, 0);
    fr_link2(i,4) = length(bb_PF_max2)/L;
    disp('end polya filter a=4')
    
    [bb_PF_12, p_PF_12]     = PF(W_shuf, 1, alpha, apr_lvl, 0);
    fr_link2(i,5) = length(bb_PF_12)/L;
    disp('end polya filter a=1')
    
    [bb_PF_ML2, p_PF_ML2, aml]   = PF(W_shuf, aml, alpha, apr_lvl, 0);
    fr_link2(i,6) = length(bb_PF_ML2)/L;
    disp('end polya filter a_{ML}')
    
end

%%% 3) FIND  JACCARD SIMILARITIES IN THE 2 BACKBONES %%%%%%%%%%%%%%%%%%%%%%

% convert edge lists in single indeces lists 
b_DF     = sub2ind(size(W), bb_DF(:,1), bb_DF(:,2));
b_HF     = sub2ind(size(W), bb_HF(:,1), bb_HF(:,2));
b_PF_min = sub2ind(size(W), bb_PF_min(:,1), bb_PF_min(:,2));
b_PF_max = sub2ind(size(W), bb_PF_max(:,1), bb_PF_max(:,2));
b_PF_1   = sub2ind(size(W), bb_PF_1(:,1), bb_PF_1(:,2));
b_PF_ML  = sub2ind(size(W), bb_PF_ML(:,1), bb_PF_ML(:,2));

b_DF2     = sub2ind(size(W), bb_DF2(:,1), bb_DF2(:,2));
b_HF2     = sub2ind(size(W), bb_HF2(:,1), bb_HF2(:,2));
b_PF_min2 = sub2ind(size(W), bb_PF_min2(:,1), bb_PF_min2(:,2));
b_PF_max2 = sub2ind(size(W), bb_PF_max2(:,1), bb_PF_max2(:,2));
b_PF_12   = sub2ind(size(W), bb_PF_12(:,1), bb_PF_12(:,2));
b_PF_ML2  = sub2ind(size(W), bb_PF_ML2(:,1), bb_PF_ML2(:,2));

wint_DF = intersect(b_DF, b_DF2);
wun_DF = union(b_DF, b_DF2);
jacc_DF = length(wint_DF) / length(wun_DF);

wint_HF = intersect(b_HF, b_HF2);
wun_HF = union(b_HF, b_HF2);
jacc_HF = length(wint_HF) / length(wun_HF);

wint_PF_min = intersect(b_PF_min, b_PF_min2);
wun_PF_min = union(b_PF_min, b_PF_min2);
jacc_PF_min = length(wint_PF_min) / length(wun_PF_min);

wint_PF_max = intersect(b_PF_max, b_PF_max2);
wun_PF_max = union(b_PF_max, b_PF_max2);
jacc_PF_max = length(wint_PF_max) / length(wun_PF_max);

wint_PF_1 = intersect(b_PF_1, b_PF_12);
wun_PF_1 = union(b_PF_1, b_PF_12);
jacc_PF_1 = length(wint_PF_1) / length(wun_PF_1);

wint_PF_ML = intersect(b_PF_ML, b_PF_ML2);
wun_PF_ML = union(b_PF_ML, b_PF_ML2);
jacc_PF_ML = length(wint_PF_ML) / length(wun_PF_ML);


%%% FIND DEGREES AND STRENGTHS OF NODES IN BB BEFORE RESHUFFLING %%%%%%%%%%

W_DF = sparse(bb_DF(:,1), bb_DF(:,2), W(b_DF));
W_DF = sparse([full(W_DF) zeros(size(W_DF, 1), size(W,1)-size(W_DF,2))]);
W_DF = sparse([full(W_DF); zeros(size(W,2)-size(W_DF,1), size(W_DF,2))]);

k_in_DF  = full(sum(W_DF)); % Degree sequence (in-degree) (sum over column)
s_in_DF  = full(sum(W_DF)); % Stregth degree (in-stregth) 
k_out_DF = full(sum(W_DF'>0)); % Degree sequence (in-degree) (sum over column)
s_out_DF = full(sum(W_DF')); % Stregth degree (in-stregth) 
% fraction w / w average (w av. = s/k)
wei_DF = W(b_DF) ./ s_in_DF(bb_DF(:,2))' .*  k_in_DF(bb_DF(:,2))'; 

%%%%

W_HF = sparse(bb_HF(:,1), bb_HF(:,2), W(b_HF));
W_HF = sparse([full(W_HF) zeros(size(W_HF, 1), size(W,1)-size(W_HF,2))]);
W_HF = sparse([full(W_HF); zeros(size(W,2)-size(W_HF,1), size(W_HF,2))]);

k_in_HF  = full(sum(W_HF)); % Degree sequence (in-degree) (sum over column)
s_in_HF  = full(sum(W_HF)); % Stregth degree (in-stregth) 
k_out_HF = full(sum(W_HF'>0)); % Degree sequence (in-degree) (sum over column)
s_out_HF = full(sum(W_HF')); % Stregth degree (in-stregth) 
% fraction w / w average (w av. = s/k)
wei_HF = W(b_HF) ./ s_in_HF(bb_HF(:,2))' .*  k_in_HF(bb_HF(:,2))'; 

%%%%

W_PF_min = sparse(bb_PF_min(:,1), bb_PF_min(:,2), W(b_PF_min));
W_PF_min = sparse([full(W_PF_min) zeros(size(W_PF_min, 1), size(W,1)-size(W_PF_min,2))]);
W_PF_min = sparse([full(W_PF_min); zeros(size(W,2)-size(W_PF_min,1), size(W_PF_min,2))]);

k_in_PF_min  = full(sum(W_PF_min)); % Degree sequence (in-degree) (sum over column)
s_in_PF_min  = full(sum(W_PF_min)); % Stregth degree (in-stregth) 
k_out_PF_min = full(sum(W_PF_min'>0)); % Degree sequence (in-degree) (sum over column)
s_out_PF_min = full(sum(W_PF_min')); % Stregth degree (in-stregth)
% fraction w / w average (w av. = s/k)
wei_PF_min = W(b_PF_min) ./ s_in_PF_min(bb_PF_min(:,2))' .*  k_in_PF_min(bb_PF_min(:,2))'; 

%%%%

W_PF_max = sparse(bb_PF_max(:,1), bb_PF_max(:,2), W(b_PF_max));
W_PF_max = sparse([full(W_PF_max) zeros(size(W_PF_max, 1), size(W,1)-size(W_PF_max,2))]);
W_PF_max = sparse([full(W_PF_max); zeros(size(W,2)-size(W_PF_max,1), size(W_PF_max,2))]);

k_in_PF_max  = full(sum(W_PF_max)); % Degree sequence (in-degree) (sum over column)
s_in_PF_max  = full(sum(W_PF_max)); % Stregth degree (in-stregth) 
k_out_PF_max = full(sum(W_PF_max'>0)); % Degree sequence (in-degree) (sum over column)
s_out_PF_max = full(sum(W_PF_max')); % Stregth degree (in-stregth)
% fraction w / w average (w av. = s/k)
wei_PF_max = W(b_PF_max) ./ s_in_PF_max(bb_PF_max(:,2))' .*  k_in_PF_max(bb_PF_max(:,2))'; 

%%%%

W_PF_1 = sparse(bb_PF_1(:,1), bb_PF_1(:,2), W(b_PF_1));
W_PF_1 = sparse([full(W_PF_1) zeros(size(W_PF_1, 1), size(W,1)-size(W_PF_1,2))]);
W_PF_1 = sparse([full(W_PF_1); zeros(size(W,2)-size(W_PF_1,1), size(W_PF_1,2))]);

k_in_PF_1  = full(sum(W_PF_1)); % Degree sequence (in-degree) (sum over column)
s_in_PF_1  = full(sum(W_PF_1)); % Stregth degree (in-stregth) 
k_out_PF_1 = full(sum(W_PF_1'>0)); % Degree sequence (in-degree) (sum over column)
s_out_PF_1 = full(sum(W_PF_1')); % Stregth degree (in-stregth)
% fraction w / w average (w av. = s/k)
wei_PF_1 = W(b_PF_1) ./ s_in_PF_1(bb_PF_1(:,2))' .*  k_in_PF_1(bb_PF_1(:,2))'; 


%%%%

W_PF_ML = sparse(bb_PF_ML(:,1), bb_PF_ML(:,2), W(b_PF_ML));
W_PF_ML = sparse([full(W_PF_ML) zeros(size(W_PF_ML, 1), size(W,1)-size(W_PF_ML,2))]);
W_PF_ML = sparse([full(W_PF_ML); zeros(size(W,2)-size(W_PF_ML,1), size(W_PF_ML,2))]);

k_in_PF_ML  = full(sum(W_PF_ML)); % Degree sequence (in-degree) (sum over column)
s_in_PF_ML  = full(sum(W_PF_ML)); % Stregth degree (in-stregth) 
k_out_PF_ML = full(sum(W_PF_ML'>0)); % Degree sequence (in-degree) (sum over column)
s_out_PF_ML = full(sum(W_PF_ML')); % Stregth degree (in-stregth)
% fraction w / w average (w av. = s/k)
wei_PF_ML = W(b_PF_ML) ./ s_in_PF_ML(bb_PF_ML(:,2))' .*  k_in_PF_ML(bb_PF_ML(:,2))'; 



%%% FIND DEGREES AND STRENGTHS OF NODES IN BB AFTER RESHUFFLING %%%%%%%%%%%

W_DF2 = sparse(bb_DF2(:,1), bb_DF2(:,2), W_shuf(b_DF2));
W_DF2 = sparse([full(W_DF2) zeros(size(W_DF2, 1), size(W,1)-size(W_DF2,2))]);
W_DF2 = sparse([full(W_DF2); zeros(size(W,2)-size(W_DF2,1), size(W_DF2,2))]);

k_in_DF2  = full(sum(W_DF2)); % Degree sequence (in-degree) (sum over column)
s_in_DF2  = full(sum(W_DF2)); % Stregth degree (in-stregth) 
k_out_DF2 = full(sum(W_DF2'>0)); % Degree sequence (in-degree) (sum over column)
s_out_DF2 = full(sum(W_DF2')); % Stregth degree (in-stregth) 
% fraction w / w average (w av. = s/k)
wei_DF2 = W_shuf(b_DF2) ./ s_in_DF2(bb_DF2(:,2))' .*  k_in_DF2(bb_DF2(:,2))'; 

%%%%

W_HF2 = sparse(bb_HF2(:,1), bb_HF2(:,2), W_shuf(b_HF2));
W_HF2 = sparse([full(W_HF2) zeros(size(W_HF2, 1), size(W,1)-size(W_HF2,2))]);
W_HF2 = sparse([full(W_HF); zeros(size(W,2)-size(W_HF2,1), size(W_HF2,2))]);

k_in_HF2  = full(sum(W_HF2)); % Degree sequence (in-degree) (sum over column)
s_in_HF2  = full(sum(W_HF2)); % Stregth degree (in-stregth) 
k_out_HF2 = full(sum(W_HF2'>0)); % Degree sequence (in-degree) (sum over column)
s_out_HF2 = full(sum(W_HF2')); % Stregth degree (in-stregth) 
% fraction w / w average (w av. = s/k)
wei_HF2 = W_shuf(b_HF2) ./ s_in_HF2(bb_HF2(:,2))' .*  k_in_HF2(bb_HF2(:,2))'; 


%%%%

W_PF_min2 = sparse(bb_PF_min2(:,1), bb_PF_min2(:,2), W_shuf(b_PF_min2));
W_PF_min2 = sparse([full(W_PF_min2) zeros(size(W_PF_min2, 1), size(W,1)-size(W_PF_min2,2))]);
W_PF_min2 = sparse([full(W_PF_min2); zeros(size(W,2)-size(W_PF_min2,1), size(W_PF_min2,2))]);

k_in_PF_min2  = full(sum(W_PF_min2)); % Degree sequence (in-degree) (sum over column)
s_in_PF_min2  = full(sum(W_PF_min2)); % Stregth degree (in-stregth) 
k_out_PF_min2 = full(sum(W_PF_min2'>0)); % Degree sequence (in-degree) (sum over column)
s_out_PF_min2 = full(sum(W_PF_min2')); % Stregth degree (in-stregth)
% fraction w / w average (w av. = s/k)
wei_PF_min2 = W_shuf(b_PF_min2) ./ s_in_PF_min2(bb_PF_min2(:,2))' .*  k_in_PF_min2(bb_PF_min2(:,2))'; 

%%%%

W_PF_max2 = sparse(bb_PF_max2(:,1), bb_PF_max2(:,2), W_shuf(b_PF_max2));
W_PF_max2 = sparse([full(W_PF_max2) zeros(size(W_PF_max2, 1), size(W,1)-size(W_PF_max2,2))]);
W_PF_max2 = sparse([full(W_PF_max2); zeros(size(W,2)-size(W_PF_max2,1), size(W_PF_max2,2))]);

k_in_PF_max2  = full(sum(W_PF_max2)); % Degree sequence (in-degree) (sum over column)
s_in_PF_max2  = full(sum(W_PF_max2)); % Stregth degree (in-stregth) 
k_out_PF_max2 = full(sum(W_PF_max2'>0)); % Degree sequence (in-degree) (sum over column)
s_out_PF_max2 = full(sum(W_PF_max2')); % Stregth degree (in-stregth)
% fraction w / w average (w av. = s/k)
wei_PF_max2 = W_shuf(b_PF_max2) ./ s_in_PF_max2(bb_PF_max2(:,2))' .*  k_in_PF_max2(bb_PF_max2(:,2))'; 

%%%%

W_PF_12 = sparse(bb_PF_12(:,1), bb_PF_12(:,2), W_shuf(b_PF_12));
W_PF_12 = sparse([full(W_PF_12) zeros(size(W_PF_12, 1), size(W,1)-size(W_PF_12,2))]);
W_PF_12 = sparse([full(W_PF_12); zeros(size(W,2)-size(W_PF_12,1), size(W_PF_12,2))]);

k_in_PF_12  = full(sum(W_PF_12)); % Degree sequence (in-degree) (sum over column)
s_in_PF_12  = full(sum(W_PF_12)); % Stregth degree (in-stregth) 
k_out_PF_12 = full(sum(W_PF_12'>0)); % Degree sequence (in-degree) (sum over column)
s_out_PF_12 = full(sum(W_PF_12')); % Stregth degree (in-stregth)
% fraction w / w average (w av. = s/k)
wei_PF_12 = W_shuf(b_PF_12) ./ s_in_PF_12(bb_PF_12(:,2))' .*  k_in_PF_12(bb_PF_12(:,2))'; 

%%%%

W_PF_ML2 = sparse(bb_PF_ML2(:,1), bb_PF_ML2(:,2), W_shuf(b_PF_ML2));
W_PF_ML2 = sparse([full(W_PF_ML2) zeros(size(W_PF_ML2, 1), size(W,1)-size(W_PF_ML2,2))]);
W_PF_ML2 = sparse([full(W_PF_ML2); zeros(size(W,2)-size(W_PF_ML2,1), size(W_PF_ML2,2))]);

k_in_PF_ML2  = full(sum(W_PF_ML2)); % Degree sequence (in-degree) (sum over column)
s_in_PF_ML2  = full(sum(W_PF_ML2)); % Stregth degree (in-stregth) 
k_out_PF_ML2 = full(sum(W_PF_ML2'>0)); % Degree sequence (in-degree) (sum over column)
s_out_PF_ML2 = full(sum(W_PF_ML2')); % Stregth degree (in-stregth)
% fraction w / w average (w av. = s/k)
wei_PF_ML2 = W_shuf(b_PF_ML2) ./ s_in_PF_ML2(bb_PF_ML2(:,2))' .*  k_in_PF_ML2(bb_PF_ML2(:,2))'; 





%%% 6) PLOT THE FRACTION OF LINKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);

num = linspace(1, 6, 6);

semilogy(num, fr_link, 'diamond', 'MarkerSize', 20, 'color', '#6D69DF', 'LineWidth',2)
hold on
for i=1:Nsims
    semilogy(num, fr_link2(i,:), '.', 'MarkerSize', 20, 'LineWidth',2)
end

xticks(num)
xticklabels({'DF', 'HF', 'PF MIN', 'PF MAX', 'PF 1', 'PF ML'})
axis([0 7 1e-4 1]) % set axis limits - first 2 for the x second 2 for the y 

title('Reshuffling test on WIOT network')
xlabel('Filter')
ylabel('Fraction of links retained in the backbone')
legend('before', 'after')
set(gca, 'FontName', 'latex')

% savefig(gcf, 'PLOT/reshuffling_WIOT.fig')


%%% 6) PLOT THE HISTOGRAM OF DEGREE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
histogram(log(wei_DF), 30,'Normalization', 'pdf', 'FaceColor', 'blue', 'EdgeColor', 'blue', 'EdgeAlpha',0.7,'FaceAlpha',0.5 )
hold on
histogram(log(wei_DF2), 30, 'Normalization', 'pdf', 'FaceColor', 'red', 'EdgeColor', 'red','EdgeAlpha',0.7,'FaceAlpha',0.5)
title('\textbf{Disparity filter}', 'Interpreter','latex', 'FontWeight','bold')
xlabel('$\mathbf{\frac{w}{\langle w \rangle }}$', 'Interpreter', 'latex', 'FontSize', 15)
ylabel('\textbf{frequency}', 'Interpreter','latex')
legend('before', 'after', 'Interpreter', 'latex', 'Location','northwest')
set(gca, 'FontName', 'latex', 'TickLabelInterpreter', 'latex')

savefig(gcf, 'PLOT/w-wmed_DF.fig')

figure
histogram(log(wei_HF), 30,'Normalization', 'pdf', 'FaceColor', 'blue', 'EdgeColor', 'blue', 'EdgeAlpha',0.7,'FaceAlpha',0.5 )
hold on
histogram(log(wei_HF2), 30, 'Normalization', 'pdf', 'FaceColor', 'red', 'EdgeColor', 'red','EdgeAlpha',0.7,'FaceAlpha',0.5)
title('\textbf{Hypergeometric filter}', 'Interpreter','latex', 'FontWeight','bold')
xlabel('$\mathbf{\frac{w}{\langle w \rangle }}$', 'Interpreter', 'latex', 'FontSize', 15)
ylabel('\textbf{frequency}', 'Interpreter','latex')
legend('before', 'after', 'Interpreter', 'latex', 'Location','northwest')
set(gca, 'FontName', 'latex', 'TickLabelInterpreter', 'latex')

savefig(gcf, 'PLOT/w-wmed_HF.fig')

figure
histogram(log(wei_PF_min), 30,'Normalization', 'pdf', 'FaceColor', 'blue', 'EdgeColor', 'blue', 'EdgeAlpha',0.7,'FaceAlpha',0.5 )
hold on
histogram(log(wei_PF_min2), 30, 'Normalization', 'pdf', 'FaceColor', 'red', 'EdgeColor', 'red','EdgeAlpha',0.7,'FaceAlpha',0.5)
title('\textbf{Polya filter with $a_{min}$}', 'Interpreter','latex', 'FontWeight','bold')
xlabel('$\mathbf{\frac{w}{\langle w \rangle }}$', 'Interpreter', 'latex', 'FontSize', 15)
ylabel('\textbf{frequency}', 'Interpreter','latex')
legend('before', 'after', 'Interpreter', 'latex', 'Location','northwest')
set(gca, 'FontName', 'latex', 'TickLabelInterpreter', 'latex')

savefig(gcf, 'PLOT/w-wmed_PF_min.fig')


figure
histogram(log(wei_PF_max), 30,'Normalization', 'pdf', 'FaceColor', 'blue', 'EdgeColor', 'blue', 'EdgeAlpha',0.7,'FaceAlpha',0.5 )
hold on
histogram(log(wei_PF_max2), 30, 'Normalization', 'pdf', 'FaceColor', 'red', 'EdgeColor', 'red','EdgeAlpha',0.7,'FaceAlpha',0.5)
title('\textbf{Polya filter with $a_{max}$}', 'Interpreter','latex', 'FontWeight','bold')
xlabel('$\mathbf{\frac{w}{\langle w \rangle }}$', 'Interpreter', 'latex', 'FontSize', 15)
ylabel('\textbf{frequency}', 'Interpreter','latex')
legend('before', 'after', 'Interpreter', 'latex', 'Location','northwest')
set(gca, 'FontName', 'latex', 'TickLabelInterpreter', 'latex')

savefig(gcf, 'PLOT/w-wmed_PF_max.fig')


figure
histogram(log(wei_PF_1), 30,'Normalization', 'pdf', 'FaceColor', 'blue', 'EdgeColor', 'blue', 'EdgeAlpha',0.7,'FaceAlpha',0.5 )
hold on
histogram(log(wei_PF_12), 30, 'Normalization', 'pdf', 'FaceColor', 'red', 'EdgeColor', 'red','EdgeAlpha',0.7,'FaceAlpha',0.5)
title('\textbf{Polya filter with $a = 1$}', 'Interpreter','latex', 'FontWeight','bold')
xlabel('$\mathbf{\frac{w}{\langle w \rangle }}$', 'Interpreter', 'latex', 'FontSize', 15)
ylabel('\textbf{frequency}', 'Interpreter','latex')
legend('before', 'after', 'Interpreter', 'latex', 'Location','northwest')
set(gca, 'FontName', 'latex', 'TickLabelInterpreter', 'latex')

savefig(gcf, 'PLOT/w-wmed_PF_1.fig')

figure
histogram(log(wei_PF_ML), 30,'Normalization', 'pdf', 'FaceColor', 'blue', 'EdgeColor', 'blue', 'EdgeAlpha',0.7,'FaceAlpha',0.5 )
hold on
histogram(log(wei_PF_ML2), 30, 'Normalization', 'pdf', 'FaceColor', 'red', 'EdgeColor', 'red','EdgeAlpha',0.7,'FaceAlpha',0.5)
title('\textbf{Polya filter with $a_{ML}$}', 'Interpreter','latex', 'FontWeight','bold')
xlabel('$\mathbf{\frac{w}{\langle w \rangle }}$', 'Interpreter', 'latex', 'FontSize', 15)
ylabel('\textbf{frequency}', 'Interpreter','latex')
legend('before', 'after', 'Interpreter', 'latex', 'Location','northwest')
set(gca, 'FontName', 'latex', 'TickLabelInterpreter', 'latex')

savefig(gcf, 'PLOT/w-wmed_PF_ML.fig')


