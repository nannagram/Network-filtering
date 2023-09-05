clear all 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file is aimed to approximate the behaviour of the p-value for the
%%% hypergeometric filter in the regime N >> Na >> Nb. 
%%% The main assumption is that we can get rid of the sum considering it as
%%% an average among some values. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate an array of 200 values for N 
N = logspace(log10(200), 5, 200);

% keep fixed Na, Nb, w in the correct regime 
Na = 70; Nb = 30;
w = 10; 

% initialize the array for the real and the approximated p-values 
pr = [];
pa = [];

% loop over N
for i=1:length(N)
    
    % array of values in the sum 
    valx = linspace(w, Nb, Nb-w+1);
    % array of values of the function for each term in the sum 
    val = [];

    %  loop over the values in the sum 
    for j=1:length(valx)
        % compute the real value 
        val(j) = nchoosek(Na,valx(j))*nchoosek(round(N(i))-Na, Nb-valx(j) );
    end
   
    % compute the mean value
    mval = mean(val);

    % find the w that has the closest value to the mean
    [c index] = min(abs(val-mval));
    w_approx = valx(index);

    % compute the function for w_approx
    m = nchoosek(Na,w_approx)*nchoosek(round(N(i))-Na, Nb-w_approx );
    
%     semilogy(valx, val)
%     hold on 

    % compute the real p-value and append it 
    p_real = hygecdf(w-1, round(N(i)), Na, Nb, 'upper');
    pr = [pr, p_real];

    % compute the approximated p-value and append it 
    p_approx = (nchoosek(round(N(i)),Nb))^-1 * (Nb-w) * m;
    pa = [pa, p_approx];
end

% plot the p-values in function of N
loglog(N,pr, 'color',[162/255 158/255 254/255], 'LineWidth',1.5)
hold on
loglog(N,pa, '--','color', [105/255 95/255 255/255], 'LineWidth',1.5)
hold on

% plot the power-law fit 
loglog(N, 10^30 * N.^(-w_approx), '-','color', [0 0 0], 'LineWidth',0.5)

legend('real', 'approximated', 'fit', 'interpreter', 'latex')
xlabel('S', 'Interpreter', 'latex')
ylabel('p-value', 'Interpreter','latex')




