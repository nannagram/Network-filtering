clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file plots the degree distribution of an ER network
%%% In the following we check that for large N the degree distribution of
%%% the ER model is approximately Gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = 1000; % Number of nodes
p = 0.2; % Connectivity parameter
Nsims = 1; % Number of adjacency matrices of the ER network to generate

deg = []; % Array to store all degree sequences

for ns = 1:Nsims

    A = ER_network(N,p); % Generating one instance of an ER network's adjacency matrix

    k = full(sum(A)); % Computing the degree sequence : sum() sums over the columns, while full create an array of dimension N filling with 0 the elements missing 

    deg = [deg; k']; % Storing the degree sequence (it concatenates the simulations in adiacent columns)

end


%%% Plotting histogram vs Gaussian distribution

histogram(k,50,'Normalization','pdf')

hold on

m = N*p;
s = sqrt(N*p);
x = linspace(m-5*s,m+5*s,1000);
y = exp(-(x-m).^2/(2*s^2))/sqrt(2*pi*s^2);

bin = [];
x = min(k):max(k);
for i= 1:length(x)
    bin(i) = nchoosek((N-1), round(x(i))) * p^(x(i)) * (1-p)^(N-1-x(i));
end


plot(x,y,'r','LineWidth',2)

xlabel('$k$','Interpreter','latex')
ylabel('$p(k)$','Interpreter','latex')
set(gca,'FontSize',20)
