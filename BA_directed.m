%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS CODE GENERATE A DIRECTED BARABASI-ALBERT NETWORK
%%% It uses the copying mechanisms to generate new links directed toward
%%% already existing nodes, and for each link it choose randomly if
%%% directing it from the new node or towards it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% THINGS YET TO IMPLEMENT:
    % A fitness model
    % The possibility to generate links between two already existing nodes
    % A boost for new links
    % Some sort of ageing 
    % ...

clear all
close all

N = 3000;   % Number of final nodes
m0 = 6;     % Starting fully connected nodes
m = 3;      % Number of new links added at each step
p = 1;   % probability to connect to the randomly chosen link (1-p is the prob to copy one of the links of the chosen node)
alfa = 0.5; % probability of generating a link outgoing from the new node


%%% UNDIRECTED NETWORK
%A = BA_network(N,m);

%%% DIRECTED NETWORK
A = BA_Dnet(N,m0, m, p, alfa);

save('saveBA.mat', 'A');


%%% PLOT THE DEGREE DISTRIBUTION
k_in = full(sum(A));
k_out = full(sum(A'));

% incoming degree distribution
h_in = histogram(k_in,100,'Normalization','pdf');
x_in = h_in.BinEdges;
x_in = x_in(2:end);
y_in = h_in.Values;

f = find(y_in == 0);
x_in(f) = []; y_in(f) = [];

% outgoing degree distribution
h_out = histogram(k_out, 100, 'Normalization','pdf');
x_out = h_out.BinEdges;
x_out = x_out(2:end);
y_out = h_out.Values;

f = find(y_out == 0);
x_out(f) = []; y_out(f) = [];

subplot(2,1,1)
loglog(x_in,y_in,'ob','MarkerSize',6,'MarkerFaceColor','b')
hold on
% power-law fit
% the coefficients are slightly different
c_in = polyfit(log(x_in),log(y_in),1);
loglog(x_in, 50*exp(c_in(2))*x_in.^-3, 'b-','LineWidth',1.5)

subplot(2,1,2)
loglog(x_out,y_out,'ob','MarkerSize',6,'MarkerFaceColor','r')
hold on 

%power-law fit 
% slightly different coeffients 
c_out = polyfit(log(x_out),log(y_out),1);
loglog(x_out, 50*exp(c_out(2))*x_out.^-3, 'r-','LineWidth',1.5)

xlabel('$k$','Interpreter','latex')
ylabel('$p(k)$','Interpreter','latex')
set(gca,'FontSize',22)


G = digraph(A);
plot(G, 'Layout','force')

