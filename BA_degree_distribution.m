clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file generates a Barabasi-Albert network and then it plots its 
%%% degree distribution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 4; % Number of links formed by new nodes
N = 2000;
Nsims = 1;

deg = [];

for ns = 1:Nsims
    
    A = BA_network(N,m);

    deg = [deg; full(sum(A))'];

end

h = histogram(deg,100,'Normalization','pdf');

x = h.BinEdges;
x = x(2:end);
y = h.Values;

f = find(y == 0);
x(f) = []; y(f) = [];

loglog(x,y,'ob','MarkerSize',6,'MarkerFaceColor','b')
hold on

c = polyfit(log(x),log(y),1);

loglog(x, exp(c(2))*x.^c(1), 'r-','LineWidth',1.5)

xlabel('$k$','Interpreter','latex')
ylabel('$p(k)$','Interpreter','latex')
set(gca,'FontSize',22)

