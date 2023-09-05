clear all 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file aims to plot the graph representation of the US air traffic
%%% network over the USA map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/PLOT/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/PLOTPLOT/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/MAX&SAM/MAX&SAM-routines/")


%load flight network with coordinates
T = readtable("new.csv");
% col 1:  passengers
% col 2:  distance
% col 5:  origin airport ID (NUM)
% col 6:  origin airport ID (LETT)
% col 7:  origin airport city name
% col 8:  destination airport ID (NUM)
% col 9:  destination airport ID (LETT)
% col 10: destination airport city name
% col 11: year (2017)
% col 12: origin airport latitude
% col 13: origin airport longitude
% col 14: destination airport latitude
% col 15: destination airport longitude

%%% find and remove links without passengers
aux = find(T{:,1}==0);
T(aux, :) = [];

%%% IN MANIERA MOLTO MESCHINA, TOGLI I VOLI DA E PER AEROPORTI DI CUI NON
%%% CONOSCI LE COORDINATE :)
aux = find(isnan(T{:,'DEPARTURE_LAT'}));
T(aux, :) = [];

aux = find(isnan(T{:,'DESTINATION_LAT'}));
T(aux, :) = [];


%%% ANALOGUE TO THE SPARSE COMMAND
% N = max(T{:, 5});
% W = sparse(N,N);
% for i=1:size(T,1)
%     dep = T{i,5};
%     arr = T{i,8};
%     W(dep, arr) = W(dep, arr) + T{i,1};
% end


%%% CREATE THE ADIACENCY MATRIX: 
%%% NOTE: THE DIFFERENT FLIGHT FROM THE SAME AIRPORT HAVE THEIR WEIGHTS
%%% SUMMED 
W = sparse(T{:,5}, T{:,8}, T{:,1});

alpha = 0.05;    % univariate significance level
a1 = 0.1;       % Polya parameter min
a2 = 4.5;          % Polya parameter max
apr_lvl = 1e40;  % approximation level polya

L = nnz(W);      % Computing the number of links in A
alpha = alpha/L; % Bonferroni correction
alpha = 1e-4;



%%% FILTER THE NETWORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [bb, p] = disp_filter(W, alpha);
% [bb, p] = hypergeom_filter(W, alpha);
% [bb, p] = PF(W, a1, alpha,apr_lvl, 0);
[bb, p] = PF(W, a2, alpha,apr_lvl, 0);
% [bb, p] = PF(W, 1, alpha,apr_lvl, 0);
% [bb, p] = PF(W, -1, alpha, apr_lvl, 0);

disp('end filtering')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mah1 = ismember(T{:,5}, bb(:,1));
mah2 = ismember(T{:,8}, bb(:,2));
mah = mah1 + mah2;
ind = find(mah==2);
T2 = T(ind, :);

% % convert indeces of links to a one index 
link_ind = sub2ind(size(W), bb(:,1), bb(:,2));
% 
% [boh, ind] = intersect(T{:,5}, bb(:,1));
% T2 = T(ind, :);
% [boh, ind] = intersect(T2{:,8}, bb(:,2));
% T2 = T2(ind, :);


dep_node = num2cell( bb(:,1) );
arr_node = num2cell( bb(:,2) );
weights  = num2cell( W(link_ind) );
dep_id   = num2cell( zeros(length(link_ind), 1) );
arr_id   = num2cell( zeros(length(link_ind), 1) );

for i = 1:length(link_ind)
    de = bb(i,1);
    ar = bb(i,2);

    ind1 = find(T{:,5} == de);
    dep_id{i} = T{ind1(1), 6};

    ind2 = find(T{:,8} == ar);
    arr_id{i} = T{ind2(1), 9};
end

T3 = [weights dep_node dep_id arr_node arr_id];
T3 = cell2table(T3);



% load coordinates of airport
coord1 = readtable("USA_Airports_big.csv");
coord2 = readtable("USA_Airports_medium.csv");
coord3 = readtable("USA_Airports_small.csv");
% col 1: longitude
% col 2: latitude
% col 3: airport ID (number)
% col 4: airport ID (letters)
% col 5: airport name 
% col 11: altitude


%%% WRITE NODES LIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
node_list = [T2{1, "ORIGIN"}];
lat = [T2{1, "DEPARTURE_LAT"}];
long = [T2{1, 'DEPARTURE_LONG'}];

for i = 1:size(T2,1)
    % find ID of airports of origin and departure ( 3 letters )
    orig = T2{i, "ORIGIN"};

    if ismember(orig, node_list) == 0
        node_list = [node_list; orig];

        lat  = [lat; T2{i, "DEPARTURE_LAT"} ];
        long = [long; T2{i, 'DEPARTURE_LONG'} ];
    end
end
for i = 1:size(T2,1)
    % find ID of airports of origin and departure ( 3 letters )
    dest = T2{i, "DEST"};
    if ismember(dest, node_list) == 0
        node_list = [node_list; dest];

        lat  = [lat;  T2{i, "DESTINATION_LAT"} ];
        long = [long; T2{i, 'DESTINATION_LONG'} ];
    end
end


disp('end writing node list')


%%% define the limits of the image

[latlim,lonlim] = geoquadline(lat,long);
[latlim,lonlim] = bufgeoquad(latlim,lonlim,5,5);

% latlim = [8.4840   76.2849];
% lonlim = [175, 0];

%%% read an image of the map
[A,RA] = readBasemapImage("bluegreen",latlim,lonlim,25);
disp('downloaded basemap')

%%% project in mercator coordinates the airports 
[x,y] = projfwd(RA.ProjectedCRS,lat,long);

% %%% build the figure 
% figure
% mapshow(A,RA)
% hold on
% axis off
% 
% %%% plot the digraph
% G = digraph(T2{:,6}, T2{:,9}, T2{:,1});
% 
% lineWidth = G.Edges.Weight/max(G.Edges.Weight)*5;
% 
% plot(G,XData=x,YData=y,LineWidth=lineWidth, NodeFontSize=11, MarkerSize=1.5, ...
%     EdgeAlpha=0.4,ArrowSize=3, NodeColor="#FA4D62",EdgeColor="#625DFB")


figure
mapshow(A, RA)
hold on 
axis off

G2 = digraph(T3{:,3}, T3{:,5}, full(T3{:,1}));
G2 = reordernodes(G2,node_list);

lineWidth2 = G2.Edges.Weight/max(G2.Edges.Weight)*5;

plot(G2,XData=x,YData=y,LineWidth=lineWidth2, NodeFontSize=7, MarkerSize=1.5, ...
    EdgeAlpha=0.7,ArrowSize=3, NodeColor="#FA4D62",EdgeColor="#625DFB")



% G2 = digraph(bb, W(link_ind));
% 
% lineWidth = G.Edges.Weight/max(G.Edges.Weight)*5;
% 
% plot(G2,XData=x,YData=y,LineWidth=lineWidth, NodeFontSize=11, MarkerSize=1.5, ...
%     EdgeAlpha=0.4,ArrowSize=3, NodeColor="#FA4D62",EdgeColor="#625DFB")

% gx = geoaxes('Basemap','topographic');
% 
% hold on 
% % geoscatter requires first the latitude than longitude 
% geoscatter(coord1{:,"Y"}, coord1{:,"X"}, 10, 'green', 'filled','o')
% geoscatter(coord2{:,"Y"}, coord2{:,"X"}, 10, 'blue','filled','o')
% geoscatter(coord3{:,"Y"}, coord3{:,"X"}, 10, 'red','filled','o')
