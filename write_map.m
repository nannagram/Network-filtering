%%% THIS FILE CAN BE USED TO LOAD COORDINATES OF THE AIRPORTS IN THE USA,
%%% THEN LOAD THE AERIAL TRAFFIC NETWORK AND ADD TO EACH AIRPORT IN THIS
%%% TABLE THE RELATIVE COORDINATES. 
%%% IT WRITES A NEW .CSV FILE
%%% THERE IS ALSO THE GEOPLOT OF AIRPORTS OVER THE WORLD MAP

clear all
close all

addpath("/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/network/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/PLOT/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/PLOTPLOT/", ...
        "/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/MAX&SAM/MAX&SAM-routines/")


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

% load flights network
T = readtable("T_T100D_MARKET_ALL_CARRIER.csv");

% initialize coloumns of airports' coordinates
% In the form of a cell array,filled with zeros
dep_lat  = num2cell( zeros(size(T,1),1) );
dep_long = num2cell( zeros(size(T,1),1) );
arr_lat  = num2cell( zeros(size(T,1),1) );
arr_long = num2cell( zeros(size(T,1),1) );

% loop over all the flights in the table
for i = 1:size(T,1)
    % find ID of airports of origin and departure ( 3 letters )
    orig = T{i, "ORIGIN"};
    dest = T{i,"DEST"};

    % find indeces of the same DEPARTURE airport in the list with the
    % coordinates, first check the bigger ones
    ind1 = find(strcmp(coord1{:,"FAA_ID"}, orig));

    if ind1 > 0
        % add latitude and longitude of dep airport 
        dep_lat{i} =  coord1{ind1,"Y"};
        dep_long{i} = coord1{ind1,"X"};

    else
        % check if it's in the middle ones
        ind2 = find(strcmp(coord2{:,"FAA_ID"}, orig));

        if ind2 > 0
            dep_lat{i} =  coord2{ind2,"Y"};
            dep_long{i} = coord2{ind2,"X"};

        else
            % lastly check if it's in the smaller ones
            ind3 = find(strcmp(coord3{:,"FAA_ID"}, orig));
            dep_lat{i} =  coord3{ind3,"Y"};
            dep_long{i} = coord3{ind3,"X"};
        end
    end
    
    % find indeces of the DESTINATION airport in the list with the
    % coordinates, first check the bigger ones (than, same as above)
    dest1 = find(strcmp(coord1{:,"FAA_ID"}, dest));

    if dest1 > 0
        arr_lat{i} =  coord1{dest1,"Y"};
        arr_long{i} = coord1{dest1,"X"};

    else
        dest2 = find(strcmp(coord2{:,"FAA_ID"}, dest));

        if dest2 > 0
            arr_lat{i} =  coord2{dest2,"Y"};
            arr_long{i} = coord2{dest2,"X"};

        else
            dest3 = find(strcmp(coord3{:,"FAA_ID"}, dest));
            arr_lat{i} =  coord3{dest3,"Y"};
            arr_long{i} = coord3{dest3,"X"};
        end
    end
end

% Add 4 new columns to the end of the table
numOfColumn = size(T, 2);

T.(numOfColumn+1) = dep_lat;
T.Properties.VariableNames{numOfColumn+1} = 'DEPARTURE_LAT'; % Change column name 

T.(numOfColumn+2) = dep_long;
T.Properties.VariableNames{numOfColumn+2} = 'DEPARTURE_LONG'; % Change column name 

T.(numOfColumn+3) = arr_lat;
T.Properties.VariableNames{numOfColumn+3} = 'DESTINATION_LAT'; % Change column name 

T.(numOfColumn+4) = arr_long;
T.Properties.VariableNames{numOfColumn+4} = 'DESTINATION_LONG'; % Change column name 

% Write to CSV file
writetable(T, 'network/new.csv')

gx = geoaxes('Basemap','topographic');

hold on 
% geoscatter requires first the latitude than longitude 
geoscatter(coord1{:,"Y"}, coord1{:,"X"}, 10, 'green', 'filled','o')
geoscatter(coord2{:,"Y"}, coord2{:,"X"}, 10, 'blue','filled','o')
geoscatter(coord3{:,"Y"}, coord3{:,"X"}, 10, 'red','filled','o')