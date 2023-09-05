%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function is aimed to evaluate the lagrangian multiplier for a direct
%%% weighted network in the configuration model (non enhanced) and the save
%%% them to a file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol] = Lag_mult(A)

    A = round(A);

    % number of nodes
    N = length(A(:,1));
    % number of links
    L = nnz(A);

    % compute the x and y array of Lagrange Multipliers for output and input
    % strengths
    sol = MAXandSAM('DWCM', A, [], [], 1e-2, 0);

    %reshape in the form of a N x 2 matrix
    sol = reshape(sol, N, 2);

    savefile='/Users/ariannaarmanetti/Desktop/TESI/Codici/FILTERS/filter-algorithms/MAX&SAM/MAX&SAM-routines/solxy.mat';
    save(savefile,'sol')

end 
