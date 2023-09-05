function [b,pvalues] = hypergeom_filter(A,alpha)
% this function returns the list of p-values computed and the edge list 
% of the backbone validated by the hypergeometric filter for a directed network


    A = round(A);

    k_in  = full(sum(A > 0));    % IN-Degree sequence
    k_out = full(sum(A' > 0));   % OUT-Degree sequence
    s_in  = full(sum(A));        % IN-Strength sequence
    s_out = full(sum(A'));       % OUT-Strength sequence

    N_in = sum(s_in);            % total incoming strength 
    N_out = sum(s_out);          % total outgoing strength 
    N = round(N_in + N_out);     % total strength

    [ind1,ind2] = find(A > 0); % Finding indices of non-zero entries in A (i.e., links)

    b = []; % Empty array to store links in backbone
    pvalues = []; % Empty array to store p-values

    for i = 1:length(ind1) % Loop on links
        if mod(i, floor(length(ind1)/length(A)*100)) < 1e-6
            disp(i/length(ind1)*100)
        end
        
        w = A(ind1(i),ind2(i)); % Weight on current link
    
        %p = hypergeom_cdf(w, N, s(ind1(i)), s(ind2(i))) % compute p-value
        p = hygecdf(w-1,N,round(s_out(ind1(i))),round(s_in(ind2(i))),'upper');
        pvalues = [pvalues; p];

        % If the p-value falls below the significance level in input, the
        % corresponding link is stored in the backbone
        if p < alpha
           b = [b; ind1(i) ind2(i)];  
        end
    end 

end

%%% MY IMPLEMENTATION OF THE HYPERGEOMETRIC CUMULATIVE %%%%%%%%%%%%%%%%%%%%

function p = hypergeom_cdf(x, M, K, N)
% this function compute the p-value for a hypergeometric distribution,
% applying the stirling approximation to M 

L = 5; %treshold after which apply the Stirling correction to the factorials 

% inizializzo c per i=0
c = exp( K*(log(K)-1) + (M-K)*(log(M-K)-1) + N*(log(N)-1) + ...
         (M-K)*(log(M-K)-1) - K*(log(K)-1) - N*(log(N)-1) - ...
         (M-K-N)*(log(M-K-N)-1) - M*(log(M)-1));

if K > L            % check su k
    if N > L        % check su N 
        if x > L    % check su x
            % itero da i = 1 a L-1
            for i = 1:L-1
                c = c + exp(K*(log(K)-1) + (M-K)*(log(M-K)-1) + N*(log(N)-1) + ...
                (M-K)*(log(M-K)-1) - log(factorial(i)) - log(factorial(K-i)) - log(factorial(N-i)) - ...
                (M-K-N+i)*(log(M-K-N+i)-1) - M*(log(M)-1));
            end
            % itero da i = L a x
            for i = 1:L-1
                c = c + exp(K*(log(K)-1) + (M-K)*(log(M-K)-1) + N*(log(N)-1) + ...
                (M-K)*(log(M-K)-1) - i*(log(i)-1) - log(factorial(K-i)) - log(factorial(N-i)) - ...
                (M-K-N+i)*(log(M-K-N+i)-1) - M*(log(M)-1));
            end
        else % caso x < L, N,K>L
            for i = 1:x-1
                c = c + exp(K*(log(K)-1) + (M-K)*(log(M-K)-1) + N*(log(N)-1) + ...
                (M-K)*(log(M-K)-1) - log(factorial(i)) - log(factorial(K-i)) - log(factorial(N-i)) - ...
                (M-K-N+i)*(log(M-K-N+i)-1) - M*(log(M)-1));
            end
        end

    else %CASO N < L
        if x > L    % check su x
            % itero da i = 1 a L-1
            for i = 1:L-1
                c = c + exp(K*(log(K)-1) + (M-K)*(log(M-K)-1) + log(factorial(N)) + ...
                (M-K)*(log(M-K)-1) - log(factorial(i)) - log(factorial(K-i)) - log(factorial(N-i)) - ...
                (M-K-N+i)*(log(M-K-N+i)-1) - M*(log(M)-1));
            end
            % itero da i = L a x
            for i = 1:L-1
                c = c + exp(K*(log(K)-1) + (M-K)*(log(M-K)-1) + log(factorial(N)) + ...
                (M-K)*(log(M-K)-1) - i*(log(i)-1) - log(factorial(K-i)) - log(factorial(N-i)) - ...
                (M-K-N+i)*(log(M-K-N+i)-1) - M*(log(M)-1));
            end
        else % caso x < L, N,K>L
            for i = 1:x-1
                c = c + exp(K*(log(K)-1) + (M-K)*(log(M-K)-1) + log(factorial(N)) + ...
                (M-K)*(log(M-K)-1) - log(factorial(i)) - log(factorial(K-i)) - log(factorial(N-i)) - ...
                (M-K-N+i)*(log(M-K-N+i)-1) - M*(log(M)-1));
            end
        end
    end

else % CASO K < L
    if N > L        % check su N 
        if x > L    % check su x
            % itero da i = 1 a L-1
            for i = 1:L-1
                c = c + exp(log(factorial(K)) + (M-K)*(log(M-K)-1) + N*(log(N)-1) + ...
                (M-K)*(log(M-K)-1) - log(factorial(i)) - log(factorial(K-i)) - log(factorial(N-i)) - ...
                (M-K-N+i)*(log(M-K-N+i)-1) - M*(log(M)-1));
            end
            % itero da i = L a x
            for i = 1:L-1
                c = c + exp(log(factorial(K)) + (M-K)*(log(M-K)-1) + N*(log(N)-1) + ...
                (M-K)*(log(M-K)-1) - i*(log(i)-1) - log(factorial(K-i)) - log(factorial(N-i)) - ...
                (M-K-N+i)*(log(M-K-N+i)-1) - M*(log(M)-1));
            end
        else % caso x < L, K<L, N>L
            for i = 1:x-1
                c = c + exp(log(factorial(K)) + (M-K)*(log(M-K)-1) + N*(log(N)-1) + ...
                (M-K)*(log(M-K)-1) - log(factorial(i)) - log(factorial(K-i)) - log(factorial(N-i)) - ...
                (M-K-N+i)*(log(M-K-N+i)-1) - M*(log(M)-1));
            end
        end

    else %CASO N < L
        if x > L    % check su x
            % itero da i = 1 a L-1
            for i = 1:L-1
                c = c + exp(log(factorial(K)) + (M-K)*(log(M-K)-1) + log(factorial(N)) + ...
                (M-K)*(log(M-K)-1) - log(factorial(i)) - log(factorial(K-i)) - log(factorial(N-i)) - ...
                (M-K-N+i)*(log(M-K-N+i)-1) - M*(log(M)-1));
            end
            % itero da i = L a x
            for i = 1:L-1
                c = c + exp(log(factorial(K))+ (M-K)*(log(M-K)-1) + log(factorial(N)) + ...
                (M-K)*(log(M-K)-1) - i*(log(i)-1) - log(factorial(K-i)) - log(factorial(N-i)) - ...
                (M-K-N+i)*(log(M-K-N+i)-1) - M*(log(M)-1));
            end
        else % caso x < L, N,K < L 
            for i = 1:x-1
                c = c + exp(log(factorial(K)) + (M-K)*(log(M-K)-1) + log(factorial(N)) + ...
                (M-K)*(log(M-K)-1) - log(factorial(i)) - log(factorial(K-i)) - log(factorial(N-i)) - ...
                (M-K-N+i)*(log(M-K-N+i)-1) - M*(log(M)-1));
            end
        end
    end
end


%for i = 0:x-1
%    c = c + exp( log(factorial(K)) + (M-K)*(log(M-K)-1) + log(factorial(N)) + ...
%        (M-K)*(log(M-K)-1) - log(factorial(i)) - log(factorial(K-i)) - log(factorial(N-i)) - ...
%        (M-K-N+i)*(log(M-K-N+i)-1) - M*(log(M)-1));
%end

p = 1 - c;

end
