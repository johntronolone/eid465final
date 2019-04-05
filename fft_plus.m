function [X, f] = fft_plus(x, fs, n)

    % regular fft and normalize
    if(n)
        X_old = fft(x, n)/fs;
    else
        X_old = fft(x)/fs; %/length(x);
    end
    
    % shift fft values greater than fs/2
    X = zeros(1, length(X_old));
    
    %left side
    X(1:ceil(length(X_old)/2)-1) = X_old(floor(length(X_old)/2)+2:n);
    
    %right side
    X(ceil(length(X_old)/2):end) = X_old(1:floor(length(X_old)/2)+1);
    
    % create freq vector
    df = fs/length(x);
    f = df*(-length(X)/2+1:length(X)/2);
    
end