function [x_hat] = ifft_plus(X)

    X_mod = zeros(1, length(X));
    
    X_mod(1:length(X)/2) = X(length(X)/2+1:end);
    X_mod(length(X)/2+1:end) = X(1:length(X)/2);
    
    x_hat_trans = ifft(X_mod);
    x_hat = length(x_hat_trans)*x_hat_trans;

end