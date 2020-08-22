function D = mdct(X)
% forward MDCT of the matrix X

% Determine block size and shape
N = size(X,1)/2;
k = (0:N-1)';
n = 0:2*N-1;

% Determine coefficients
F = cos((pi/N*(k + 0.5))*(n + 0.5 + N/2));
D = F*X;

end 

