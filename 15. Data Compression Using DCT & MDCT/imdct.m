function X = imdct(D)
% Inverse MDCT of the matrix D

% Determine block size and shape
N = size(D,1);
k = (0:N-1)';
n = 0:2*N-1;

% Compute coefficients
F = cos((pi/N*(k + 0.5))*(n + 0.5 + N/2));
X = (1/N)*F'*D;

end