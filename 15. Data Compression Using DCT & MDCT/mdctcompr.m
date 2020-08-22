function [y, ra] = mdctcompr(x, N, rth, win, beta)
% Applied MDCT compression to the input signal x
% N =   DCT block length
% rth = threshold value
% win and beta are inputs to the pbwin function, which 
% generates Princen-Bradley windows
% win = 0 generates a rectangular window
% win = 1 generates a sine window
% win = 2 generates a vorbis window
% win = 3 generates a KBD window
% beta = Kaiser shape parameter (for KBD window only)
% 
X = buffer(x, 2*N, N, 'nodelay');
M = size(X,2);
w = (pbwin(N, win, beta))';      
W = repmat(w,1,M);
D = mdct(W.*X);

% DCT Method-2
Dthr = rth*max(max(abs(D)));          % ratio times the max value in D
I = find(abs(D) < Dthr);            % Indices of coefficients to be discarded
C = D;                              % Copy of D
C(I) = 0;                           % Discard coefficients below Dthr
ra = 1 - length(I)/(N*M);           % Actual compression ratio

Y = W.*imdct(C);

% Overlap-add function adjusted for blocks of length 2N
N = size(Y,1)/2;    % Y has columns of length 2*N
M = size(Y,2);      % M columns
L = N*M + N;                    % Extended length of y
y = zeros(L,1);                  % Initialize output vector

n = (1:2*N)';                   % Vectorize columns

for m = 0:M-1                     % Number of blocks (columns)         
    y(m*N + n) = y(m*N + n) + Y(:,m+1);  % add vectors with overlap of N = 50%
end
% End of adjusted overlap-add function

y = y(1:length(x));     % Make output length the same as the input length

end