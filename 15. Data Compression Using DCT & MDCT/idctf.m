function X = idctf(D)
% Computes the IDCT of D

[N,M] = size(D);    % column and row lengths

k = (0:2*N-1)';    % Vectorize Ck_ext index k

% Construct normalizing, s vector
s0 = sqrt(N);   
sk = sqrt(N/2)*ones(N-1,1);
s = [s0;sk];

S = repmat(s,1,M);  % Scale D to Ck
Ck = S.*D;

Ck_ext = [Ck; zeros(1,M); -flipud(Ck(2:N,:))]; % Extend Ck 

% Reconstruct Yk from Ck_ext
weights = 2*exp(j*pi.*k/(2*N));                 
WEIGHTS = repmat(weights,1,M);

Yk = WEIGHTS.*Ck_ext; 

Y = real(ifft(Yk));     % Take the inverse fft

X = Y(1:N,:);           % Keep the first N values of each column

end
