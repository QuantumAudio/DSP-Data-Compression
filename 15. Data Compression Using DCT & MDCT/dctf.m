function D = dctf(X)
% forward DCT of the matrix X

[N,M] = size(X);

Y = [X; flipud(X)]; % Construct y column vectors 

Yk = fft(Y);    % Take fft of each column

% Define weighting coefficients of DFT
k = (0:2*N-1)';                     % Vectorize DFT column index k
weights = (1/2)*exp(-j*pi*k/(2*N)); % Repeat for M columns 
WEIGHTS = repmat(weights,1,M); 
% Calculate C_ext as DFT of Yk
Ck_ext = WEIGHTS.*Yk;                        
% Keep first N elements of each column
Ck = Ck_ext(1:N,:);    

% Construct normalizing s vector
s0 = sqrt(N);   
sk = sqrt(N/2)*ones(N-1,1);
s = [s0;sk];
% Repeat for M columns
S = repmat(s,1,M);
            
% Normalize output matrix of DCT coefficients
D = real(Ck./S);

end