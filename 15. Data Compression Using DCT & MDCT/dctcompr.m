function [y, ra] = dctcompr(x, N, r, method);
% Implements the DCT compression on input signal x using 1 of 2 methods.
% If method = 1, a fixed number of values are kept per column
% If method = 2, values above a certain threshold are kept.
% N = number of input samples per block
% r = compression ration (between 0 and 1)

X = buffer(x,N);
D = dct(X);
M = size(D,2);

if method == 1          % Keep a fixed number of values per column
    Nr = round(r*N);     % Number of values kept in each frame            
    ra = Nr/N;          % actual compression ratio
    [~,Ir] = sort(abs(D), 'descend'); % Sort columns descending
    Ir = Ir(1:Nr,:);                   % Keep highest Nr values
    C = zeros(size(D));                 % kept DCT coefficients
    for m = 1:M
        Dr(:,m) = D(Ir(:,m),m);         % Nr x M matrix of sorted coeficients
        C(Ir(:,m),m) = Dr(:,m);         % C = new DCT coefficients
    end
end
   
if method == 2          % Keep values above a certain threshold
    Dthr = r*max(max(abs(D)));          % ratio times the max value in D
    I = find(abs(D) < Dthr);            % Indices of coefficients to be discarded
    C = D;                              % Copy of D
    C(I) = 0;                           % Discard coefficients below Dthr
    ra = 1 - length(I)/(N*M);           % Actual compression ratio
end
    
Y = idct(C);              % Take inverse dct of kept coefficients
y = Y(:);                 % Concatenate columns
y = y(1:length(x));       % Make output the same length as the input

end
