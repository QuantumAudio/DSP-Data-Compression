% Normalized DCT Matrix A for arbitray value of N
A = @(n) sqrt(2/N)*[ones(1,N)/sqrt(2); cos(pi*(1:N-1)'*(((0:N-1) + 1/2)/N))];
% Define a random matrix X
X = rand(1024,100);

%% Built-in dct Implementation
tic
D = dct(X);
Xdct = idct(D);
t_dct = toc;

%% dctf Implementation
tic 
D = dctf(X);    % forward DCT
Xdctf = idctf(D);
t_dctf = toc;

%% Matrix Implementation
tic
N = 1024;
n = 1:N;
A_dct = A(n);
D = A_dct*X;
Xmat = A_dct'*D;
t_matrix = toc;

%% Comparing the Different Implementations
format long
disp('X = ');
disp(X(1:10,1:3));
disp('Xdct = ');
disp(Xdct(1:10,1:3));
disp('Xdctf = ');
disp(Xdctf(1:10,1:3));
disp('Xmat = ');
disp(Xmat(1:10,1:3));

XdctErr = abs(X - Xdct);
disp('XdctError = ');
disp(XdctErr(1:10,1:3));

XdctfErr = abs(X - Xdctf);
disp('XdctfError = ');
disp(XdctfErr(1:10,1:3));

XmatErr = abs(X - Xmat);
disp('XmatError = ');
disp(XmatErr(1:10,1:3));

% Example Input Signal
L = 40; t = (0:L-1)/L;
x = sin(10*t.^2) + 2*t;
N = 10;

%% Method-1, r = 0.4, r_actual = 0.4
method = 1;
r = 0.4;
[y, ra] = dctcompr(x,N,r,1);
close all
plot(t,x,'r-',t,y,'b.'); title('method-1, N = 10, r = 0.4')
xlabel('t'),legend('original','recovered')

%% Method-2, r = 0.014, r_actual = 0.4
method = 2;
r = 0.014;
[y, ra] = dctcompr(x,N,r,2);
close all
plot(t,x,'r-',t,y,'b.'); title('method-2, N = 10, r_{thr} = 0.014')
xlabel('t'),legend('original','recovered')

%% dctcompr, method-1 on audio signal
[x, fs] = audioread('flute2.wav'); % 4 sec audio sample
x = x(:,1); % as a column
method = 1;
r = 0.2;
[y, ra] = dctcompr(x,N,r,1);
y = y/max(abs(y));
% Listen to results
soundsc(y,fs);
% Output results
audiowrite('dctcompression.wav', [x', zeros(1,fs),y'],fs);

%% 2D-DCT on an image file
X = imread('cameraman.tif');    % read image, 256x256 matrix
D = dct2(X);                    % compute its 2D-DCT

Dmax = max(max(abs(D)));        % Dmax = 30393.4687
Dth = 10;                       % select a threshold
rth = Dth/Dmax;                 % with threshold factor, rth = 3.2902e-04
C = D;
C(abs(D)<Dth) = 0;              % compressed DCT
ra = length(find(C))/prod(size(C)); % actual compression ratio,  ra = 26617/65536 = 0.4061
Y = idct2(C);                       % inverse 2D-DCT
figure; imshowpair(X,Y,'montage')   % display images side by side

% DCT Basis Functions
% number of basis functions
% N = 8
N = 16;
% for plotting
k = 0:N-1;
n = k';
% basis function parameters
s0 = sqrt(N);   
sk = sqrt(N/2)*ones(N-1,1);
s = [s0;sk];
% generate basis functions
Fk = (1./s).*cos((pi*k/N).*(n + 0.5)); % Each column of Fk is a basis function of n

% Plot basis functions
for l = 1:N
subplot(6,4,l)
stairs(Fk(:,l));
xlim([0,N]);
end

%% Generating Princen-Bradley Windows
N = 4; beta = 5;
rectangular = pbwin(N,0,beta)
sine = pbwin(N,1,beta)
vorbis = pbwin(N,2,beta)
KBD = pbwin(N,3,beta)
%% Verifying Princen-Bradley Condition
rectangular(1:N).^2 + rectangular(N+1:2*N).^2
sine(1:N).^2 + sine(N+1:2*N).^2
vorbis(1:N).^2 + vorbis(N+1:2*N).^2
KBD(1:N).^2 + KBD(N+1:2*N).^2
%% Plot Results
N = 128, beta = 20;
n = 0:2*N-1;
close all
plot(n,sine,'r',n,vorbis,'g',n,KBD,'b'); axis([0,256,0,1.5]);
legend('sine','vorbis','KBD'); grid on; xlabel('n')
title('Princen-Bradley windows, 2N = 256, \beta = 20');

% Adusting rth until ra is in the range [0.15, 0.20] 
[x, fs] = audioread('flute2.wav'); % 4 sec audio sample
win = 1;
beta = 5;
rth = 0.0003;
N = 1024;
[y, ra] = mdctcompr(x, N, rth, win, beta);

%% rth = 0.005
t = 0:0.01:0.99;    % 100 time instants in the interval [0,1)
x = sin(10*t.^2) + 2*t;        % signal samples
N = 20; rth = 0.005;
win = 3; beta = 15;
[y, ra] = mdctcompr(x, N, rth, win, beta);
close all
plot(t,x,'r-',t,y,'b'); title('2N = 40, r_{thr} 0.005, r_a = 0.225')
xlabel('t'), legend('original', 'compressed'), grid on

%% rth = 0.05
t = 0:0.01:0.99;    % 100 time instants in the interval [0,1)
x = sin(10*t.^2) + 2*t;        % signal samples
N = 20; rth = 0.05;
win = 3; beta = 15;
[y, ra] = mdctcompr(x, N, rth, win, beta);
close all
plot(t,x,'r-',t,y,'b'); title('2N = 40, r_{thr} 0.050, r_a = 0.150')
xlabel('t'), legend('original', 'compressed'), grid on

%% TDAC Illustration
x = (1:28)';    % Example Input Signal
N = 4;          % overlap
beta = 6;       % KBD Parameter
win = 3;        % KBD window


X = buffer(x, 2*N, N, 'nodelay');
M = size(X,2);
w = (pbwin(N, win, beta))';      
W = repmat(w,1,M);
D = mdct(W.*X);
Y = W.*imdct(D);

% generate values to print table
y1 = [Y(:,1); zeros(20,1)];
y2 = [zeros(4,1);Y(:,2);zeros(16,1)];
y3 = [zeros(8,1);Y(:,3);zeros(12,1)];
y4 = [zeros(12,1);Y(:,4);zeros(8,1)];
y5 = [zeros(16,1);Y(:,5);zeros(4,1)];
y6 = [zeros(20,1);Y(:,6)];
y = y1 + y2 + y3 + y4 + y5 + y6;


fprintf('   x     y1       y2       y3       y4       y5       y6       y\n');
fprintf('--------------------------------------------------------------------\n');
fprintf('%3.0f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f \n',[x';y1';y2';y3';y4';y5';y6';y']);
