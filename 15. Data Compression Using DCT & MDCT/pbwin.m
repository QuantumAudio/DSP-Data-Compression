function w = pbwin(N, type, beta)
% Generates Princen-Bradley windows
% type = 0 generates a rectangular window
% type = 1 generates a sine window
% type = 2 generates a vorbis window
% type = 3 generates a KBD window
% N = DCT block length
% beta = Kaiser shape parameter (for KBD window only)

n = 0:2*N-1;

if type == 0
    w = ones([1,2*N]);    % Generate a rectangular window
end

if type == 1
    w = sqrt(2)*sin((pi/(2*N))*(n + 0.5)); % Generate sine window
end

if type == 2
    w = sqrt(2)*sin((pi/2)*sin((pi/(2*N))*(n + 0.5)).^2); % Generate vorbis window
end


if type == 3
    k = 0:N;
    % Compute zeroth order Bessel function of the first kind
    fk = besseli(0, beta*sqrt(1-((k - N/2)/(N/2)).^2)); 

    S = sum(fk);

    for n = 0:N-1
        wf(n+1) = sqrt((2/S)*sum(fk(1:(n+1))));    % Compute forward half of window
    end
    
    w = [wf, fliplr(wf)];                       % Concatenate forward and reverse halves
end

end