function CSD = calcCSD_classic(EVP)
% EVP is in (chan x ms) i.e. (32 x 301)


totchan = size(EVP,1);

% electrode contact spacing in mm:
d = .1; % This is the electrode spacing

%padarray would also work here for simplicity.
EVPvak = padarray(EVP,[1 0],NaN,'replicate');
 
% Calculate CSD
%a -- chan above
%b -- f(x) voltage
%c -- chan below
for i = 2:totchan+1  %loop through channels of padded LFP (i.e 2-33)
    a = EVPvak(i-1,:); %output is (ms x trial)
    b = EVPvak(i,:); % This is your target channel
    c = EVPvak(i+1,:);
    chOut = i-1;
     CSD(chOut,:) = (a - 2*b + c)./(d^2);
end

%See page 3 schroeder et al (1998) for more detail.