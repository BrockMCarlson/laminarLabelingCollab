function CSD = calcCSD_classic(LFP)
% LFP is in (ch x ms x trial)


totchan = size(LFP,1);
maxtr = size(LFP,3);

% electrode contact spacing in mm:
d = .1; % This is the electrode spacing

%padarray would also work here for simplicity.
LFPvak(1,:,:) = LFP(1,:,:);
LFPvak(2:totchan+1,:,:) = LFP;
LFPvak(totchan+2,:,:) = LFP(end,:,:);
 
% Calculate CSD
%a -- chan above
%b -- f(x) voltage
%c -- chan below
for i = 2:totchan+1  %loop through channels of padded LFP (i.e 2-33)
    a = LFPvak(i-1,:,:); %output is (ms x trial)
    b = LFPvak(i,:,:); % This is your target channel
    c = LFPvak(i+1,:,:);
    chOut = i-1;
     CSD(chOut,:,:) = (a - 2*b + c)./(d^2);
end

%See page 3 schroeder et al (1998) for more detail.