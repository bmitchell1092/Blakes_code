function [sdf, sua, tm] = spk2sdf(spk,Fs,st)

% st and spk are in sampels
% Fs is sampeling frequency
% all outputs are at 1Hz
% MAC

if nargin < 3 || isempty(st)
    st = 0;
end

% Defin Kernal (k must sum to 1)
k    = jnm_kernel( 'psp', 20 );

% spk resolution is 1000kHz (1spk/ms)
r   = Fs / 1000; 
st  = round(st/r);
spk = unique(round(spk/r)) - st;

% setup time 
if ~any(spk < 0)
    tm = [0:max(spk)];
    sua = zeros(size(tm));
    sua(spk + 1) = 1;
else
    tm = [min(spk):max(spk)];
    sua = zeros(size(tm));
    sua(spk - min(spk) + 1) = 1;
end

% pad sua if needed (must be longer than kernal)
keep = ones(size(sua));
if length(sua) < length(k)
    padln = length(k) - length(sua) + 1;
    keep  = [keep zeros(1,padln)];
    sua   = [sua zeros(1,padln)];
end

% convolve and trim
sdf = doConv(sua,k) * 10^3; % spikes per second
sdf = sdf(logical(keep));
sua = sua(logical(keep));




%k = normpdf(fix([-3.5*binstd:binstd*3.5]),0,binstd); % GAUSSIAN KERNEL
%k = poisspdf(fix([0:binstd*3]),binstd); % Poisson Kernel, binstd sets mean and std