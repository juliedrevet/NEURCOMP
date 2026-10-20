function [xpeak_jk,ypeak_jk] = jackknife_peak(Y,xvec,imax)
% jackknife_peak  jackknifed latency and amplitude of the first peak of a time course
%   [xpeak_jk,ypeak_jk] = jackknife_peak(Y,xvec,imax)
%   for each participant left out, time (xpeak_jk) and value (ypeak_jk) of
%   the first peak after 0 of the mean time course of the other participants:
%   first local maximum between 0 and xvec(imax) reaching at least half of
%   the maximum on this range, refined by spline interpolation around the
%   peak sample
%
%   Y      participants x time samples
%   xvec   time of each sample (s)
%   imax   last sample of the search range (index in xvec): the range ends at tmax = xvec(imax)   [numel(xvec)]
%
%   xpeak_jk, ypeak_jk   participants x 1
%   jackknifed latency = mean(xpeak_jk); statistics on xpeak_jk or ypeak_jk
%   must be corrected for the jackknife (t/(n-1), F/(n-1)^2; Ulrich & Miller, 2001)
% [adapted by claude]

if nargin < 3
    imax = numel(xvec);
end

ndata = size(Y,1);
xpeak_jk = nan(ndata,1);
ypeak_jk = nan(ndata,1);
for idata = 1:ndata
    % leave one out
    yavg = mean(Y(setdiff(1:ndata,idata),:),1);
    [xpeak_jk(idata),ypeak_jk(idata)] = get_peak(xvec,yavg,imax);
end

end

function [xpeak,ypeak] = get_peak(xvec,yavg,imax)
% time and value of the first local maximum of yavg after 0 (up to
% xvec(imax)) reaching at least half of the maximum on this range, refined
% by spline interpolation around the peak sample
irange = find(xvec >= 0 & (1:numel(xvec)) <= imax);
y = yavg(irange);
ilocmax = find(y(2:end-1) > y(1:end-2) & y(2:end-1) >= y(3:end))+1; % local maxima
ilocmax = ilocmax(y(ilocmax) >= 0.5*max(y));
if isempty(ilocmax) % no local maximum: maximum on the range
    [~,ipeak] = max(y);
else
    ipeak = ilocmax(1);
end
imax = irange(ipeak);
fun = @(t)-interp1(xvec(irange),yavg(irange),t,'spline');
xpeak = fminbnd(fun,xvec(max(imax-1,irange(1))),xvec(min(imax+1,irange(end))));
ypeak = -fun(xpeak);
end
