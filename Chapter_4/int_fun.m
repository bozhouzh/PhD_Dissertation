function F = int_fun(T0mhat, x0, optTmhat, silenceWindow)

nOptT = size (x0, 1);
assert(numel(optTmhat) ~= nOptT,'The number are not compatible!')

timesData = silenceWindow(1):1e-7:silenceWindow(2);
UAtTimesData = interp1(T0mhat(:,1), T0mhat(:,2), timesData);

for iT = 1:nOptT
    a = x0(iT, 1);
    t0 = x0(iT, 2);
    nP=size(optTmhat, 2);
    TITShiftef=zeros(nP, 2);
    TITShiftef(:,:) = optTmhat(iT,:,:);
    TITShiftef(:,1) = TITShiftef(:,1) + t0;
    UAtTimesData = UAtTimesData + a*interp1(TITShiftef(:,1), TITShiftef(:,2), timesData);
end

F = trapz(timesData, UAtTimesData.^2);

end
