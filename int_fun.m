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







% % % 
% % % 
% % % 
% % % function F = int_fun(x0)
% % % 
% % % a = x0(1);
% % % t0 = x0(2);
% % % 
% % % global TU0 TU45
% % % 
% % % 
% % % TU45Shiftef = TU45;
% % % TU45Shiftef(:,1) = TU45(:,1) + t0;
% % % 
% % % timesData = 4.5e-4:1e-5:5.5e-4;
% % % 
% % % TU0AtTimesData = interp1(TU0(:,1), TU0(:,2), timesData);
% % % TU45AtTimesData = interp1(TU45Shiftef(:,1), TU45Shiftef(:,2), timesData);
% % % 
% % % % figure(1)
% % % % plot(timesData, TU0AtTimesData, 'b')
% % % % hold on
% % % % plot(timesData, a*TU22_5AtTimesData, 'r')
% % % 
% % % 
% % % 
% % % 
% % % F = trapz(timesData, (TU0AtTimesData + a*TU45AtTimesData).^2);
% % % 
% % % 
% % % end