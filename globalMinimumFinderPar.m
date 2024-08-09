clear;
clc;

format long

parpool('local', 12)

TU0 = importdata('R_Main+3C1+2C2.txt');
TU45 = importdata('R_C2.txt');
%TU90 = importdata('R_C2.txt');

nPoints = size(TU0, 1);

silenceWindow=[4.7e-4,4.84e-4];
windowWidth = silenceWindow(2)-silenceWindow(1);
powerUsed = 10; % variable
NumTrialPoints = 2^powerUsed;
nOptT = 1;
optTmhats = {'TU45'};
nameSim = strcat('f500khx_',num2str(nOptT),'TU45_power',num2str(powerUsed));
%initialValues = 'ZERO';n_run = 1;
initialValues = 'Random';n_run = 25;

xSols = zeros(n_run,nOptT,2);
fSols = zeros(n_run,1);
reducPs = zeros(n_run,nOptT);

rng default % For reproducibility
gs = GlobalSearch('Display', 'off', 'NumTrialPoints', NumTrialPoints);

optTmhat = zeros(nOptT, nPoints, 2);

for iHat = 1:nOptT
    switch optTmhats{iHat}
        case 'TU45'
            optTmhat(iHat,:,:) = TU45;
        case 'TU90'
            optTmhat(iHat,:,:) = TU90;
        otherwise
            error('TUmHat has not been provided')
    end
end

lb = repmat([-3, 0], nOptT, 1);
ub = repmat([3, silenceWindow(1)], nOptT, 1);

F0 = int_fun(TU0, repmat([0, 0], nOptT, 1), optTmhat, silenceWindow);

tic
parfor irun = 1:n_run
    fprintf("%d from %d \n",irun, n_run)
    switch initialValues
        case 'ZERO'
            x0 = repmat([0, 0], nOptT, 1);
        case 'Random'
            x0 = zeros(nOptT, 2);
            for i = 1:nOptT
                x0(i, 1) = 6*rand-3;
                x0(i, 2) = silenceWindow(1)*rand;
            end
        otherwise
                error('initialValues has not been provided')
    end

    func = @(x0) int_fun(TU0, x0, optTmhat, silenceWindow);


    problem = createOptimProblem('fmincon','x0',x0,...
        'objective',func,'lb', lb,'ub', ub);

    [xSols(irun,:,:), fSols(irun)] = runopt(gs, problem);
    
end

toc

p = gcp;
delete(p)
delete(gcp('nocreate'));

U0 = sqrt(F0);
Uopt = sqrt((fSols));
reducPs = abs(U0-Uopt)/U0*100;

[fmin, idx]=min(fSols);
fprintf("Minimum Potential Function: %e  \n", fmin)
for i=1:nOptT
    fprintf("Optimized x(%d): (%e, %e)  \n", i, xSols(idx,i,:))
end
fprintf("Reduced Percentage: %e  \n", reducPs(idx))

% save(strcat(nameSim,'.mat'), "fSols")
% save(strcat(nameSim,'.mat'), "xSols", '-append')

function [xSols, fSols] = runopt(gs, problem)
    [xSols, fSols] = run(gs, problem);
end