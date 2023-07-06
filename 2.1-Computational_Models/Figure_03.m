%% Figure 03
% Confusion Matrix and Parameter Recovery
%% Simulating Data
numFac   = 5;
numTrial = 50;
numProf  = 4;
ansScale = [1 10];
noiseLvl = 1;
factors  = repelem((1:numFac)',numTrial,1);
% Create some data with function
data = sim_data(numFac, numTrial, ansScale, noiseLvl);
RP   = data + randn(numFac*numTrial,1)*0.1;
% Plot to see if it works
figure; plot(data); xlabel('Trials'); ylabel('Answers'); title('Simulated response')

%% Confusion Matrix



%% Parameter Recovery
nRuns = 100;
nMod  = 3;
setParam = randi([200,800],nRuns,2)/1000;
simData  = zeros(numFac*numTrial,nRuns,nMod);
recParam = zeros(nRuns,nMod+1);
options  = optimset('Display', 'final', 'MaxIter', 1e3);

for iRun = 1:nRuns
    [~, simData(:,iRun,1)] = model_01(data,setParam(iRun,1));
    model = @(param) model_01(simData(:,iRun,1),param);
    recParam(iRun,1) = fminsearch(model, .05, options);
    
    [~, simData(:,iRun,2)] = model_02(data,factors,setParam(iRun,1));
    model = @(param) model_02(simData(:,iRun,2),factors,param);
    recParam(iRun,2) = fminsearch(model, .05, options);
    
    [~, simData(:,iRun,3)] = model_03(data,RP,setParam(iRun,:));
    model = @(param) model_03(simData(:,iRun,3),RP,param);
    recParam(iRun,3:4) = fminsearch(model, [.05,.5], options);
end

figure;
for iMod = 1:nMod
    subplot(1,3,iMod);
    scatter(setParam(:,1),recParam(:,iMod));
    xlabel('Set Parameters'); ylabel('Recovered Parameters')
end

%% Simulating Data
% Function that creates some data
function data = sim_data(nFac,nTrial,scale,noise)
    data = [];
    for iFac = 1:nFac
        pickAvg = randi(scale);
        facData = repelem(pickAvg,nTrial,1) + randn(nTrial,1)*noise;
        data    = [data; facData]; %#ok<AGROW>
    end
    % Cut-off at the answer scales
    data(data < scale(1)) = scale(1);
    data(data > scale(2)) = scale(2);
end

%% Models
% Regular Rescorla-Wagner model
function [model_fit, modelData] = model_01(data, p)
    modelData = zeros(length(data), 1);
    % Run over the data
    for iD = 1:length(data)
        % Calculate the Prediction Error (PE)
        PE = data(iD) - modelData(iD);
        % Update the next estimate based on the PE and the alpha
        modelData(iD + 1) = modelData(iD) + p(1) * PE;
    end
    modelData = modelData(2:end);
    model_fit = calcFit(data,modelData);
    if p > 1 || p < 0
        model_fit = 10e4;
    end
end

% Factor Model
function [model_fit, modelData] = model_02(data, factor, p)
    modelData = zeros(length(data), 1);
    facCheck  = 1;
    resetVal  = 5;
    % Run over the data
    for iD = 1:length(data)
        % The first time there's a new factor give it the resetVal
        if factor(iD) == facCheck
            modelData(iD) = resetVal;
            facCheck = facCheck + 1;
        end
        % Calculate the Prediction Error (PE)
        PE = data(iD) - modelData(iD);
        % Update the next estimate based on the PE and the alpha
        modelData(iD + 1) = modelData(iD) + p(1) * PE;
    end
    modelData = modelData(2:end);
    model_fit = calcFit(data,modelData);
    if p > 1 || p < 0
        model_fit = 10e4;
    end
end

% Rescorla-Wagner model with the Reference Point (RP)
function [model_fit, modelData] = model_03(data, RP, p)
    modelData = zeros(length(data), 1);
    for iD = 1:length(data)
        % Calculate the Prediction Error (PE)
        PE = data(iD) - modelData(iD);
        % Update the next estimate based on the PE and the alpha
        modelData(iD + 1) = (p(2) * RP(iD)) + ((1-p(2)) * modelData(iD) + p(1) * PE);
    end
    modelData = modelData(2:end);
    model_fit = calcFit(data,modelData);
    if p(1) > 1 || p(1) < 0
        model_fit = 10e4;
    end
end

% Calculates the SSE
function fit = calcFit(data,modelData)
    fit = sum((data - modelData).^2);
end

% Calculate the BIC
function BIC = calcBIC(nTrials, nParam, fVal)
    BIC = (nTrials * log(fVal/nTrials)) + (nParam * log(nTrials));
end