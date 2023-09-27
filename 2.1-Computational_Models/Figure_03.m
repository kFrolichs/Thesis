%% Figure 03
% Confusion Matrix and Parameter Recovery
%% Simulating Data
numFac   = 5;
numTrial = 50;
totTrial = numFac*numTrial;
ansScale = [1 10];
noiseLvl = .3;
factors  = repelem((1:numFac)',numTrial,1);
% Create some data with function
data = sim_data(numFac, numTrial, ansScale, noiseLvl);
RP   = data + randn(numFac*numTrial,1)*0.1;
% Plot to see if it works
figure; plot(data); xlabel('Trials'); ylabel('Answers'); title('Simulated response')

%% Parameter Recovery
nRuns = 100;
nMod  = 3;
setParam = randi([200,800],nRuns,2)/1000;
simData  = zeros(totTrial,nRuns,nMod);
recParam = zeros(nRuns,nMod+1);
options  = optimset('Display', 'final', 'MaxIter', 1e3);

for iRun = 1:nRuns
    [~, simData(:,iRun,1)] = model_01([data,data],setParam(iRun,1));
    sDat = simData(:,iRun,1) + randn(totTrial,1)*noiseLvl;
    model = @(param) model_01([sDat,data],param);
    recParam(iRun,1) = fminsearch(model, .05, options);
    
    [~, simData(:,iRun,2)] = model_02([data,data],factors,setParam(iRun,1));
    sDat = simData(:,iRun,2) + randn(totTrial,1)*noiseLvl;
    model = @(param) model_02([sDat,data],factors,param);
    recParam(iRun,2) = fminsearch(model, .05, options);
    
    [~, simData(:,iRun,3)] = model_03([data,data],RP,setParam(iRun,:));
    sDat = simData(:,iRun,3) + randn(totTrial,1)*noiseLvl;
    model = @(param) model_03([sDat,data],RP,param);
    recParam(iRun,3:4) = fminsearch(model, [.05,.5], options);
end

%% Plot the results for the Parameter Recovery
fh1     = figure; fh1.Position = [375,528,593,571];
let_pos = [-.1, 1.1];
let     = {'a','b','c'};

for iMod = 1:nMod
    subplot(2,3,iMod);
    scatter(setParam(:,1),recParam(:,iMod));
    xlabel('Set Parameters'); ylabel('Recovered Parameters')
    title(['Model 0' num2str(iMod)]); ylim([0 1])
    text(.3,.9,['\rho: ' num2str(corr(setParam(:,1),recParam(:,iMod)))])
    
    text(let_pos(1), let_pos(2), let{iMod},'FontSize',20,'units','normalized')
end
sgtitle('Parameter Recovery')

%% Confusion Matrix
getBIC = zeros(nRuns,nMod,nMod);
getMin = zeros(nRuns,nMod);

for iMod = 1:nMod
    for iRun = 1:nRuns
        % Load the simulated data from each respective model
        sDat = simData(:,iRun,iMod) + randn(totTrial,1)*noiseLvl;
        
        % Run the simulated data on each model
        % Model 01
        model = @(param) model_01([sDat,data],param);
        [~,fVal] = fminsearch(model, .05, options);
        getBIC(iRun,1,iMod) = calcBIC(totTrial, 1, fVal);
        % Model 02
        model = @(param) model_02([sDat,data],factors,param);
        [~,fVal] = fminsearch(model, .05, options);
        getBIC(iRun,2,iMod) = calcBIC(totTrial, 1, fVal);
        % Model 03
        model = @(param) model_03([sDat,data],RP,param);
        [~,fVal] = fminsearch(model, [.05, .5], options);
        getBIC(iRun,3,iMod) = calcBIC(totTrial, 2, fVal);
        
        % Determine what model got the best BIC score
        [~,getMin(iRun,iMod)] = min(getBIC(iRun,:,iMod));
    end
end

% Calculate the Confusion Matric
CM = zeros(nMod);
for iRow = 1:nMod
    for iCol = 1:nMod
        CM(iRow,iCol) = sum(getMin(:,iRow) == iCol)/nRuns;
    end
end

%% Plot the results for the Confusion Matrix
subplot(2,3,[4,5])
imagesc(CM); colormap('gray'); colorbar;
xticks(1:3); yticks(1:3); xticklabels({'Model 01','Model 02','Model 03'})
yticklabels({'Model 01','Model 02','Model 03'});
ylabel('Simulated'); xlabel('Recovered'); title('Confusion Matrix')

text(let_pos(1), let_pos(2), 'd','FontSize',20,'units','normalized')

%% Save the Figure
saveas(figure(2), 'C:\Users\frolichs\Documents\Koen\Projects\Thesis\Figures\Chapter_02\Figure_03.png')
cropPlot('C:\Users\frolichs\Documents\Koen\Projects\Thesis\Figures\Chapter_02\Figure_03.png')
%% Functions used within this file
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
        PE = data(iD,2) - modelData(iD);
        % Update the next estimate based on the PE and the alpha
        modelData(iD + 1) = modelData(iD) + p(1) * PE;
    end
    modelData = modelData(2:end);
    model_fit = calcFit(data(:,1),modelData);
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
        PE = data(iD,2) - modelData(iD);
        % Update the next estimate based on the PE and the alpha
        modelData(iD + 1) = modelData(iD) + p(1) * PE;
    end
    modelData = modelData(2:end);
    model_fit = calcFit(data(:,1),modelData);
    if p > 1 || p < 0
        model_fit = 10e4;
    end
end

% Rescorla-Wagner model with the Reference Point (RP)
function [model_fit, modelData] = model_03(data, RP, p)
    modelData = zeros(length(data), 1);
    for iD = 1:length(data)
        % Calculate the Prediction Error (PE)
        PE = data(iD,2) - modelData(iD);
        % Update the next estimate based on the PE and the alpha
        modelData(iD + 1) = (p(2) * RP(iD)) + ((1-p(2)) * modelData(iD) + p(1) * PE);
    end
    modelData = modelData(2:end);
    model_fit = calcFit(data(:,1),modelData);
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