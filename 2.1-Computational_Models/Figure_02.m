%% 2.1.1.3 Applying the model to real data
close all
fh3 = figure; fh3.Position = [680,249,810,849];
colGrn = [.466 .674 .188];
colRed   = [.635 .078 .184];
%% Create some data for visualization
data = [3,8,9,10,6,6,2,9,5,3,9,8,9,3,7,7,2,5,3,8];
participant_data = data + sind(1:18:360)*3; % Add a trend for fun
participant_data(participant_data > 10) = 10;
participant_data(participant_data < 0)  = 0;

bad_fit_data = randi([1,10],1,20) + sind(180:18:539)*3; % Random data with opposite trend
bad_fit_data(bad_fit_data > 10) = 10;
bad_fit_data(bad_fit_data < 0) = 0;

good_fit_data = participant_data + randn(1,20)*3;
good_fit_data(good_fit_data > 10) = 10;
good_fit_data(good_fit_data < 0) = 0;

% Participant data with bad fitting data (A)
subplot(4,5,[1,2]); hold on
plot(participant_data,'k','linewidth',2); ylim([0, 11]); xlim([1,20]);
ylabel('Rating'); plot(bad_fit_data,'color',colRed,'linewidth',2);
title('Simulations with bad fit'); text(-1,14,'\bf\fontsize{18}a')
drawnow(); pause(.5);

% Participant data with good fitting data (B)
subplot(4,5,[3,4]); hold on
plot(participant_data,'k','linewidth',2); ylim([0, 11]); xlim([1,20]);
ylabel('Rating'); plot(good_fit_data,'color',colGrn,'linewidth',2);
title('Simulations with good fit'); text(-1,14,'\bf\fontsize{18}b')
drawnow(); pause(.5);

% Create the legend for A and B
legPlot = subplot(4,5,5); plot(1,nan,'k','linewidth',2); hold on;
plot(1,nan,'color',colRed,'linewidth',2); plot(1,nan,'color',colGrn,'linewidth',2);
set(legPlot,'Visible','off'); l = legend({'Participant', 'Bad Fit', 'Good Fit'});
title(l,'Simulated'); l.FontSize = 10;

% Participant data with bad fitting data (SSE) (C)
subplot(4,5,[6,7]); hold on
plot(participant_data,'k','linewidth',2); ylim([0, 11]); xlim([1,20]);
ylabel('Rating'); xlabel('Trial'); plot(bad_fit_data,'color',colRed,'linewidth',2);
% title('Participant data with bad fitting model'); 
text(-1,14,'\bf\fontsize{18}c')
plotSSE(participant_data, bad_fit_data);

% Participant data with good fitting data (SSE) (D)
subplot(4,5,[8,9]); hold on
plot(participant_data,'k','linewidth',2); ylim([0, 11]); xlim([1,20]);
ylabel('Rating'); xlabel('Trial'); plot(good_fit_data,'color',colGrn,'linewidth',2);
% title('Participant data with good fitting model');
text(-1,14,'\bf\fontsize{18}d')
plotSSE(participant_data, good_fit_data);

% Create the legend for C and D
legPlot = subplot(4,5,10); plot(1,nan,'k','linewidth',2); hold on;
plot(1,nan,'color',colRed,'linewidth',2); plot(1,nan,'color',colGrn,'linewidth',2);
plot(1,nan,'color',[0 .447 .741],'linewidth',2);
set(legPlot,'Visible','off'); l = legend({'Participant', 'Bad Fit', 'Good Fit','Abs. Diff'});
title(l,'Simulated'); l.FontSize = 10;

%% How parameters change the models 'behavior' (E)
subplot(4,5,11:14); hold on
calcPE   = @(expectation, outcome) outcome - expectation;
newModel = @(alpha, expectation, PE, gamma, RP) (gamma*RP) + ((1-gamma) * (expectation + alpha * PE));
alphaVal = [.05,.1,.2,.3];
gammaVal = [.5,.1,.9,.3];
RP       = data; % Assumes RP are exactly like data
expec    = zeros(4,21);

for iVal = 1:4
    for iCalc = 1:20
        PE = calcPE(expec(iCalc), data(iCalc));
        expec(iVal,iCalc+1) = newModel(alphaVal(iVal), expec(iCalc),PE,gammaVal(iVal),RP(iCalc));
    end
end
plot(data,'k:','linewidth',3); plot(expec','linewidth',2);
xlim([1,20]); ylim([0,11]);
title('Effects of different parameters on model behavior');
text(0.1,13,'\bf\fontsize{18}e')

% The legend for section E
legPlot2 = subplot(4,5,15); hold on
plot(1,nan,'k:','linewidth',3); plot(1,nan,'linewidth',2); plot(1,nan,'linewidth',2);
plot(1,nan,'linewidth',2); plot(1,nan,'linewidth',2);
set(legPlot2,'Visible','off'); 
l2 = legend('Participant',['\alpha: ' num2str(alphaVal(1)) ', \gamma: ' num2str(gammaVal(1))], ...
    ['\alpha: ' num2str(alphaVal(2)) ', \gamma: ' num2str(gammaVal(2))],...
    ['\alpha: ' num2str(alphaVal(3)) ', \gamma: ' num2str(gammaVal(3))],...
    ['\alpha: ' num2str(alphaVal(4)) ', \gamma: ' num2str(gammaVal(4))]);
title(l2, 'Parameters'); l2.FontSize = 10; l2.Position = [.75941878213,.344687010432,.171604933451,.131919902191];
drawnow();

%% Model Fitting at work (F)
% First plot the legend
legPlot3 = subplot(4,5,20); hold on
plotCol = [217,217,217;189,189,189;150,150,150;99,99,99;37,37,37]/256;
plot(1,nan,'k','linewidth',3); plot(1,nan,'color',plotCol(1,:),'linewidth',2);
set(legPlot3,'visible','off');
l3 = legend('Participant','Model w. param'); title(l3,'fminsearch');
l3.FontSize = 10; l3.Position = [.760653348707,.18744319238,.180246910084,.07008244815];

subplot(4,5,16:19); hold on
xlim([1,20]); ylim([0,11]);
text(0.2,13,'\bf\fontsize{18}f')

% This is the regular data + the RP (which is the data with some randomly distributed noise)
newData = [data; data+randn(1,20)*5];
param   = [.05, .8];
model   = @(param) fitModel(newData,param);
plot(data, 'k','linewidth',3);
title('fminsearch visualized')
global countPlot
countPlot = [];

% Have to rewrite our model slightly to accomodate for the fminsearch function
% Needs parameters in one matrix (see function below)
options = optimset('Display', 'final', 'MaxIter', 1e3);
[estimates, fval,exitflag,output] = fminsearch(model, param, options);
% Plotting happens inside the fitModel function

t1 = text(16,8,['\fontsize{12}\alpha: ' num2str(estimates(1))]);
t2 = text(16,10,['\fontsize{12}\gamma: ' num2str(estimates(2))]);

%% Functions used for this plotting
function plotSSE(data1, data2)
    thisSSE = 0;

    for iCheck = 1:20
        % Plot a line between each point
        plot([iCheck, iCheck],[data1(iCheck), data2(iCheck)],'color',[0 .447 .741],'linewidth',2)
        % Add textual information about each position and the accumulated SSE
        % Absolute Difference, Squared Difference, and SSE
        absDiff = round(abs(data1(iCheck)-data2(iCheck)),2);
        sqDiff  = round(abs(data1(iCheck)-data2(iCheck))^2,2);
        thisSSE  = thisSSE + sqDiff;

        t1 = text(20, 8, ['Abs Diff: ' num2str(absDiff)]);
        t2 = text(20, 6, ['Sqr Diff: ' num2str(sqDiff)]);
        t3 = text(20, 4, ['SSE: ' num2str(thisSSE)]);
        drawnow(); pause(.003);
        delete([t1,t2,t3])
    end
    title(['SSE: ' num2str(thisSSE)])
end

% This is the model in a function that returns the fit between data and model
function fitVal = fitModel(data, param)
    global countPlot
    plotCol = [217,217,217;189,189,189;150,150,150;99,99,99;37,37,37]/256; %#ok<NASGU>
    learnVec   = zeros(20,1);
    learnAndRP = zeros(1,20);
    
    for i = 1:20
        RP = data(1,i) - learnVec(i);
        learnVec(i+1) = learnVec(i) + param(1) * RP;
        learnAndRP(i) = param(2) * data(2,i) + (1-param(2)) * data(1,i);
    end
    fitVal = sum((data(1,:) - learnAndRP).^2);
    
    if param(1) <= eps || param(2) <= eps || param(1) >= 1 || param(2) >= 1
        fitVal = 10e6;
    end
    
    % Plot the last 5 results as the function is running
    countPlot = [countPlot; learnAndRP];
    for iPlot = 1:size(countPlot,1)
        eval(['plot' num2str(iPlot) ' = plot(countPlot(' num2str(iPlot) ',:),"color",plotCol(' num2str(iPlot) ',:),"linewidth",2);'])
    end
    t1 = text(16,8,['\fontsize{12}\alpha: ' num2str(param(1))]);
    t2 = text(16,10,['\fontsize{12}\gamma: ' num2str(param(2))]);
    drawnow(); %pause(.1);
    
    % Remove the lines
    for iPlot = 1:size(countPlot,1)
        eval(['delete(plot' num2str(iPlot) ');']);
    end
    % Remove the last line
    if size(countPlot,1) > 4
        countPlot(1,:) = [];
    end
    % Remove the text with current parameter values
    delete([t1,t2])
end
