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

% How parameters change the models 'behavior' (E)
subplot(4,5,11:14); hold on


% Model Fitting at work (F)
subplot(4,5,16:19); hold on

%% 2.1.1.4 Model Checks

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
        drawnow(); pause(.3);
        delete([t1,t2,t3])
    end
%     title({titleNam,['SSE: ' num2str(thisSSE)]})
    title(['SSE: ' num2str(thisSSE)])
end
