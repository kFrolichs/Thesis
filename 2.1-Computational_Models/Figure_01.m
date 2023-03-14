% This code is used to generate all the figures in Chapter 2.1
% In an attempt to explain our computational modelling rationale
% ~Koen Frolichs
%% 2.1.1.1 The Rescorla-Wagner learning rule
% Setting up the simple conditioning example
alpha       = .05; % Learning rate
% expec_t1    = 0;   % The dog expects no food at the beginning
timesteps   = 10;  % The number of trials in the experiment
reward      = 1;   % The dog will always be rewarded
expec_time  = zeros(1, timesteps+1); % To keep track of the dog's expectation over time
                                     % Starts with no expectation i.e., 0
close all
fh1 = figure; fh1.Position = [680,172,966,926]; % If the figure is placed weirdly remove this Position

%% %%%% Part A %%%%
subplot(2,2,1);
xlim([-.5,2.5]); xticks(0:2); ylim([-.1,1.1]); yticks([-.1,0,1]); hold on
xlabel('timesteps'); ylabel('Expectation for reward')
title('Change in food expectation with continuous rewards')
text(-.75,1.25,'\fontsize{15}\bfA')
% First expectation
scatter(0,0,'k','filled'); text(-.4,-.05,'Expectation t=0'); drawnow(); pause(.5)
% Calculate the PE
plot([0,0],[0,1],'color', [0 0.4470 0.7410],'Linewidth',2); scatter(0,1,[],[0 0.4470 0.7410],'filled')
text(.1,.5,'PE = 1','color', [0 0.4470 0.7410]); drawnow(); pause(.5)
% Update for the next timestep
plot([0,1],[0,0.05],'k','Linewidth',2); yticks([-.1,0,0.05,1]); scatter(1,.05,'k','filled')
text(.7, 0.02,'Expectation t=1'); drawnow(); pause(.5)
% Calculate next PE
plot([1,1],[0.05,1],'color', [0 0.4470 0.7410],'Linewidth',2); scatter(1,1,[],[0 0.4470 0.7410],'filled')
text(1.1,.5,'PE = .95','color', [0 0.4470 0.7410]); drawnow(); pause(.5)
% Update a last time
plot([1,2],[.05,.0975],'k','Linewidth',2); yticks([-.1,0,.05,.0975,1]); scatter(2,.0975,'k','filled')
text(1.6, 0.07,'Expectation t=2'); drawnow(); pause(.5)

%% %%%% Part B %%%%
% Create two anonymous functions to do the calculations
calcPE  = @(expectation, outcome) outcome - expectation; % To calculate the Prediction Error
RW_rule = @(expectation, PE, learning_rate) expectation + PE * learning_rate; % To calculate the RW updating

% Run the simulation
for iT = 1:timesteps
    % Retreive the expectation at this timepoint
    expectation = expec_time(iT);
    % Calculate the PE
    PE = calcPE(expectation, reward);
    % Calculate the new expectation at the next timepoint
    expec_time(iT+1) = RW_rule(expectation, PE, alpha);
end
% Plot the expectations over time

subplot(2,2,2); xlim([1,11]); xticks(0:11); ylim([0,.5]); hold on
xlabel('timesteps'); ylabel('Expectation for reward')
title('Change in food expectation with continuous rewards')
text(-.5,.55,'\fontsize{15}\bfB')

% Animate the plot so you see the updating over time
pause(1)
for iPlot = 1:timesteps
    plot(expec_time(1:iPlot+1),'k','Linewidth',3);
    pause(.2); drawnow()
end

%% %%%% Part C %%%%
rew = [1 1 0 1 0 0 0 1 1 0 0];
xAx = 2:11;
subplot(2,2,3); hold on
xlim([1,11]); ylim([-.1,1.1]); text(-2,0,'No Reward'); text(-1.5,1,'Reward')
ylabel('Expectation for reward'); xlabel('timesteps'); xticks(1:11); xticklabels(0:10)
scatter(xAx,rew(1:10),'k','filled'); scatter(1,0,[],[0.8500 0.3250 0.0980],'filled')
expec = zeros(1,11);
text(-.1,1.2,'\fontsize{15}\bfC')

for iPlot = 1:timesteps+1
    expec(iPlot + 1) = expec(iPlot) + ((rew(iPlot) - expec(iPlot))*.3);
    plot(expec(1:iPlot),'k','Linewidth',3);
    drawnow(); pause(.35);
end

%% %%%% Part D %%%%
% This is largely a repeat from B just with some different learning rates
alphaVal = [0.05, 0.1, 0.2, 0.3, 0.4];
expec_alpha  = zeros(length(alphaVal), timesteps+1);
% Run the simulation
for iAlpha = 1:length(alphaVal)
    for iT = 1:timesteps
        % Retreive the expectation at this timepoint
        expectation = expec_alpha(iAlpha,iT);
        % Calculate the PE
        PE = calcPE(expectation, reward);
        % Calculate the new expectation at the next timepoint
        expec_alpha(iAlpha,iT+1) = RW_rule(expectation, PE, alphaVal(iAlpha));
    end
end

subplot(2,2,4); xlim([1,11]); xticks(0:11); ylim([0,1.05]); hold on
xlabel('timesteps'); ylabel('Expectation for reward')
text(-.5,1.15,'\fontsize{15}\bfD')

for iPlot = 1:length(alphaVal)
    plot(expec_alpha(iPlot,:),'Linewidth',3);
    title(['Learning rate: ' num2str(alphaVal(iPlot))])
    drawnow(); pause(.35);
end
title({'Change in food expectation with continuous rewards','for different learning rates'})
leg = legend(string(alphaVal),'Location', 'Northwest');
title(leg,"Alpha's")
