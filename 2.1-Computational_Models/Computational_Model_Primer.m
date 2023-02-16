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
fh1 = figure; subplot(2,2,1);
plot(expec_time); xlabel('timesteps'); ylabel('Expectation for reward')
title('Change in food expectation with continuous rewards')

text(-1,.5,'\fontsize{15}\bfA')

%% 2.1.1.2 Increasing the complexity of the models

%% 2.1.1.3 Applying the model to real data

%% 2.1.1.4 Model Checks