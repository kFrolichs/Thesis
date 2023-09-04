%% Figure 01
% Shows the creation of a toy dataset

%% Create Toy Data
% Settings
ans_scale    = [1 10];
adjustment   = 1;
ans_adjusted = [ans_scale(1) + adjustment, ans_scale(2) - adjustment];
num_factors  = 5;
item_factors = 10;
num_par      = 100;

%% Create Behavioral Data
% Get an average value for each factor
all_parts = zeros(50,100);
for iPar = 1:num_par
    factor_avg = randi(ans_adjusted, num_factors, 1);
    % Get all items with this average factor value
    all_items  = repelem(factor_avg, item_factors, 1);
    all_items = all_items + randi([0,2],50,1)-1;
    % Add it to the participants
    all_parts(:,iPar) = all_items;
end

% Create some randomly distributed noise
noise_lvl  = 1;
all_noise  = randn(num_factors*item_factors,num_par)*noise_lvl;
% Add the noise to the participants' answers
all_parts = all_parts + all_noise;
% Cut off above and below the answer scale
all_parts(all_parts < ans_scale(1)) = ans_scale(1);
all_parts(all_parts > ans_scale(2)) = ans_scale(2);

%% Conceptual model
% The Coarse Big-5 model assumes every item within the factor is the same
model_coarse = repelem(eye(num_factors),item_factors,item_factors);
% The Fine model takes into account each items similarity to all other items
model_fine   = corr(all_parts');

figure; subplot(1,2,1); imagesc(model_coarse); subplot(1,2,2); imagesc(model_fine)
colormap('gray')
%% Neural

%% Calculate Data RDMs

%% Compare Model and Data RDMs