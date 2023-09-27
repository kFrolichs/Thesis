%% Figure 01
% Shows the creation of a toy dataset
fh1 = figure;
fh1.Position = [365,271,875,827];

%% Create Toy Data
% Settings
ans_scale    = [1 10];
adjustment   = 1;
ans_adjusted = [ans_scale(1) + adjustment, ans_scale(2) - adjustment];
num_factors  = 5;
item_factors = 10;
num_par      = 100;
tot_items    = num_factors * item_factors;
num_part     = 30;
let_pos      = [-.2, 1.15];

%% Create Behavioral Data
% Get an average value for each factor
all_parts = zeros(50,100);
for iPar = 1:num_par
    factor_avg = randi(ans_adjusted, num_factors, 1);
    % Get all items with this average factor value
    all_items  = repelem(factor_avg, item_factors, 1);
    % Add a bit of variety to the items [-1, 1]
    all_items = all_items + randi([0,2],50,1)-1;
    % Add it to the participants
    all_parts(:,iPar) = all_items;
end

% Create some randomly distributed noise
noise_lvl  = 2;
all_noise  = randn(tot_items,num_par) * noise_lvl;
% Add the noise to the participants' answers
all_parts = all_parts + all_noise;
% Cut off above and below the answer scale
all_parts(all_parts < ans_scale(1)) = ans_scale(1);
all_parts(all_parts > ans_scale(2)) = ans_scale(2);

%% Conceptual model
% The Coarse Big-5 model assumes every item within the factor is the same
model_coarse   = repelem(eye(num_factors),item_factors,item_factors);
% The Fine model takes into account each items similarity to all other items
model_fine     = corr(all_parts');

% Turn them into dissimilarity models
modelRDMfine   = 1 - model_fine;
modelRDMcoarse = 1 - model_coarse;

%% Show the lower triangles for both models
% We only care about the lower triangle
getLow = tril(ones(tot_items));
fac_label = {'One','Two','Three','Four','Five'};

% Plot the similarity matrix
ax1 = subplot(3,3,1); imagesc(model_fine, 'AlphaData', getLow);
title('Similarity Matrix'); xticks(5:10:50); xticklabels(fac_label);
yticks(5:10:50); yticklabels(fac_label); xlabel('Factors'); ylabel('Factors')
text(let_pos(1), let_pos(2), 'a','FontSize',20,'units','normalized');
ax1.Box = 'off';

% Plot the Coarse Model set the upper triangle to transparent
ax2 = subplot(3,3,2); imagesc(modelRDMcoarse, 'AlphaData', getLow);
title('Coarse Model'); xticks(5:10:50); xticklabels(fac_label);
yticks(5:10:50); yticklabels(fac_label); xlabel('Factors'); ylabel('Factors')
text(let_pos(1), let_pos(2), 'b','FontSize',20,'units','normalized')
ax2.Box = 'off';

% plot the Fine Model
ax3 = subplot(3,3,3); imagesc(modelRDMfine, 'AlphaData', getLow);
title('Fine Model'); xticks(5:10:50); xticklabels(fac_label);
yticks(5:10:50); yticklabels(fac_label); xlabel('Factors'); ylabel('Factors')
text(let_pos(1), let_pos(2), 'c','FontSize',20,'units','normalized')
ax3.Box = 'off';

% Change the colors
colormap(ax1,magma); colormap(ax2, viridis); colormap(ax3, viridis);

%% Neural
% Settings
brain_size     = [20 30]; % Assume a small brain ;)
[neural_coarse, neural_fine] = deal(zeros([brain_size, tot_items, num_part]));
factor_select  = [1:10:50; 10:10:50]';
factor_lower   = logical(tril(ones(item_factors),-1));
total_voxel    = brain_size(1) * brain_size(2);
% Pick some clean neurons
number_neurons = 3;
pick_neuron    = randperm(total_voxel,number_neurons);
get_clean_neu  = zeros(tot_items, number_neurons);

% I will use the conceptual models we created to create "neural activity"
% patterns that will give us results that we expect
% Will look like we present each of the 'trait words' in a row
for iPar = 1:num_part
    count_item = 1;
    for iFac = 1:num_factors
        factor_avg = rand(brain_size);
        for iItm = 1:item_factors
            % Select a correlation value from the matrix
            % First select values from the right factor
            factor_values = model_fine(factor_select(iFac,1):factor_select(iFac,2), ...
                factor_select(iFac,1):factor_select(iFac,2));
            % Take the lower correlation values without the self-correlations
            all_correlation = factor_values(factor_lower);
            % Take a random value from this
            get_value = randi(length(all_correlation));
            % Multiply this random value times the factor_average "neural activity"
            activity_clean = factor_avg * all_correlation(get_value);
            neural_fine(:,:,count_item, iPar) = activity_clean;

            count_item = count_item + 1;
        end
    end
end
% Add some noise to the data
neural_fine_noise = neural_fine + randn(size(neural_fine))*.2;

%% Show some "neurons'" response magnitudes to the 50 stimuli
% Setting
[row_n, col_n] = ind2sub(brain_size,pick_neuron);

sample_neurons = zeros(tot_items, number_neurons);
for iNeu = 1:number_neurons
    sample_neurons(:,iNeu) = squeeze(neural_fine(row_n(iNeu), col_n(iNeu), :, 1)); 
end

% Plot Settings
tLen   = tot_items;
xPatch = [0 tLen tLen 0];
nPatch = 5;
sPatch = tLen/nPatch;
xPidx  = 0:sPatch:tLen;
yP     = [0 0 max(sample_neurons,[],'all') max(sample_neurons,[],'all')];

subplot(3,3,[4 5]);
plot(sample_neurons,'linewidth',2); hold on
% Add a patch for visibility
for iPat = 1:2:6
    xP = [xPidx(iPat) xPidx(iPat+1) xPidx(iPat+1) xPidx(iPat)];
    patch(xP,yP,'k', 'FaceAlpha', .1, 'EdgeColor','none')
end
% l = legend({'#1','#2','#3'},'location','northwestoutside');
l = legend({'#1','#2','#3'},'location','northeastoutside');
title(l,'Neurons')
xticks(5:10:50); xlabel('Items'); ylabel('Response Magnitude')
title('Example Neuronal Activity per trait word');
text(-.1, let_pos(2), 'd','FontSize',20,'units','normalized')
ylim([0 max(sample_neurons,[],'all')])
xticklabels({'Factor 1','Factor 2','Factor 3','Factor 4','Factor 5',})

%% Calculate Data RDMs
% GLM Settings
X       = [ones(tot_items,1), eye(tot_items)]; % Only add an intercept
glmCode = @(X,y) ((X'*X)^-1)*X'*y;
betas  = zeros(size(neural_fine_noise));

% Loop over the voxels and calculate the beta's in a standard GLM
for iPar = 1:num_part
    beta = zeros(tot_items);
    for iVoxX = 1:brain_size(1)
        for iVoxY = 1:brain_size(2)
             y = squeeze(neural_fine_noise(iVoxX,iVoxY,:,iPar));
             b = glmCode(X,y);
             betas(iVoxX,iVoxY,:,iPar) = b(2:end);
        end
    end
end

%% Calculate the data RDMs
dataRDMs = zeros(tot_items, tot_items, num_part);

for iPar = 1:num_part
    beta_reshape       = reshape(betas(:,:,:,iPar), total_voxel, tot_items);
    dataRDMs(:,:,iPar) = 1 - corr(beta_reshape,'rows','pairwise');
end

ax6 = subplot(3,3,6); imagesc(dataRDMs(:,:,1)); colormap(ax6, magma);
title('Example Data RDM'); xticks(5:10:50); xticklabels(fac_label);
yticks(5:10:50); yticklabels(fac_label); xlabel('Factors'); ylabel('Factors')
text(let_pos(1), let_pos(2), 'e','FontSize',20,'units','normalized')

%% Compare Model and Data RDMs
par_results = zeros(num_part,3); % First fine, Second coarse, Third random
lower_tri   = logical(tril(ones(tot_items),-1));
rand_model  = abs(randn(tot_items));

for iPar = 1:num_part
    % Take the lower triangle (without the identity)
    data = dataRDMs(:,:,iPar);
    % Fine Model
    par_results(iPar,1) = corr(data(lower_tri),modelRDMfine(lower_tri),'type','Spearman');
    % Coarse Model
    par_results(iPar,2) = corr(data(lower_tri),modelRDMcoarse(lower_tri),'type','Spearman');
    % Random Model
    par_results(iPar,3) = corr(data(lower_tri),rand_model(lower_tri),'type','Spearman');
end
% Fisher F-to-z on the coefficients
par_results_Z = atanh(par_results);

%% Helper function to cleanly plot a Boxplot together with the individual data points
% Random model
% subplot(3,3,7);
% plotBoxScat(par_results_Z(:,3),'colors',[0 .447 .741])
% % ylim([min(par_results_Z(:,2))*.95 max(par_results_Z(:,2))*1.05])
% title('Random Model'); xticks([]); ylabel('Correlation Coefficients')
% text(let_pos(1), let_pos(2), 'f','FontSize',20,'units','normalized')
% 
% % Coarse Model
% subplot(3,3,8);
% plotBoxScat(par_results_Z(:,2),'colors',[.850 .325 .098])
% % ylim([min(par_results_Z(:,2))*.95 max(par_results_Z(:,2))*1.05])
% title('Coarse Model'); xticks([]); ylabel('Correlation Coefficients')
% text(let_pos(1), let_pos(2), 'g','FontSize',20,'units','normalized')
% 
% % Fine Model
% subplot(3,3,9);
% plotBoxScat(par_results_Z(:,1),'colors',[.929 .694 .125])
% % ylim([min(par_results_Z(:,1))*.95 max(par_results_Z(:,1))*1.05])
% title('Fine Model'); xticks([]); ylabel('Correlation Coefficients')
% text(let_pos(1), let_pos(2), 'h','FontSize',20,'units','normalized')

subplot(3,3,[7 8 9])
plotBoxScat(par_results_Z(:,[3,2,1]))
xticklabels({'Random Model','Coarse Model','Fine Model'})
ylabel('Correlation Coefficients'); title('Results');
text(-.05,let_pos(2),'f','FontSize',20,'units','normalized')

%% Add arrows to the plot
% These coordinates are hard-coded. If the plot size changes these will be incorrect
arh1 = annotation('arrow',[0.630285714285717,0.697142857142857],[0.707496977025392,0.652962515114873]);
arh2 = annotation('arrow',[0.885142857142862,0.884571428571429],[0.677267230955259,0.633615477629988]);
arh3 = annotation('arrow',[0.690857142857148,0.622857142857143],[0.408827085852477,0.332527206771463]);
arh4 = annotation('arrow',[0.882857142857148,0.882285714285715],[0.374969770253928,0.331318016928658]);

arh1.LineWidth = 5; arh2.LineWidth = 5; arh3.LineWidth = 5; arh4.LineWidth = 5;

arh1.Color = [.850 .325 .098]; arh3.Color = [.850 .325 .098];
arh2.Color = [.929 .694 .125]; arh4.Color = [.929 .694 .125];

% Save the figure
% saveas(figure(1), 'C:\Users\frolichs\Documents\Koen\Projects\Thesis\Code\Thesis\2.2-Representational_Similarity_Analysis\Figure_01.png')