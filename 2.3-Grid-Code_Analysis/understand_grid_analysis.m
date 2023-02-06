% We want to analyze whether people have grid-like activity in mPFC
% That is, a sixfold symmetric activity pattern
figure;
plot(sind(6*(1:360))); ylabel('Activity'); xlabel('Movement in angle');
xlim([0 360])
% The highest neural activity would be at the 6 peaks of this plot

%% The 'issue' is that every person has a different 'angle' or 'phase'
% This phase needs to be estimated first (first half of the data)

% How the shift works:
% So you 'push-back' the function by advancing the degrees by the shift
% e.g., When you shift by 30° you calculate sin(31) instead of sin(1) etc..
fh = figure;
fh.WindowState = 'Maximized';
% Sine
subplot(2,1,1); hold on
plot(sind((1:360)),'linewidth',2)
for iShift = 1:5
    plot(sind((1:360)+ 30*iShift))
end
xlim([1 360])
l = legend( string((0:5)*30));
title(l,'Shift by:'); title('Sine')

% Cosine
subplot(2,1,2); hold on
plot(cosd((1:360)),'linewidth',2)
for iShift = 1:5
    plot(cosd((1:360)+ 30*iShift))
end
xlim([1 360])
l = legend( string((0:5)*30));
title(l,'Shift by:'); title('Cosine')
suptitle('Phase shift for Sine and Cosine')

% So in order to have the peak at the correct angle you have to subtract
% from Cosine (Because Cosine starts at the peak)

%% How is this shift estimated from the neural data?
% This is done by looking at how the neural activity reacts to the 
% Cosine and Sine of the angles made during the experiment (i.e., their betas)
% The combination of these betas can give us the main grid angle

% Let's create some random beta's for Cos and Sin (10 participants) [-1 1]
betas   = 2*rand(10,2)-1;
% And lets see the difference between atan and atan2
gridAng = atand(betas(:,2) ./ betas(:,1));
gridAngle = atan2d(betas(:,2), betas(:,1));
disp('    Beta Cos  Beta Sin Grid Angle')
disp([betas gridAng])

% In a experiment where we care about 0-360 this would be the peak of activity
% Regular atan
fh = figure; fh.WindowState = 'Maximized'; subplot(2,1,1); hold on
for iP = 1:5
    plot(cosd( (1:360) - gridAng(iP)), 'linewidth', 2)
end
l = legend( string(gridAng(1:5)));
title(l,'Peak at:')
title('atand')

% atan2
subplot(2,1,2); hold on
for iP = 1:5
    plot(cosd( (1:360) - gridAngle(iP)), 'linewidth', 2)
end
l = legend( string(gridAngle(1:5)));
title(l,'Peak at:')
title('atan2d')
txt = 'When their input is the same \newline';
txt = [txt 'There is a difference between atan and atan2 \newline'];
txt = [txt 'Because atan is bounded [-90, 90] and atan2 [-180, 180]'];
suptitle(txt)
% So basically the combination of Sin and Cosine determines the Main Grid
% Angle (phi)

% The difference between atan and atan2 matters because their answers can
% differ by 180°. Or 30° when divided by 6.
%% How does this work for a 6-fold signal
% Will only use atan2d from now on
gridAngle = (atan2d(betas(:,2), betas(:,1)))/6;
disp('    Beta Cos  Beta Sin Grid Angle')
disp([betas gridAngle])

fh = figure; fh.WindowState = 'Maximized'; hold on
plot(cosd( 6*(1:360) ), 'linewidth', 3, 'color','k')
for iP = 1:5
    plot(cosd( 6*((1:360) - gridAngle(iP)) ), 'linewidth', 2)
end
l = legend( ['Baseline'; string(gridAngle(1:5))]);
title(l,'Peak at:')
title('atan2d')

%% Once more for fake participants
% Let's say we have participants with known grid-angles 15, 30, and 45°
% The Parametric Modulators for Sine and Cosine reflect how much each of
% these contribute to the potential grid-angle.

% The 'complexity' is the 6-fold symmetry (which means you can only have a
% maximum offset of 59°). This is achieved by using the following PM's
% Cos(6*theta(t)) & Sin(6*theta(t)), which compresses the cos and sin waves
% into having 6 peaks.

% But why can you still get values bigger than this 59°?

% figure
% plot(cosd(6*(1:360))*-1); hold on
% plot(sind(1:360))

y = repmat([zeros(20,1); ones(20,1)],3,1)*3 + randn(120,1);
X = repmat([zeros(20,1); ones(20,1)],3,1); % This would be my regressor
X = [X, randn(120,1)]; % Second one just noise (to check what would happen)
% This is Matlab's GLM
mdl = fitglm(X,y);
% This is just in code
glmCode = @(X,y) ((X'*X)^-1)*X'*y;
betas = glmCode(X,y);
% Both give the same beta's

%% Now how do I add a Parametric Modulator to this?
y = repelem([3 0 -1 0 2 0 -1],20)';
y = y + randn(length(y),1)*.5;
xInt = ones(length(y),1);
xReg = repelem([1 0 1 0 1 0 1],20)';
xPM1 = repelem([3 0 -1 0 2 0 -1],20)';
xPM2 = repelem(randi([0 5],1,7),20)';
X = [xInt, xReg, xReg.*xPM1 xReg.*xPM2];

beta = glmCode(X,y);
% mdl = fitglm(X,y)

%% Create a single voxel's output to stimuli of 1-360° with a phase-shift
glmCode  = @(X,y) ((X'*X)^-1)*X'*y; % To calculate the beta's
phShift  = -8; % Change this to change participants 'Angle'
noiseLvl = .3;
expAng   = 1:360;

% -== Create Single Voxel's response ==-
y = cosd( 6*(expAng - phShift)) + randn(1,360)*noiseLvl;
y2 = cosd( 6*(expAng - phShift)) + randn(1,360)*noiseLvl;
fh = figure; fh.WindowState = 'Maximized';

% ######################## Plotting #######################################
% Stimulus Response
subplot(4,5,[6 7])
plot(y, 'LineWidth', 2);
title(['Neural Response to stimuli with angles [1-360°], phase shift: ' num2str(phShift) '°'])
ylabel('Magnitude Neural Activity'); xlabel('Trial/ Stimulus Angle'); xlim([1 360])
text(100, max(y)*3, 'GLM ONE', 'clipping', 'off', 'FontSize', 50)
% Design Matrix Stimulus Response
subplot(4,5,[3 8])
imagesc(y'); colormap('gray'); xticks(1); xticklabels('Neural Activity');
set(gca,'xaxisLocation','top'); ylabel('Trial/ Stimulus Angle')
% #########################################################################

% -== Create the design matrix ==-
xReg = ones(360,1);
xPM1 = cosd(6*expAng)';
xPM2 = sind(6*expAng)';
X    = [xReg, xReg.*xPM1, xReg.*xPM2];

% ######################## Plotting #######################################
% Design Matrix
subplot(4,5,[4 5 9 10]); imagesc(X); colormap('gray'); ylabel('Trial')
xticks(1:3); xticklabels({'Regressor','PM cos','PM sin'}); set(gca,'xaxisLocation','top')
% #########################################################################

% -== Calculate the beta's ==- (Could also use MATLAB's fitglm() )
betas = glmCode(X,y');
disp('Recovered betas:');
disp(['Cosine: ' num2str(betas(2))]);
disp(['Sine  : ' num2str(betas(3))]);
disp('')

% -== Calculate the main Grid Angle ==-
gridAngle = atan2d(betas(3), betas(2))/6;
disp(['Calculated Grid Angle: ' num2str(gridAngle)])
disp(['Actual Grid Angle: ' num2str(phShift)])
disp('')

% -== Calculate the 2nd GLM ==- where we check the consistency of the Grid Angle
% Uses a 'different' set of reponses
xPM3 = cosd(6* (expAng - gridAngle))';
X    = [xReg, xReg.*xPM3];
% x2 should be fit significantly (check output)
mdl = fitglm(X,y2') %#ok<NOPTS>

% ######################## Plotting #######################################
% Second Data set: Neural Response + Parametric Modulator
subplot(4,5,[16 17])
plot(y2, 'LineWidth',2)
hold on; plot(xPM3, 'LineWidth',2); xlim([1 360]);
legend('Neural Response','Fit Response (PMod)')
xlabel('Trial'); ylabel('Response')
title(['Calculated Grid Angle: ' num2str(gridAngle) ' vs Actual Grid Angle: ' num2str(phShift)])
text(100, max(y)*3, 'GLM TWO', 'clipping', 'off', 'FontSize', 50)

% Second Data set: Design Matrix
subplot(4,5,[13 18])
imagesc(y2'); colormap('gray'); xticks(1); xticklabels('Neural Activity');
set(gca,'xaxisLocation','top'); ylabel('Trial/ Stimulus Angle')
subplot(4,5,[14 15 19 20])
imagesc(X); colormap('gray'); ylabel('Trial/ Stimulus Angle')
xticks(1:2); xticklabels({'Regressor','PM GA'}); set(gca,'xaxisLocation','top')
% #########################################################################


%% Needs variables from above section to function
% Find the aligned & misaligned trials (first do for 6 each)
gaAlign = [(phShift - 15) + (0:30:390); (phShift + 15) + (0:30:390)];
allTrials = cell(size(gaAlign,2),1);

for iTrial = 1:size(expAng,2)
    % Angle and Onset are the same in this case
    trial = expAng(iTrial);
    idx = find(trial(1) >= gaAlign(1,:) & trial(1) < gaAlign(2,:));
    if ~isempty(idx)
        allTrials{idx,1} = [allTrials{idx,1}; trial];
    end
end

if ~isempty(allTrials{13,1})
    allTrials{1,1} = [allTrials{1,1}; allTrials{13,1}];
end
if ~isempty(allTrials{14,1})
    allTrials{2,1} = [allTrials{2,1}; allTrials{14,1}];
end

%% Create and calculate the regressors and pmods for aligned vs misaligned
fh = figure;
fh.WindowState = 'Maximized';
pause(.5)
% -== Create the design matrix ==-
dm = zeros(size(expAng,2),24);
dmIdx = 1:2:24;
xIdx  = string(0:30:330);
xlab  = cell(1,24);
for iT = 1:12
    % Regressor
    dm(allTrials{iT,1},dmIdx(iT)) = ones(size(allTrials{iT,1},1),1);
    xlab{dmIdx(iT)} = ['Reg ' num2str(iT)];
    % Parametric Modulator
    dm(allTrials{iT,1},dmIdx(iT)+1) = cosd(6*(allTrials{iT,1}-phShift));
    xlab{dmIdx(iT)+1} = ['pMod: ' xIdx{iT}];
end
subplot(4,4,5)
text(-20,175,'Regress & PMOD', 'FontSize', 15)
imagesc(dm); colormap('gray')
% xticks(1:24); xticklabels(xlab); xtickangle(45); set(gca, 'xaxislocation','top')

% Plot and show differences
cols = repmat([234 69 14; 14 69 234]/256,6,1);
betas = glmCode(dm,y2');

% Parametric Modulators & Regressor
subplot(4,4,6)
b = bar(betas); xticklabels(string(0:30:330)); ylabel('Beta'); title('Regressor & Parametric Modulator')
b.FaceColor = 'flat';
b.CData(2:2:24,:) = cols;
b.CData(1:2:24,:) = repmat([80 80 80]/256,12,1);
text(-10,1,'Regressor','Color',[80 80 80]/256,'FontWeight','Bold')
text(-10,.5,'Aligned','Color',[234 69 14]/256,'FontWeight','Bold')
text(-10,0,'Misaligned','Color',[14 69 234]/256,'FontWeight','Bold')
% Regressor
% figure; b = bar(betas(1:2:24)); xticklabels(string(0:30:330)); ylabel('Beta'); title('Regressor')
% b.FaceColor = 'flat';
% b.CData = cols;

%% And if we do just the regressors??
dm = zeros(size(expAng,2),12);
dmIdx = 1:2:24;
xIdx  = string(0:30:330);
xlab  = cell(1,12);
for iT = 1:12
    % Regressor
    dm(allTrials{iT,1},iT) = ones(size(allTrials{iT,1},1),1);
    xlab{iT} = ['Reg ' num2str(iT)];
end
subplot(4,4,1)
text(2,-175,'Design Matrix', 'FontSize', 15)
text(-7,175,'Regressor', 'FontSize', 15)
imagesc(dm); colormap('gray')
% xticks(1:24); xticklabels(xlab); xtickangle(45); set(gca, 'xaxislocation','top')
% Calculate the betas
betas = glmCode(dm,y2');

% Regressor
subplot(4,4,2)
text(2,2,'Resulting Betas', 'FontSize', 15)
b = bar(betas); xticklabels(string(0:30:330)); ylabel('Beta'); title('Regressor')
b.FaceColor = 'flat';
b.CData = cols;
text(-5,.25,'Aligned','Color',[234 69 14]/256,'FontWeight','Bold')
text(-5,-.25,'Misaligned','Color',[14 69 234]/256,'FontWeight','Bold')

%% Orthogonalize the Parametric Modulators
dm = zeros(size(expAng,2),24);
dmIdx = 1:2:24;
xIdx  = string(0:30:330);
xlab  = cell(1,24);
for iT = 1:12
    % Regressor
    re = ones(size(allTrials{iT,1},1),1);
    dm(allTrials{iT,1},dmIdx(iT)) = ones(size(allTrials{iT,1},1),1);
    xlab{dmIdx(iT)} = ['Reg ' num2str(iT)];
    % Parametric Modulator
    pm = cosd(6*(allTrials{iT,1}-phShift));
    dm(allTrials{iT,1},dmIdx(iT)+1) = pm - dot(pm,re)*re/norm(re.^2);
    xlab{dmIdx(iT)+1} = ['pMod: ' xIdx{iT}];
end
subplot(4,4,9)
text(-20,175,'Regress & Orth PMOD', 'FontSize', 15)
imagesc(dm); colormap('gray')
% xticks(1:24); xticklabels(xlab); xtickangle(45); set(gca, 'xaxislocation','top')

% Plot and show differences
cols = repmat([234 69 14; 14 69 234]/256,6,1);
betas = glmCode(dm,y2');

% Parametric Modulators
subplot(4,4,10)
% b = bar(betas(2:2:24)); xticklabels(string(0:30:330)); ylabel('Beta'); title('Parametric Modulator: Orthogonalized')
% b.FaceColor = 'flat';
% b.CData = cols;
% Regressor
% figure; b = bar(betas(1:2:24)); xticklabels(string(0:30:330)); ylabel('Beta'); title('Regressor: PM orthogonalized')
% b.FaceColor = 'flat';
% b.CData = cols;

b = bar(betas); xticklabels(string(0:30:330)); ylabel('Beta'); title('Regressor and Orthogonalized Parametric Modulator')
b.FaceColor = 'flat';
b.CData(2:2:24,:) = cols;
b.CData(1:2:24,:) = repmat([80 80 80]/256,12,1);
text(-10,2.5,'Regressor','Color',[80 80 80]/256,'FontWeight','Bold')
text(-10,0,'Aligned','Color',[234 69 14]/256,'FontWeight','Bold')
text(-10,-2.5,'Misaligned','Color',[14 69 234]/256,'FontWeight','Bold')

%% Orthogonalize the Regressor
dm = zeros(size(expAng,2),24);
dmIdx = 1:2:24;
xIdx  = string(0:30:330);
xlab  = cell(1,24);
for iT = 1:12
    % Regressor
    pm = cosd(6*(allTrials{iT,1}-phShift));
    re = ones(size(allTrials{iT,1},1),1);
    re_orth = re - dot(re,pm)*pm/norm(pm.^2);
    dm(allTrials{iT,1},dmIdx(iT)) = re_orth;
    xlab{dmIdx(iT)} = ['Reg ' num2str(iT)];
    % Parametric Modulator
    dm(allTrials{iT,1},dmIdx(iT)+1) = pm;
    xlab{dmIdx(iT)+1} = ['pMod: ' xIdx{iT}];
end
subplot(4,4,13)
text(-20,175,'Orth Regress & PMOD', 'FontSize', 15)
imagesc(dm); colormap('gray')
% xticks(1:24); xticklabels(xlab); xtickangle(45); set(gca, 'xaxislocation','top')

% Plot and show differences
cols = repmat([234 69 14; 14 69 234]/256,6,1);
betas = glmCode(dm,y2');

% Parametric Modulators
subplot(4,4,14)
% figure; b = bar(betas(2:2:24)); xticklabels(string(0:30:330)); ylabel('Beta'); title('Parametric Modulator: Orthogonalized Regressor')
% b.FaceColor = 'flat';
% b.CData = cols;
% Regressor
% figure; b = bar(betas(1:2:24)); xticklabels(string(0:30:330)); ylabel('Beta'); title('Regressor: Orthogonalized')
% b.FaceColor = 'flat';
% b.CData = cols;

b = bar(betas); xticklabels(string(0:30:330)); ylabel('Beta'); title('Orthogonalized Regressor & Parametric Modulator')
b.FaceColor = 'flat';
b.CData(2:2:24,:) = cols;
b.CData(1:2:24,:) = repmat([80 80 80]/256,12,1);
text(-10,1,'Regressor','Color',[80 80 80]/256,'FontWeight','Bold')
text(-10,.5,'Aligned','Color',[234 69 14]/256,'FontWeight','Bold')
text(-10,0,'Misaligned','Color',[14 69 234]/256,'FontWeight','Bold')

%% Now do it for the aligned vs misaligned still
gaAlign = [(phShift - 15) + (0:30:390); (phShift + 15) + (0:30:390)];
ali     = gaAlign(:,1:2:14);
mis     = gaAlign(:,2:2:14);
allTrials = cell(2,1);

for iTrial = 1:size(expAng,2)
    % Angle and Onset are the same in this case
    trial = expAng(iTrial);
    idxAli = find(trial(1) >= ali(1,:) & trial(1) < ali(2,:));
    % Shouldn't have overlap but just to be sure
%     idxMis = find(trial(1) >= mis(1,:) & trial(1) < mis(2,:));
    if ~isempty(idxAli)
        allTrials{1,1} = [allTrials{1,1}; trial];
    else
        allTrials{2,1} = [allTrials{2,1}; trial];
    end
end

%% Only the regressor
dm = zeros(size(expAng,2),2);

% Create the Design Matrix
dm(allTrials{1,1},1) = ones(length(allTrials{1,1}),1);
dm(allTrials{2,1},2) = ones(length(allTrials{1,1}),1);
subplot(4,4,3)
text(1,-150,'Design Matrix', 'FontSize', 15)
imagesc(dm); colormap('gray')
xticks(1:2); xticklabels({'Aligned','Misaligned'}); xtickangle(15); set(gca, 'xaxislocation','top')

% Plot and show differences
cols = [234 69 14; 14 69 234]/256;
betas = glmCode(dm,y2');

% Regressor
subplot(4,4,4)
text(.8,1,'Resulting Betas', 'FontSize', 15)
b = bar(betas); xticklabels({'Aligned','Misaligned'}); ylabel('Beta'); title('Regressor')
b.FaceColor = 'flat';
b.CData = cols;

%% Regressor and Parametric Modulator
dm = zeros(size(expAng,2),4);

% Create the Design Matrix
% Aligned
dm(allTrials{1,1},1) = ones(length(allTrials{1,1}),1);
dm(allTrials{1,1},2) = cosd(6*(allTrials{1,1}-phShift));
% Misaligned
dm(allTrials{2,1},3) = ones(length(allTrials{2,1}),1);
dm(allTrials{2,1},4) = cosd(6*(allTrials{2,1}-phShift));
subplot(4,4,7)
imagesc(dm); colormap('gray')
xticks(1:4); xticklabels({'Reg Aligned','pMod Aligned','Reg Misaligned','pMod Misaligned'}); xtickangle(15); set(gca, 'xaxislocation','top')

% Plot and show differences
cols = [234 69 14; 14 69 234]/256;
betas = glmCode(dm,y2');

% Regressor
subplot(4,4,8)
b = bar(betas); xticklabels({'Reg Aligned','pMod Aligned','Reg Misaligned','pMod Misaligned'})
ylabel('Beta'); title('Regressor & Pmod'); 
b.FaceColor = 'flat';  xtickangle(15)
b.CData(1:2:4,:) = repmat([80 80 80]/256,2,1);
b.CData(2:2:4,:) = cols;

%% Regressor and Orthogonalized Parametric Modulator
orthX =@(x,y) x - dot(x,y)*y/norm(y.^2);
% re_orth = re - dot(re,pm)*pm/norm(pm.^2);
dm = zeros(size(expAng,2),4);

% Create the Design Matrix
% Aligned
re = ones(length(allTrials{1,1}),1);
pm = cosd(6*(allTrials{1,1}-phShift));
dm(allTrials{1,1},1) = re;
dm(allTrials{1,1},2) = orthX(pm,re);
% Misaligned
re = ones(length(allTrials{2,1}),1);
pm = cosd(6*(allTrials{2,1}-phShift));
dm(allTrials{2,1},3) = re;
dm(allTrials{2,1},4) = orthX(pm,re);
subplot(4,4,11)
imagesc(dm); colormap('gray')
xticks(1:4); xticklabels({'Reg Aligned','pMod Aligned','Reg Misaligned','pMod Misaligned'}); xtickangle(15); set(gca, 'xaxislocation','top')

% Plot and show differences
cols = [234 69 14; 14 69 234]/256;
betas = glmCode(dm,y2');

% Regressor
subplot(4,4,12)
b = bar(betas); xticklabels({'Reg Aligned','pMod Aligned','Reg Misaligned','pMod Misaligned'})
ylabel('Beta'); title('Regressor & Orthogonalized Pmod'); 
b.FaceColor = 'flat';  xtickangle(15)
b.CData(1:2:4,:) = repmat([80 80 80]/256,2,1);
b.CData(2:2:4,:) = cols;

%% Orthogonalized Regressor and Parametric Modulator
orthX =@(x,y) x - dot(x,y)*y/norm(y.^2);
% re_orth = re - dot(re,pm)*pm/norm(pm.^2);
dm = zeros(size(expAng,2),4);

% Create the Design Matrix
% Aligned
re = ones(length(allTrials{1,1}),1);
pm = cosd(6*(allTrials{1,1}-phShift));
dm(allTrials{1,1},1) = orthX(re,pm);
dm(allTrials{1,1},2) = pm;
% Misaligned
re = ones(length(allTrials{2,1}),1);
pm = cosd(6*(allTrials{2,1}-phShift));
dm(allTrials{2,1},3) = orthX(re,pm);
dm(allTrials{2,1},4) = pm;
subplot(4,4,15)
imagesc(dm); colormap('gray')
xticks(1:4); xticklabels({'Reg Aligned','pMod Aligned','Reg Misaligned','pMod Misaligned'}); xtickangle(15); set(gca, 'xaxislocation','top')

% Plot and show differences
cols = [234 69 14; 14 69 234]/256;
betas = glmCode(dm,y2');

% Regressor
subplot(4,4,16)
b = bar(betas); xticklabels({'Reg Aligned','pMod Aligned','Reg Misaligned','pMod Misaligned'})
ylabel('Beta'); title('Orthogonalized Regressor & Pmod'); 
b.FaceColor = 'flat';  xtickangle(15)
b.CData(1:2:4,:) = repmat([80 80 80]/256,2,1);
b.CData(2:2:4,:) = cols;