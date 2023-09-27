% Figure 02
close all
% Needs some variables from Figure 01.
% So run that first.
let_pos    = [-.1, 1.1];

fh1 = figure;
fh1.Position = [680,53,941,1045];

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
    if ~isempty(idxAli)
        allTrials{1,1} = [allTrials{1,1}; trial];
    else
        allTrials{2,1} = [allTrials{2,1}; trial];
    end
end
dm = zeros(size(expAng,2),2);

% Create the Design Matrix
dm(allTrials{1,1},1) = ones(length(allTrials{1,1}),1);
dm(allTrials{2,1},2) = ones(length(allTrials{1,1}),1);
subplot(3,3,1)
text(1,-150,'Design Matrix', 'FontSize', 15)
imagesc(dm); colormap('gray')
xticks(1:2); xticklabels({'Aligned','Misaligned'}); xtickangle(15);
set(gca, 'xaxislocation','top'); ylabel('Trials/Slices')
text(let_pos(1), let_pos(2), 'a','FontSize',20,'units','normalized')

% Show on neural data what happens
subplot(3,3,2); plot(y2, 'LineWidth',2)
hold on; xlim([1 360]);
xlabel('Trial'); ylabel('Response')
% Color the regressors on here
getRegressors = (-30:30:360)+8;
for iPat = 1:length(getRegressors)-1
    xP = [getRegressors(iPat) getRegressors(iPat+1) getRegressors(iPat+1) getRegressors(iPat)];
    if mod(iPat,2)
        patch(xP,yP,[1 0 0],'FaceAlpha', .2, 'EdgeColor','none')
    else
        patch(xP,yP,[0 0 1],'FaceAlpha', .2, 'EdgeColor','none')
    end
end
title('Trials per regressor')
text(let_pos(1), let_pos(2),'b','FontSize',20,'units','normalized')

% Plot and show differences
cols = [234 69 14; 14 69 234]/256;
betas = glmCode(dm,y2');
% Regressor
subplot(3,3,3)
text(.8,1,'Resulting Betas', 'FontSize', 15)
b = bar(betas); xticklabels({'Aligned','Misaligned'}); ylabel('Beta');
title('Two Regressors');
b.FaceColor = 'flat';
b.CData = cols;
text(let_pos(1), let_pos(2), 'c','FontSize',20,'units','normalized')

%% Another 2nd GLM that looks for all aligned vs misaligned trials
gaAlign = [(gridAngle - 15) + (0:30:390); (gridAngle + 15) + (0:30:390)];
allTrials = cell(size(gaAlign,2),1);

for iTrial = 1:length(expAng)
    trial = expAng(iTrial);
    idx = find(trial >= gaAlign(1,:) & trial < gaAlign(2,:));
    if ~isempty(idx)
        allTrials{idx,1} = [allTrials{idx,1}; trial];
    end
end

% Put the trials outside the regular 0-360 in the correct spot
if ~isempty(allTrials{13,1})
    allTrials{1,1} = [allTrials{1,1}; allTrials{13,1}];
end
if ~isempty(allTrials{14,1})
    allTrials{2,1} = [allTrials{2,1}; allTrials{14,1}];
end

cols = repmat([234 69 14; 14 69 234]/256,6,1);
dm = zeros(size(expAng,2),12);
dmIdx = 1:2:24;
xIdx  = string(0:30:330);
xlab  = cell(1,12);
for iT = 1:12
    % Regressor
    dm(allTrials{iT,1},iT) = ones(size(allTrials{iT,1},1),1);
    xlab{iT} = ['Reg ' num2str(iT)];
end
subplot(3,3,4)
imagesc(dm); colormap('gray')
text(let_pos(1), let_pos(2), 'd','FontSize',20,'units','normalized')
% Calculate the betas
betas = glmCode(dm,y2');

% Show on neural data what happens
subplot(3,3,5); plot(y2, 'LineWidth',2)
hold on; xlim([1 360]);
xlabel('Trial'); ylabel('Response')
% Color the regressors on here
getRegressors = (-30:30:360)+8;
for iPat = 1:length(getRegressors)-1
    xP = [getRegressors(iPat) getRegressors(iPat+1) getRegressors(iPat+1) getRegressors(iPat)];
    if mod(iPat,2)
        patch(xP,yP,[1 0 0],'FaceAlpha', .2, 'EdgeColor','k')
    else
        patch(xP,yP,[0 0 1],'FaceAlpha', .2, 'EdgeColor','k')
    end
end
title('Trials per regressor')
text(let_pos(1), let_pos(2),'e','FontSize',20,'units','normalized')

% Regressor
subplot(3,3,6)
text(2,2,'Resulting Betas', 'FontSize', 15)
b = bar(betas); xticklabels(string(0:60:360)); ylabel('Beta'); title('Twelve Regressor')
b.FaceColor = 'flat';
b.CData = cols;
xticks(0:2:12); xlabel('Angles')
text(let_pos(1), let_pos(2), 'f','FontSize',20,'units','normalized')

%% -== Calculate the 2nd GLM ==- where we check the consistency of the Grid Angle
% Uses a 'different' set of reponses
xPM3 = cosd(6* (expAng - gridAngle))';
X    = [xReg, xReg.*xPM3];
% x2 should be fit significantly (check output)
% mdl = fitglm(X,y2') %#ok<NOPTS>

betas = glmCode(X,y2');

% Plotting
% Second Data set: Design Matrix
subplot(3,3,7)
imagesc(X); colormap('gray'); ylabel('Trial/ Stimulus Angle')
xticks(1:2); xticklabels({'Intercept','PM Grid-Angle'}); set(gca,'xaxisLocation','top')
title('Design Matrix: GLM 2 - One PM')
text(let_pos(1), let_pos(2), 'g','FontSize',20,'units','normalized')

% Second Data set: Neural Response + Parametric Modulator
subplot(3,3,8); plot(y2, 'LineWidth',2)
hold on; plot(xPM3, 'LineWidth',2); xlim([1 360]);
xlabel('Trial'); ylabel('Response')
title('\color[rgb]{0 .447 .741}Neural \color[rgb]{.85 .325 .098}Fitted')
for iPat = 1:2:nPatch
    xP = [xPidx(iPat) xPidx(iPat+1) xPidx(iPat+1) xPidx(iPat)];
    patch(xP,yP,'k', 'FaceAlpha', .1, 'EdgeColor','none')
end
% legend('Neural Response','Fit Response (PMod)')
text(let_pos(1), let_pos(2), 'h','FontSize',20,'units','normalized')

subplot(3,3,9)
bar(betas(2)); ylabel('beta'); xticks([]); xlabel('')
title('One Parametric Modulator')
text(let_pos(1), let_pos(2), 'i','FontSize',20,'units','normalized')