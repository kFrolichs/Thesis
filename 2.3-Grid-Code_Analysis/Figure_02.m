fh2 = figure;
fh2.Position = [653,210,941,585];

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%  @@@@@@@@@@@@@@@@@@@@@ Creation of the Toy data @@@@@@@@@@@@@@@@@@@@@@@@@
%  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% Creating the Toy Data
% Create a single voxel's output to stimuli of 1-360° with a phase-shift
phShift  = -8; % Change this to change participants 'Angle'
noiseLvl = .3;
expAng   = 1:360;

% -== Create Single Voxel's response ==-
y  = cosd( 6*(expAng - phShift)) + randn(1,360)*noiseLvl;
y2 = cosd( 6*(expAng - phShift)) + randn(1,360)*noiseLvl;
% fh1 = figure;
% fh1.Position = [680,233,997,865];

% Plotting
subplot(3,4,[1 2 3]);
plot(y, 'LineWidth', 2);
title('Idealized neural response');
text(150,1.75,['phase shift: ' num2str(phShift) '°'],'fontsize',18)
ylabel('Magnitude Neural Activity'); xlabel('Trial/ Stimulus Angle');
xlim([1 360]); ylim([-2,2])

% Add a patch for visibility
tLen   = length(expAng);
nPatch = 6;
sPatch = tLen/nPatch;
xPidx  = 0:sPatch:tLen;
yP     = [-2 -2 2 2];

for iPat = 1:2:nPatch
    xP = [xPidx(iPat) xPidx(iPat+1) xPidx(iPat+1) xPidx(iPat)];
    patch(xP,yP,'k', 'FaceAlpha', .1, 'EdgeColor','none')
end
text(let_pos(1), let_pos(2), 'a','FontSize',20,'units','normalized')

%% Design Matrix Stimulus Response
subplot(3,4,4);
imagesc(y'); colormap('gray'); xticks(1); xticklabels('Idealized Neural Response');
set(gca,'xaxisLocation','top'); ylabel('Trial/ Stimulus Angle'); title('Different View')
text(let_pos(1)*1.1, let_pos(2)*1.1, 'b','FontSize',20,'units','normalized')

%% Create the design matrix
glmCode = @(X,y) ((X'*X)^-1)*X'*y; % To calculate the beta's
xReg    = ones(360,1);
xPM1    = cosd(6*expAng)';
xPM2    = sind(6*expAng)';
X       = [xReg, xReg.*xPM1, xReg.*xPM2];

% Plotting
subplot(3,4,[5 6 9 10]);
imagesc(X); colormap('gray'); ylabel('Trial')
xticks(1:3); xticklabels({'Intercept','PM Cosine','PM Sine'}); %set(gca,'xaxisLocation','top')
title('Design Matrix: GLM 1')
text(let_pos(1), let_pos(2), 'c','FontSize',20,'units','normalized')

%% Calculate the beta's (Could also use MATLAB's fitglm() )
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

% Empty plot, only do text
ax1 = subplot(3,4,[7 8]);
plot(nan,nan);
set(ax1,'Visible','off');

text(0.2,.4,{['\color[rgb]{.85 .325 .098}Cosine \beta: ' num2str(betas(2))],...
    ['\color[rgb]{.929 .694 .125}Sine \beta: ' num2str(betas(3))],...
    ['\color[rgb]{0 .447 .741}Calculated \phi: ' num2str(gridAngle)]},'FontSize',14)

%% Show the phase shifted
dispNewAngle = cosd(6* (expAng - gridAngle))';
subplot(3,4,[11 12]);
plot(dispNewAngle,'linewidth',2)
xlabel('Trial'); ylabel('Expected Activity'); title('Expected activity pattern from combined \beta''s')
xlim([0 360]); xticks(0:60:360)
hold on
plot(cosd((1:360)*6)*0.6726, 'linewidth',2)
plot(sind((1:360)*6)*-.74879,'linewidth',2)

for iPat = 1:2:nPatch
    xP = [xPidx(iPat) xPidx(iPat+1) xPidx(iPat+1) xPidx(iPat)];
    patch(xP,yP,'k', 'FaceAlpha', .1, 'EdgeColor','none')
end

text(let_pos(1), let_pos(2), 'd','FontSize',20,'units','normalized')
