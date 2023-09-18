% Figure 01
% Shows how the angle influences the response magnitude
% Creation of the toy data
% GLM 1: finding the grid-angle

%% Grid-Cells and angles

%% Create Hexagonally Spaced Gaussians
numSteps   = 50;
numRepeats = 6;
x          = linspace(-3,3,numSteps);
y          = x';
[X,Y]      = meshgrid(x,y);
let_pos    = [-.1, 1.1];

% Repeat X a number of times
X = repmat(X,1,numRepeats);
% shift a row half the steps of one to the right
Xshift = circshift(X,round(numSteps / 2),2);
% Stick 'm together and repeat to match column length
Xtot = repmat([X; Xshift],numRepeats/2,1);
% Create Y of same size
Ytot = repmat(Y,numRepeats,numRepeats);

z = exp(-(Xtot.^2 + Ytot.^2) / 8);

% Remove a bit of the plot so the center of a Gaussian is on [0, 0]
z(1:numSteps/2,:) = [];
z(:,1:numSteps/2) = [];
z = z(1:end-1,1:end-1);
grid = z;

%% Get the positions for animation
totDeg = 360;
steps = round(size(grid,1) / 2.3);
stepsTot = 1:steps;

[x, y] = deal(ones(steps,totDeg));

for iDeg = 1:totDeg
    x(:,iDeg) = cosd(iDeg) * stepsTot + 125; % +steps to center it
    y(:,iDeg) = sind(iDeg) * stepsTot + 150;
end

%% Calculate the total 'activation' for every angled line
angleVal = zeros(totDeg, steps);

for iDeg = 1:totDeg
%     for iLine = 1:gridSize
    for iLine = 1:steps
        xVal = round(x(iLine,iDeg));
        yVal = round(y(iLine,iDeg));
        if xVal == 0; xVal = 1; end
        if yVal == 0; yVal = 1; end
        angleVal(iDeg, iLine) = z(yVal,xVal);
    end
end
angleSum = sum(angleVal,2);
% Has to be reversed because the axis on the plots differ
angleSum = fliplr(angleSum');

%% Animate the plot
fh1 = figure;
fh1.Position = [680,53,941,1045];
subplot(5,4,[1 2 5 6]);
imagesc(grid); colormap(magma); hold on
yMins = [40, 90];

for iDeg = 1:totDeg
%     subplot(1,2,1)
    subplot(5,4,[1 2 5 6]);
    lp1 = line(x(:,iDeg,1),y(:,iDeg,1), 'linewidth', 3, 'color', 'k');
    axis([0 size(grid,1) 0 size(grid,1)]);
    title(['Current Angle: ',num2str(iDeg),'°']);
    set(gca,'XTick',[], 'YTick', [],'ZTick',[]) % Removes the ticks
    
%     subplot(1,2,2)
    subplot(5,4,[3 4 7 8]);
    plot(1:iDeg,angleSum(1:iDeg),'color','k','linewidth',3)
    axis([0 360 yMins])
    xlabel('Trial Angle')
    ylabel('Response Magnitude')
    set(gca, 'YTick', [])

    drawnow % Animate it
    delete(lp1)
end

text(let_pos(1), let_pos(2), 'b','FontSize',20,'units','normalized')
subplot(5,4,[1 2 5 6]);
text(let_pos(1), let_pos(2), 'a','FontSize',20,'units','normalized')

pause(1)
    
%% Finish with a static image for the thesis
% Add a patch in both the plots
angles = [-15 + (0:30:360); 15 + (0:30:360)];
angles(:,1) = [];
angles(2,end) = 360;

subplot(5,4,[1 2 5 6]); colormap('gray')
title('Aligned and Misaligned trials')
for iPatch = 1:12
    endpoints = angles(:,iPatch);
    xPos = x([1, end],endpoints);
    yPos = y([1, end],endpoints);
    
    if mod(iPatch,2) == 1
        col = [0 0 1];
    else
        col = [1 0 0];
    end
    patch([xPos(1,1) xPos(2,1) xPos(2,2) xPos(1,2)], [yPos(1,1) yPos(2,1) yPos(2,2) yPos(1,2)], col, 'FaceAlpha', .3, 'EdgeColor','none')
end
% Add the last blue one [0 - 15]
patch([125 244 240 126], [150 150 181 150], [1 0 0], 'FaceAlpha', .3, 'EdgeColor','none')

% Add the patch for the left figure
subplot(5,4,[3 4 7 8]);
title('Response magnitude based on trial angle')
angles = [15:30:360, 360, 0, 15];
for iPatch = 1:12
    if mod(iPatch,2) == 1
        col = [0 0 1];
    else
        col = [1 0 0];
    end
    patch([angles(iPatch) angles(iPatch+1) angles(iPatch+1) angles(iPatch)], [yMins(1) yMins(1) yMins(2) yMins(2)], col, ...
        'FaceAlpha', .2, 'EdgeColor','none')
end
patch([0 15 15 0], [yMins(1) yMins(1) yMins(2) yMins(2)], col, 'FaceAlpha', .2, 'EdgeColor','none')

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
subplot(5,4,[9 10 11]);
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
text(let_pos(1), let_pos(2), 'c','FontSize',20,'units','normalized')

%% Design Matrix Stimulus Response
subplot(5,4,12);
imagesc(y'); colormap('gray'); xticks(1); xticklabels('Idealized Neural Response');
set(gca,'xaxisLocation','top'); ylabel('Trial/ Stimulus Angle'); title('Different View')
text(let_pos(1)*1.1, let_pos(2)*1.1, 'd','FontSize',20,'units','normalized')

%% Create the design matrix
glmCode = @(X,y) ((X'*X)^-1)*X'*y; % To calculate the beta's
xReg    = ones(360,1);
xPM1    = cosd(6*expAng)';
xPM2    = sind(6*expAng)';
X       = [xReg, xReg.*xPM1, xReg.*xPM2];

% Plotting
subplot(5,4,[13 14 17 18]);
imagesc(X); colormap('gray'); ylabel('Trial')
xticks(1:3); xticklabels({'Intercept','PM Cosine','PM Sine'}); %set(gca,'xaxisLocation','top')
title('Design Matrix: GLM 1')
text(let_pos(1), let_pos(2), 'e','FontSize',20,'units','normalized')

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
ax1 = subplot(5,4,[15 16]);
plot(nan,nan);
set(ax1,'Visible','off');

text(0.2,.4,{['\color[rgb]{.85 .325 .098}Cosine \beta: ' num2str(betas(2))],...
    ['\color[rgb]{.929 .694 .125}Sine \beta: ' num2str(betas(3))],...
    ['\color[rgb]{0 .447 .741}Calculated \phi: ' num2str(gridAngle)]},'FontSize',14)

%% Show the phase shifted
dispNewAngle = cosd(6* (expAng - gridAngle))';
subplot(5,4,[19 20]);
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

text(let_pos(1), let_pos(2), 'f','FontSize',20,'units','normalized')
