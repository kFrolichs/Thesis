% Plots the columns in data in a boxplot and adds a scatterplot for the
% individual datapoints
% Assumes data is ordered in the rows

function plotBoxScat(data, varargin)
    %% Settings
    scatMargin = [30, 70];
    dSize      = size(data);
    colData    = [0 .447 .741;.850 .325 .098;.929 .694 .125;.494 ...
                .184 .556;.466 .674 .188;.301 .745 .933;.635 .078 .184];
    pVal       = 0.05 / dSize(2);
    
    if ~isempty(varargin)
        % First check if colors need to be changed
        colPos = returnVarPos(varargin,'colors');
        if colPos ~= 0
            % Assume the next position after colors is the data
            vSize = size(varargin{colPos+1});
            % If it has 3 columns and rows == data columns it is color data
            % If not, ignore and pick standard color values
            if vSize(1) == dSize(2) && vSize(2) == 3
                colData = varargin{colPos+1};
            else
                warning('Colors not recognized. Will use standard colors')
            end
        end

        % Check if we want to Bonferroni Correct
        bonPos = returnVarPos(varargin,'Bonferroni');
        if bonPos ~= 0
            % Bonferroni Correction is automically applied
            if strcmpi(varargin{bonPos+1},'no')
                pVal = 0.05;
            end
        end

        % Check if we want to do one-sided t-test
        tesPos = returnVarPos(varargin,'one_sided');
    end


    %% Plotting
    % Plot the boxplot
    boxplot(data,'colors',colData(1:dSize(2),:)); hold on
    
    % Plot the scatterplot
    xVal = []; col = [];
    for iIdx = 1:dSize(2)
        xVal = [xVal randi(scatMargin,dSize(1),1)/100 ...
            + repelem(iIdx-1,dSize(1),1)]; %#ok<*AGROW>
        col  = [col; repmat(colData(iIdx,:),dSize(1),1)];
    end
    s = scatter(xVal(:),data(:));
    s.CData = col;

    %% Add ttest to see if its larger than zero
%     minVal = min(data,[],'all');
%     
%     for iStat = 1:dSize(2)
%         % t-test
%         if ~exist('tesPos','var') || tesPos == 0
%             [~,p,~,stats] = ttest(data(:,iStat),0);
%         elseif strcmpi(varargin{tesPos+1},'left')
%             [~,p,~,stats] = ttest(data(:,iStat),0,'tail','left');
%         elseif strcmpi(varargin{tesPos+1},'right')
%             [~,p,~,stats] = ttest(data(:,iStat),0,'tail','right');
%         end
%         
%         % Bonferroni corrected p-value
%         if p >= pVal
%             text(iStat*.9,minVal*.9,'\bfn.s.')
%         elseif p <= (pVal/5)
%             text(iStat,minVal*.9,'\bf**')
%         else
%             text(iStat,minVal*.9,'\bf*')
%         end
% 
%         disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p)])
%     end
    
    %% Plot Settings
    xlim([.25, dSize(2)+.25])
end

%% Helper function
function pos = returnVarPos(allVar,varName)
    pos = 0;
    for iVar = 1:numel(allVar)
        if strcmpi(allVar{iVar},varName)
            pos = iVar;
        end
    end
end