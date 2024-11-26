% Generate figure 2: validation of emotion indcution by subjective ratings.

% Load data: this variable is obtained by running 'Analysis_ModelFree.m'
    load([cd filesep 'Results' filesep 'InductionRatings'],'allRatings')
    load('participants')
    
% Settings
    distribution = 'boxchart'; %barplot with scattered data points if not "violin" or "boxchart"
    which_ratings = {'RateHappy','RateSad','RateAngry','RateFear','Mood'};
    EmotionColors =     [1.0000    0.7961    0.1373     %yellow: happy
                        0.0510    0.4588    0.9922      %blue: sad
                        0.8706    0.2196    0.0784      %red: angry
                        0.6196    0.0392    0.5216      %purple: fear
                        0.5000    0.5000    0.5000      %grey: neutral
                        0.2000    0.6000    0.2000];    %green: mood
                    
% Produce figure 2
    %ratings (fig 2 a, b)
        figure('Position',[253.5000 336 928 462])
        emo_ind1 = [1,2,5,3,4];
        emo_ind2 = [1,2,5];
        dx = 0.15; %#ok<*NASGU> %horizontal offset of the errorbarv
        titles = {'Happiness induction','Sadness induction','Neutral induction','Anger induction','Fear induction'};
        studies = {'exploratory','confirmatory'};
        rating_data = {allRatings.RateHappy,allRatings.RateSad,allRatings.RateAngry,allRatings.RateFear};
        for ii = 1:length(emo_ind1)
            emo = emo_ind1(ii);
            data = cellfun(@(x) x(:,emo),rating_data,'UniformOutput',false);
            if ismember(emo,emo_ind2); N=2; else; N = 1; end %number of studies
            for j = 1:N
                ha = subplot(2,5,ii+(j-1)*length(emo_ind1)); hold on; %box on
                data2 = cellfun(@(x) x(strcmp(participants.experiment,studies{j}),:),data,'UniformOutput',false);
                RH_Boxchart(data2,EmotionColors(1:4,:));
                ylim([0,1]); yticks(0:0.2:1); xticklabels({'H','S','A','F'}); xlabel('Rating'); xlim([0.5,1.5])
                if j == 1; ylabel({'Rating score (raw)';'Exploratory'}); else; ylabel({'Rating score (raw)';'Confirmatory'}); end
                if emo~=1 
                    ha.YAxis.Visible = 'off'; % remove y-axis
                end
            end
        end
        for ii = 1:length(titles); subplot(2,5,ii); title(titles{ii}); end
        
    %mood (fig 2 c)
        subplot(2,20,34:40); hold on
        line([0,4],[0,0],'color',0.8*ones(1,3),'LineWidth',0.5)
        data = {allRatings.Mood(strcmp(participants.experiment,studies{1}),emo_ind2),...
            allRatings.Mood(strcmp(participants.experiment,studies{2}),emo_ind2)};
        RH_Boxchart(data,EmotionColors(6,:))
        ylim([-1.5,1.5]); xlim([0.33 3.66])
        yticks(-1.5:0.5:1.5); %ylabel('Mood score [z]'); title('All inductions')
        xticks(1:3); xticklabels({'Happiness','Sadness','Neutral'}); xlabel('Induction')
        ylabel('Mood score (Z)')