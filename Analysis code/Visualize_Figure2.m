% Generate figure 2: validation of emotion indcution by subjective ratings.

% Load data: this variable is obtained by running 'Analysis_ModelFree.m'
    load([cd filesep 'Results' filesep 'InductionRatings'],'allRatings')
    load('participants')
    
% Settings
    distribution = 'violin'; %barplot with scattered data points if not "violin"
    which_ratings = {'RateHappy','RateSad','RateAngry','RateFear','Mood'};
    EmotionColors =     [1.0000    0.7961    0.1373     %yellow: happy
                        0.0510    0.4588    0.9922      %blue: sad
                        0.8706    0.2196    0.0784      %red: angry
                        0.6196    0.0392    0.5216      %purple: fear
                        0.5000    0.5000    0.5000      %grey: neutral
                        0.2000    0.6000    0.2000];    %green: mood
                    
% Produce figure 2
    %ratings (fig 2 a, b)
        figure
        emo_ind1 = [1,2,5,3,4];
        emo_ind2 = [1,2,5];
        dx = 0.15; %#ok<*NASGU> %horizontal offset of the errorbarv
        titles = {'Happiness induction','Sadness induction','Neutral induction','Anger induction','Fear induction'};
        studies = {'exploratory','confirmatory'};
        for ii = 1:length(emo_ind1)
            emo = emo_ind1(ii);
            for R = 1:4 % H, S, A, F
                if ismember(emo,emo_ind2); N=2; else; N = 1; end
                data = cell2mat(allRatings.perStudy.(which_ratings{R})(:,emo));
                for j = 1:N
                    ha = subplot(2,5,ii+(j-1)*length(emo_ind1)); hold on; %box on
                    data2 = allRatings.(which_ratings{R})(strcmp(participants.experiment,studies{j}),emo);
                    if strcmp(distribution,'violin')
                        violinplot(R*ones(size(data2)),data2,'FaceColor',EmotionColors(R,:))
                    else
                        M = data(1,j); SEM = data(2,j);
                        bar(R,M,'EdgeColor','none','LineWidth',1.5,'FaceColor',EmotionColors(R,:)); %#ok<*UNRCH>
                        errorbar(R+dx,M,SEM,'vertical','k','CapSize',0,'LineWidth',1.25)
                        scatter(R*ones(size(data2))-dx,data2,5,(EmotionColors(R,:)).^2.5,'filled')
                    end
                    ylim([0,1]); yticks(0:0.2:1); ylabel('Rating score (raw)')
                    xticks(1:4); xticklabels({'H','S','A','F'}); xlabel('Rating'); xlim([0.33,4.66])
                    if emo~=1 
                        ha.YAxis.Visible = 'off'; % remove y-axis
                    end
                end
            end 
        end
        for ii = 1:length(titles); subplot(2,5,ii); title(titles{ii}); end
            
    %mood (fig 2 c)
        subplot(2,20,34:40); hold on
        ngroups = 2;
        catwidth = min(0.8, ngroups/(ngroups + 1.5));
        if ~strcmp(distribution,'violin')
            dx = 0.05; %horizontal offset of the errorbar
            M = cell2mat(allRatings.perStudy.Mood(1,[1,2,5])');
            SEM = cell2mat(allRatings.perStudy.Mood(2,[1,2,5])');
            B = bar(M,'EdgeColor','none','FaceColor',EmotionColors(6,:));
        end
        for group = 1:2
            for cat = 1:3
                data2 = allRatings.Mood(strcmp(participants.experiment,studies{group}),emo_ind2(cat));
                x = cat - catwidth/2 + (2*group-1) * catwidth / (2*ngroups);
                if strcmp(distribution,'violin')
                    violinplot(ones(size(data2))*x,data2,'FaceColor',EmotionColors(end,:))
                else
                    y = M(cat,group);
                    e = SEM(cat,group);
                    errorbar(x+dx,y,e,'vertical','k','CapSize',0,'LineWidth',1.25)
                    scatter(x*ones(size(data2))-dx,data2,8,(EmotionColors(end,:)).^2.5,'filled')
                end
            end
        end        
        ylim([-1.5,1.5]); xlim([0.33 3.66])
        yticks(-1.5:0.5:1.5); %ylabel('Mood score [z]'); title('All inductions')
        xticks(1:3); xticklabels({'Happiness','Sadness','Neutral'}); xlabel('Induction')
        ylabel('Mood score (Z)')