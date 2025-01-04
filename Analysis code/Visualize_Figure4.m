% Reproduce figure 4: model-free choice data

% Settings
close all
load('participants')
load([cd filesep 'Results' filesep 'choiceModelFree'])
EmotionColors = [...
        1.0000    0.7961    0.1373 %H: yellow
        0.0510    0.4588    0.9922 %S: blue
        0.8706    0.2196    0.0784 %A: red
        0.6196    0.0392    0.5216 %F: purple
        0.5000    0.5000    0.5000 %N: grey
        0.2000    0.6000    0.2000]; %all: green
    
%% Panel (a) Choice rates for exploratory vs. confirmatory studies
    choicerates = choice_mdlfree.choiceRate(:,1:4)-choice_mdlfree.choiceRate(:,5); %difference w.r.t. neutral
    figure(1)
    set(gcf,'Position',[520 300.5000 581 497.5000])
    for i = 1:2 %exploratory and confirmatory
        ha = subplot(2,3,i); hold on
        line([0,4],[0,0],'color',0.8*ones(1,3),'LineWidth',0.5)
        if i == 1; ii = strcmp(participants.experiment,'exploratory');
        else; ii = strcmp(participants.experiment,'confirmatory'); 
            ha.YAxis.Visible = 'off'; % remove y-axis
        end
        data = {choicerates(ii,1),choicerates(ii,2),choicerates(ii,3),choicerates(ii,4)};
        RH_Boxchart(data,EmotionColors(1:4,:))
        ha = gca; ha.FontSize=10;
        ylim([-0.35,0.35]); yticks(-0.3:0.15:0.3); ylabel('\Delta costly choice rate','FontSize',11); yticklabels({'-30%','-15%','0%','15%','30%'});
        xticklabels({'H','S','A','F'}); xlabel('Induction'); xlim([0.5,1.5])
    end
        
%% Panel (c) Regressions of choice against mood 
    figure(1); 
    for i = 1:2 %exploratory and confirmatory
        ha = subplot(2,3,3+i); hold on
        line([0,1.5],[0,0],'color',0.8*ones(1,3),'LineWidth',0.5)
        if i == 1; ii = strcmp(participants.experiment,'exploratory');
        else; ii = strcmp(participants.experiment,'confirmatory'); 
            ha.YAxis.Visible = 'off'; % remove y-axis
        end
        data = {choice_mdlfree.beta_mood(ii,2),choice_mdlfree.beta_mood(ii,1),choice_mdlfree.beta_mood(ii,3),choice_mdlfree.beta_mood(ii,4)};
        RH_Boxchart(data,EmotionColors(6,:))
        ha = gca; ha.FontSize=10;
        ylim([-1,1]); yticks(-1:0.5:1); ylabel('Weight on mood','FontSize',11);
        xticklabels({'risk','delay','effort','all'}); xlim([0.5,1.5])
    end
    
%% Panel (d) Response Time
    % set(gcf,'Position',[520 182.5000 791 615.5000])
    data = {NaN(size(participants,1),2),...
        [choice_mdlfree.RT_perEmo_SS(:,1)-choice_mdlfree.RT_perEmo_SS(:,3),...
        choice_mdlfree.RT_perEmo_LL(:,1)-choice_mdlfree.RT_perEmo_LL(:,3)],...
       [choice_mdlfree.RT_perEmo_SS(:,2)-choice_mdlfree.RT_perEmo_SS(:,3)...
        choice_mdlfree.RT_perEmo_LL(:,2)-choice_mdlfree.RT_perEmo_LL(:,3)],...
        NaN(size(participants,1),2)}; %for visual arrangement
    figure(1); subplot(2,3,6); hold on
        line([0,3],[0,0],'color',0.8*ones(1,3),'LineWidth',0.5)
        RH_Boxchart(data,EmotionColors([1,1,2,2],:));
        xticks([1,2]); xticklabels({'uncostly','costly'}); 
        ha=gca;ha.FontSize=10; 
        ylim([-0.6,0.8]); yticks(-0.6:0.2:0.8); xlim([0.6,2.4]); ylabel('\Delta RT (Z)','FontSize',11)
        box on
        
%% Panel (b) Choice rates across studies
% Note, this figure requires the RainCloudPlot toolbox and all of its dependencies to be added to the MATLAB path:
% https://github.com/RainCloudPlots/RainCloudPlots

%Get data
    data = choicerates(:,[1,2]);
    deltadata = data(:,1)-data(:,2);
%Colors
    color_H = EmotionColors(1,:);
    color_S = EmotionColors(2,:);    
%Visualize
    hf2 = figure(2);
    hf2.Position = [2316 225 729 520];
    subplot(2,3,1); hold on
        h1 = raincloud_plot(data(:,1), 'color', color_H, 'line_width', 1.5,'density_type', 'ks');
        xlim([-0.35,0.35]); xticks(-0.4:0.1:0.4); 
        ylim([-7,50])
        camroll(90); axis off
    subplot(2,3,2); hold on
        hf3 = raincloud_plot(data(:,2), 'color', color_S, 'line_width', 1.5,'density_type', 'ks');
        ylim([-7,50])
        xlim([-0.35,0.35]); xticks(-0.4:0.1:0.4); 
        camroll(90); set(gca, 'XDir', 'reverse'); axis off
    subplot(2,3,3); cla; hold on
        [sorted,I] = sort(-deltadata,'Descend');
        norm_sorted = RH_Normalize(sorted);
        CR = data(I,:);
        %Lines
            pos_xdata = RH_RandInt(size(CR),0.3,0.4);
            neg_xdata = RH_RandInt(size(CR),0.6, 0.7);
            color_low = 0.85*ones(1,3); %light gray
            for ppt = 1:length(I) %sorted from positive to negative (dominant effect last)
                delta = CR(ppt,2)-CR(ppt,1); 
                norm_delta = norm_sorted(find(sorted==delta,1,'first'));
                color = ((norm_delta - min(norm_sorted))./(max(norm_sorted)-min(norm_sorted)))*color_low;
                Y = CR(ppt,:);
                X = [pos_xdata(I(ppt)), neg_xdata(I(ppt))];
                plot(X,Y,'Color',color,'LineWidth',1)
            end
        %Scatter dots            
            Y = CR;
            X = [pos_xdata(I), neg_xdata(I)];
            scatter(X(:,1),Y(:,1),10,'MarkerFaceColor',color_H,'MarkerEdgeColor','none','LineWidth',0.5);
            scatter(X(:,2),Y(:,2),10,'MarkerFaceColor',color_S,'MarkerEdgeColor','none','LineWidth',0.5);
        %Layout
            xlim([-0.1 1.1]); xticks([]); box on
            ylim([-0.35,0.35]); yticks(-0.3:0.1:0.3); yticklabels({'-30%','-20%','-10%','0%','+10%','+20%','+30%'});
            ha=gca;ha.FontSize=10;
            ylabel('\Delta costly choice rate','FontSize',11)
    
%% Auxiliary functions
function sem = RH_SEM(data)
    % Calculates standard error of the mean
    sem = nanstd(data) ./ sqrt(sum(~isnan(data)));
end

function [output] = RH_RandInt(sizevector,lowerlimit,upperlimit)
    % Generates a random number matrix of values between lowerlimit and upperlimit.
    output = lowerlimit + (upperlimit-lowerlimit) .* rand(sizevector);
end