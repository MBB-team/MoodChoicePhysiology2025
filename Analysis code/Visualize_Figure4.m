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
    %Summarize results
        choicerates = choice_mdlfree.choiceRate(:,1:4)-choice_mdlfree.choiceRate(:,5); %difference w.r.t. neutral
        results = struct;
        results.M_expl = nanmean(choicerates(participants.study<3,:));
        results.SEM_expl = RH_SEM(choicerates(participants.study<3,:));
        results.M_conf = nanmean(choicerates(participants.study>=3,:));
        results.SEM_conf = RH_SEM(choicerates(participants.study>=3,:));
        for study = 1:4
            results.M_perstudy(study,:) = nanmean(choicerates(participants.study==study,:));
            results.SEM_perstudy(study,:) = RH_SEM(choicerates(participants.study==study,:));
        end
    %Visualize
        figure(1)
        for i = 1:2
            subplot(2,3,i); cla; hold on
            if i == 1
                M = results.M_expl; SEM = results.SEM_expl; which_studies = 1:2; which_emotions = 1:4; which_title = 'expl';
            else
                M = results.M_conf; SEM = results.SEM_conf; which_studies = 3:4; which_emotions = 1:2; which_title = 'Confirmatory';
            end
            for ii = which_emotions
                B = bar(ii,M(ii),'EdgeColor','none','FaceColor',EmotionColors(ii,:));
                errorbar(ii, M(ii), SEM(ii), 'color', 'k','CapSize',8,'LineWidth',1)
            end
            ha = gca; ha.FontSize=10;
            if i==1
                ylim([-0.05,0.05]); yticks(-0.05:0.025:0.05); ylabel('\Delta costly choice rate','FontSize',11); xlabel('Induction'); yticklabels({'-5%','-2.5%','0%','2.5%','5%'});
            else
                ylim([-0.025,0.025]); yticks(-0.025:0.0125:0.025);  xlabel('Induction'); yticklabels({'-2.5%','-1.25%','0%','1.25%','2.5%'});
            end
            xticks(1:4); xticklabels({'H','S','A','F'}); xlim([0.33,4.66]); title(''); box on
        end
        
%% Panel (c) Regressions of choice against mood 
    figure(1); subplot(2,3,4:5); cla; hold on
    plot_betas = cell(1,2);
    studies = 1+double(participants.study>2); %1 exploratory / 2 confirmatory
    s_names = {'expl','conf'}; t_names = {'delay','risk','effort'};
    T = table;
    X = [1,6,11; 3,8,13];
    type_order = [2,1,3]; %risk/delay/effort
    for S = 1:2
        for type=1:3
            betas = choice_mdlfree.beta_mood(studies==S,type_order(type));
            swarmchart(X(S,type)*ones(size(betas)),betas,25,EmotionColors(end,:),'filled','MarkerFaceAlpha',0.66,'MarkerEdgeAlpha',0.66)
            %line(X(S,type)+[-0.33,0.33],nanmean(betas)*ones(1,2),'color','k','LineWidth',1.5)
        end
    end
    xticks(unique(X)'), xticklabels(repmat({'exp','conf'},1,3))
    ylabel('Weight on mood','FontSize',11); ylim([-0.6,1.3])
    box on
    
%% Panel (d) Response Time
    M = [nanmean([choice_mdlfree.RT_perEmo_SS(:,1)-choice_mdlfree.RT_perEmo_SS(:,3), choice_mdlfree.RT_perEmo_SS(:,2)-choice_mdlfree.RT_perEmo_SS(:,3)]);...
        nanmean([choice_mdlfree.RT_perEmo_LL(:,1)-choice_mdlfree.RT_perEmo_LL(:,3), choice_mdlfree.RT_perEmo_LL(:,2)-choice_mdlfree.RT_perEmo_LL(:,3)])];
    E = [RH_SEM([choice_mdlfree.RT_perEmo_SS(:,1)-choice_mdlfree.RT_perEmo_SS(:,3), choice_mdlfree.RT_perEmo_SS(:,2)-choice_mdlfree.RT_perEmo_LL(:,3)]);...
        RH_SEM([choice_mdlfree.RT_perEmo_LL(:,1)-choice_mdlfree.RT_perEmo_LL(:,3), choice_mdlfree.RT_perEmo_LL(:,2)-choice_mdlfree.RT_perEmo_LL(:,3)])];
    figure(1); subplot(2,3,6); hold on
        b = bar(M,'EdgeColor','none');
        b(1).FaceColor = EmotionColors(1,:);
        b(2).FaceColor = EmotionColors(2,:);
        errorbar([0.8571,1.14286,1.8571,2.1486],[M(1,1),M(1,2),M(2,1),M(2,2)],[E(1,1),E(1,2),E(2,1),E(2,2)],'color','k','CapSize',5,'LineWidth',1,'LineStyle','none')
        xticks([1,2]); xticklabels({'uncostly','costly'}); 
        ha=gca;ha.FontSize=10; 
        ylim([-0.075,0.15]); ylabel('\Delta RT (Z)','FontSize',11)
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