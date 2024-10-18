% Visualize figure 6: Prediction of choice behaviour from physiological mood proxy

% Load results
    load('physiology_analysis.mat')

%% Panel a
    labels = {'pupil dilation','skin conductance','zygomaticus','corrugator','mood score','choice rate'};
    % Average correlation coefficients (set diagonals to NaN)
        select_correlations = physiology_correlations(1:6,1:6,:);
        meanCorrelations = nanmean(select_correlations,3);        
        for i = 1:length(meanCorrelations)
            meanCorrelations(i,i) = 0;
        end
    %P-values of correlation coefficients
        phys_P = NaN(size(meanCorrelations)); 
        n = 14; %Bonferroni
        for i = 1:numel(phys_P)
            [I,J] = ind2sub(size(phys_P),i);
            [~,p] = ttest(select_correlations(I,J,:));
            phys_P(i) = p * n; %bonerroni
        end
        phys_P(phys_P>1) = 1;
    %Figure
        hf3 = figure; hold on; 
        hf3.Position = [945 522 576 322];
        i_grey = find(phys_P>0.05);
        imagesc(meanCorrelations)
        set(gca,'YDir','reverse'); c = colorbar; caxis([-0.25,0.25]) %#ok<CAXIS>
        xticks(1:length(meanCorrelations)); yticks(1:length(meanCorrelations)); xticklabels(labels); yticklabels(labels); xtickangle(30)    
        for i = 1:numel(meanCorrelations)
            [I,J] = ind2sub(size(meanCorrelations),i);
            %Draw grey square
                if ismember(i,i_grey)
                    rectangle('Position',[I-0.5,J-0.5,1 1],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
                end
            %Write txt
                if round(phys_P(i),3) == 0
                    txt_p = '(<0.001)';
                else
                    txt_p = ['(' num2str(round(phys_P(i),3)) ')'];
                end
                if I ~= J && J <= 4
                    boxtext = {num2str(round(meanCorrelations(i),3)); txt_p};
                    if ismember(i,i_grey)
                        text(I,J,boxtext,'HorizontalAlignment','center','FontSize',9,'color',[0.5,0.5,0.5])
                    else
                        text(I,J,boxtext,'HorizontalAlignment','center','FontSize',9)
                    end
                elseif I == J
                     rectangle('Position',[I-0.5,J-0.5,1 1],'FaceColor',[1,1,1],'EdgeColor','none')
                end
        end
        set(gca,'xaxisLocation','top')
        box off
        ylim([0.5,4.5]); %axis square; axis tight
        xlim([0.5,6.5])
        rectangle('Position',[0.5,0.5,1,4],'EdgeColor',[1,1,1],'FaceColor',[1,1,1])
        rectangle('Position',[1.5,1.5,1,3],'EdgeColor',[1,1,1],'FaceColor',[1,1,1])
        rectangle('Position',[2.5,2.5,1,2],'EdgeColor',[1,1,1],'FaceColor',[1,1,1])
        rectangle('Position',[3.5,3.5,1,1],'EdgeColor',[1,1,1],'FaceColor',[1,1,1])
        plot([4.5 4.5],[0.5 length(meanCorrelations)+0.5],'color',[0 0 0],'LineWidth',3)
        colormap cool

%% Panel b
    hf1 = figure; hold on
    hf1.Position = [916 520 290 395];% [680 538 402 440];
    %Data
        X = nanmean(binned_physiology.rated_mood);
        Y1 = nanmean(binned_physiology.proxy_arousal);
        E1 = RH_SEM(Y1);
        Y2 = nanmean(binned_physiology.proxy_valence);
        E2 = RH_SEM(Y2);    
    %Top: arousal
        ha1 = subplot(2,1,1); hold on
        ar_color = [224,130,20]./255;
        plot(X,Y1,'LineWidth',2.5,'color',[253 184 99]./255)
        errorbar(X,Y1,E1,'LineStyle','none','Marker','o','MarkerSize',4,'Color',ar_color,... ,'MarkerEdgeColor','ar_color',...
            'LineWidth',1.5,'CapSize',4,'MarkerFaceColor',ar_color)
        xlim([-2,2]); ha1.XAxis.Visible = 'off';
        xlabel(''); ylabel('')
    %Bottom: valence
        subplot(2,1,2); hold on
        val_color = [53,151,143]./255;
        plot(X,Y2,'LineWidth',2.5,'color',[199,234,229]./255)
        errorbar(X,Y2,E2,'LineStyle','none','Marker','o','MarkerSize',4,'Color',val_color,... ,'MarkerEdgeColor','ar_color',...
            'LineWidth',1.5,'CapSize',4,'MarkerFaceColor',val_color)
        xlim([-2,2]); ylim([-0.35,0.45]); yticks(-0.3:0.15:0.45)
        xlabel(''); ylabel('')

%% Panel c
% Make figure
    hf2 = figure; hold on; box on
    hf2.Position = [675 235 478 245];
% Color map
    moodProxyColor = [0,176,80]./255; %green
    emotion_colors = [0.0510    0.4590    0.9920; ... %blue for sad
        0.5, 0.5, 0.5; ... %grey for neutral
        1.0000    0.9    0.1]; %yellow for happy
    color_map = NaN(101,3);
    for i = 1:size(color_map,2)
        color_map(:,i) = interp1([1;51;101], emotion_colors(:,i), 1:101)';
    end
% Regressions
    betas_moodProxy = cell2mat(beta_mood_residuals(:,2));
    X = linspace(-0.3,0.3);
    Y = betas_moodProxy(:,1) + betas_moodProxy(:,2).*X;
    E = RH_SEM(Y);
    Y = nanmean(Y);
    patch([X X(end:-1:1)], [Y+E Y(end:-1:1)-E(end:-1:1)], moodProxyColor, 'facealpha', 0.33, 'Edgecolor', 'none');
    plot(X, Y, 'linestyle', '-', 'LineWidth', 1.25, 'color', moodProxyColor);
%  Binned residuals
    errorbar(nanmean(binned_physiology.proxy_mood),nanmean(binned_physiology.residuals),RH_SEM(binned_physiology.residuals),...
        'LineStyle','none','Color','k','LineWidth',1.5,'CapSize',5)
    scatter(nanmean(binned_physiology.proxy_mood),nanmean(binned_physiology.residuals),60,nanmean(binned_physiology.rated_mood),...
        'filled','MarkerEdgeColor','k')
    C = colorbar; colormap(color_map); 
    xlim([-0.35,0.35])

%% Panel d
% Figure
    hf4 = figure; hold on
    hf4.Position = [916 236 290 227];
% Visualize correlation
    betas_ratedMood = cell2mat(beta_mood_residuals(:,1));
    fit_mdl = fitglm(betas_ratedMood(:,2),betas_moodProxy(:,2));
    X = linspace(-0.075,0.175);
    [Y,Y_CI] = predict(fit_mdl,X','Alpha',0.05); 
    E = Y_CI(:,2)'- Y';
    Y = Y';
% Plot individual values
    color = [0,90,50]./255;    
    patch([X X(end:-1:1)], [Y+E Y(end:-1:1)-E(end:-1:1)], color, 'facealpha', 0.05, 'Edgecolor', color);
    plot(X, Y, 'linestyle', '-', 'LineWidth', 1, 'color', (color));
    scatter(betas_ratedMood(:,2),betas_moodProxy(:,2),20,color,'filled')
    xlim([-0.08 0.18]); box on
    xlabel(''); ylabel('')

%% Auxiliary function

function sem = RH_SEM(data)
    sem = nanstd(data) ./ sqrt(sum(~isnan(data))); %#ok<NANSTD>
end