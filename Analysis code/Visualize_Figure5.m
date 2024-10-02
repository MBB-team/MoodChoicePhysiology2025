% Setup
    data_directory = [cd filesep 'Data']; %Fill in the directory where the data is stored here
    load('participants')
% Exclusion criterion: choice rate
    choicerates = NaN(size(participants,1),1);
    for ppt = 1:size(participants,1)
        load([data_directory filesep participants.dataset{ppt} filesep 'AllData.mat']);
        choicerates(ppt) = nanmean(AllData.trialinfo.choiceLL);
    end
    i_exclude = find(choicerates<0.05 | choicerates>0.95); %do not visualize modelling results if there is almost no variation in behaviour
    
%% Panel (a) - discount functions
% Prepare
    load([cd filesep 'Results' filesep 'choiceModelBased_noMood'])
    titles = {'delay','risk','effort'};
    costlabels = {'Months to wait','Probability of losing','Effort level'};
% Visualize
    figure
    for type = [2,1,3]
        i_plot = find(type==[2,1,3]);
        subplot(1,3,i_plot); hold on
        data = cell2mat(winning_model_data.discountFunction(:,type));
        N = size(data,2);
        plot(linspace(0,N),data','color',[0.6 0.6 0.6])
        plot(nanmean(data),'color','k','LineWidth',4)
        axis([0,N,0,1]); title(titles{type});
        xlabel(costlabels{type}); xticks(0:25:N); xticklabels(0:3:12); xtickangle(0)
        if type > 1; xticklabels({'0%','25%','50%','75%','100%'}); end
        yticks(0:0.5:1); yticklabels({'0€','15€','30€'}); ylabel('Subjective value'); ylim([0 1])
        box on; axis square
    end
    
%% Panel (b) - psychometric curve of model without mood
    Draw_Psychometric_Curve(winning_model_data,data_directory,false,i_exclude)
    
%% Panel (c) - residuals regressed against mood
% Settings
    n_bins = 9;
    moodcolor = [0,176,80]./255; %green
    emotion_colors = [0.0510    0.4590    0.9920; ... %blue for sad
        0.5, 0.5, 0.5; ... %grey for neutral
        1.0000    0.9    0.1]; %yellow for happy
    color_map = NaN(101,3);
    for i = 1:size(map,2)
        color_map(:,i) = interp1([1;51;101], emotion_colors(:,i), 1:101)';
    end
% Gather regression coefficients and binned residuals
    k_mood = NaN(size(participants,1),2);
    binned_residuals = NaN(size(participants,1),n_bins);
    binned_mood = NaN(size(participants,1),n_bins);
    for ppt = 1:size(participants,1)
        if ~ismember(ppt,i_exclude)
            %data
                disp(['PPT #' num2str(ppt)])
                rated_mood = winning_model_data.ratedMood{ppt};
                residuals = winning_model_data.residuals{ppt};
            %regression
                fit_model = fitglm(rated_mood,residuals);
                k_mood(ppt,:) = fit_model.Coefficients.Estimate';
            %binning
                [sorted_mood,i_sorted] = sort(rated_mood); %sorted mood ratings
                sorted_res = residuals(i_sorted);
                n_small = floor(length(sorted_mood)/n_bins); %number of trials in small bins
                largebins = rem(length(sorted_mood),n_small); %number of large bins (n_small+1 trials)
                smallbins = n_bins - largebins; %number of small bins (n_small trials)
                %SEQUENCE OF SORTED TRIALS: note, smaller bins on the edges, asymmetric with more small bins on the lower end of DV
                    bins_n = [n_small*ones(1,ceil(smallbins/2)),(n_small+1)*ones(1,largebins),n_small*ones(1,floor(smallbins/2))];
                bins_cum = [0 cumsum(bins_n)];
                for bin = 1:n_bins
                    i_bin = bins_cum(bin)+1:bins_cum(bin+1); %trial index for given bin
                    binned_mood(ppt,bin) = mean(sorted_mood(i_bin));
                    binned_residuals(ppt,bin) = mean(sorted_res(i_bin));
                end
        end
    end
% Visualize
    figure; hold on; box on
    %Mean regression of residuals against mood ratings
        X = linspace(-2,2);
        Y = k_mood(:,1) + k_mood(:,2).*X;
        E = nanstd(Y) ./ sqrt(sum(~isnan(Y)));
        Y = nanmean(Y);
        patch([X X(end:-1:1)], [Y+E Y(end:-1:1)-E(end:-1:1)], moodcolor, 'facealpha', 0.25, 'Edgecolor', 'none');
        plot(X, Y, 'linestyle', '-', 'LineWidth', 2, 'color', moodcolor);
    %Binned residuals overlaid
        x = nanmean(binned_mood);
        y = nanmean(binned_residuals);
        err = nanstd(binned_residuals)./sqrt(sum(~isnan(binned_residuals)));
        E = errorbar(x,y,err,'k','CapSize',4,'LineStyle','none');
        E.LineWidth = 1;
        scatter(x,y,60,x,'filled','MarkerEdgeColor','k')
    %Layout
        xlabel('Mood score'); ylabel('Choice residuals')
        axis([-2.5,2.5,-0.04,0.04])
        colormap(color_map); colorbar; caxis([-1.6, 1.6])

%% Panel (d) - modelled psychometric curves with mood parameter

% Load model
    load([cd filesep 'Results' filesep 'choiceModelBased_withMood'])
    Draw_Psychometric_Curve(winning_model_data,data_directory,true,i_exclude)
    
%% Panel (e) - mood bias parameter distribution

% Note: this figure requires the RainCloudPlots toolbox and its dependencies to be installed:
% https://github.com/RainCloudPlots/RainCloudPlots/tree/master
    addpath(genpath('C:\Users\Roeland\Documents\MATLAB\toolbox\colorspace'))
    addpath(genpath('C:\Users\Roeland\Documents\MATLAB\toolbox\Colorbrewer'))
    addpath(genpath('C:\Users\Roeland\Documents\MATLAB\toolbox\Robust_Statistical_Toolbox'))
    addpath(genpath('C:\Users\Roeland\Documents\MATLAB\toolbox\RainCloudPlots'))

% Settings
    moodcolor = [0.2000    0.6000    0.2000];
    i_exclude = 52; %no variability in behaviour
    betaMood = [Inversionresults.muPhi.betaMood]';
    betaMood(i_exclude) = NaN;
    
% Plot
    figure('Position',[843 190 833 619])
    subplot(2,5,1:2); hold on
        X = 1:2:7;
        studies = participants.study>2;
        for s = 1:2
            betas  = betaMood(studies==(s-1));
            swarmchart(X(s)*ones(size(betas)),betas,25,moodcolor(end,:),'filled','MarkerFaceAlpha',0.66,'MarkerEdgeAlpha',0.66)
            line(X(s)+[-0.33,0.33],nanmean(betas)*ones(1,2),'color','k','LineWidth',1.5)
        end        
        ylim([-0.75 1.25]); yticks(-0.75:0.25:1.25); box on; title('')
        xlim([-1 7]); xticks(X); xticklabels([{'exp'},{'conf'}]); xtickangle(0); 
    subplot(2,5,5)
        h1 = raincloud_plot(betaMood, 'line_width', 1.5,'color',moodcolor);
        camroll(90); set(gca, 'XDir','reverse')
        ylim([-2 5]); xlim([-0.75 1.25]); xticks([]); yticks([]); axis off
        
%% Auxiliary functions
function Draw_Psychometric_Curve(winning_model_data,data_directory,mood_in_model,i_exclude)

%Settings
    n_bins = 9; 
    DV_range = [-5,5];
    X_plot_range = [-1,1];
    select_conditions = [1,2,5]; %Happy, sad, neutral
    emotioncolors = [1,0.756078431372549,0.107254901960784; %Yellow for happy
        0.0509803921568627,0.458823529411765,0.992156862745098]; %Blue for sad
    load('participants') %#ok<LOAD>

%Sort data
    P_LL_all = cell(size(participants,1),3);
    Y_all = cell(size(participants,1),3);
    DV_all = cell(size(participants,1),3);
    for ppt = 1:size(participants,1)
        %Load the data
            disp(['PPT #' num2str(ppt)])
            load([data_directory filesep participants.dataset{ppt} filesep 'AllData.mat']) %#ok<LOAD>
            include_choicetypes = true(1,3);   
            for type = 1:3
                if participants.study(ppt) == 1 && type == 2 %Excude probability discounting from study 1 
                    include_choicetypes(type) = 0;
                elseif participants.study(ppt) == 2 && type == 1 %Exclude delay discounting from study 2
                    include_choicetypes(type) = 0;
                end
            end
            select = ismember(AllData.trialinfo.choicetype,find(include_choicetypes)) & ismember(AllData.trialinfo.condition,[1,2,5]);     
            conditions = AllData.trialinfo.condition(select);
        %Sort per condition
            for cond = 1:length(select_conditions)
                P_LL = winning_model_data.P_LL{ppt};
                P_LL_all{ppt,cond} = P_LL(conditions==select_conditions(cond));
                Y = winning_model_data.Y{ppt};
                Y_all{ppt,cond} = Y(conditions==select_conditions(cond));
                DV = winning_model_data.DV{ppt};
                DV_all{ppt,cond} = DV(conditions==select_conditions(cond));
            end
    end %for ppt
           
%Draw
    figure; hold on; box on
    %Draw reference lines
        line([0 0],[0 1],'Color',[0.75 0.75 0.75],'LineStyle','--','LineWidth',0.3)
        line([min(X_plot_range) max(X_plot_range)],[0.5 0.5],'Color',[0.75 0.75 0.75],'LineStyle','--','LineWidth',0.3)
    %Psychometric curve per condition
        for cond = [2,1] %plot happy and sad only                
            %Bin the data
                binned_DV = NaN(size(DV_all,1),n_bins); 
                binned_P_LL = binned_DV; 
                binned_Y = binned_DV;
                for ppt = 1:size(DV_all,1)
                    %Get participant data
                        P_LL = P_LL_all{ppt,cond}; %Modelled choice probabilities
                        Y = Y_all{ppt,cond}; %Observed decisions
                        DV = DV_all{ppt,cond}; %Modelled value difference            
                    %Set out-of-bounds DV values to NaN
                        DV(DV<DV_range(1) | DV>DV_range(2)) = NaN; 
                    %Sort according to DV
                        [DV,I] = sort(DV); %Sort DV
                        P_LL = P_LL(I); %Sort P_SS according to the order of DV    
                        Y = Y(I); %Sort choices according to the order of DV    
                    %Get rid of NaNs
                        P_LL = P_LL(~isnan(DV));
                        Y = Y(~isnan(DV));
                        DV = DV(~isnan(DV));
                    %Store in bins
                        if ~isempty(DV) && ~ismember(ppt,i_exclude)
                            n_small = floor(length(DV)/n_bins); %number of trials in small bins
                            largebins = rem(length(DV),n_small); %number of large bins (n_small+1 trials)
                            smallbins = n_bins - largebins; %number of small bins (n_small trials)
                            %SEQUENCE OF SORTED TRIALS: note, smaller bins on the edges, asymmetric with more small bins on the lower end of DV
                                bins_n = [n_small*ones(1,ceil(smallbins/2)),(n_small+1)*ones(1,largebins),n_small*ones(1,floor(smallbins/2))];
                            bins_cum = [0 cumsum(bins_n)];
                            for bin = 1:n_bins
                                i_bin = bins_cum(bin)+1:bins_cum(bin+1); %trial index for given bin
                                binned_DV(ppt,bin) = mean(DV(i_bin));
                                binned_P_LL(ppt,bin) = mean(P_LL(i_bin));
                                binned_Y(ppt,bin) = mean(Y(i_bin));
                            end
                        else
                            binned_DV(ppt,:) = NaN(1,n_bins);
                            binned_P_LL(ppt,:) = NaN(1,n_bins);
                            binned_Y(ppt,:) = NaN(1,n_bins);
                        end
                end %for ppt
            %Visualize
                color = emotioncolors(cond,:);
                if ~mood_in_model %display observed behaviour
                    SEM = nanstd(binned_Y,[],1)./sqrt(sum(~isnan(binned_Y)));
                    errorbar(nanmean(binned_DV),nanmean(binned_Y,1),SEM,'o','Color',color,...
                        'MarkerSize',5,'MarkerEdgeColor',color,'MarkerFaceColor',color,'LineWidth',1.2);
                    plot(nanmean(binned_DV),nanmean(binned_Y),'color',color,'LineWidth',1.5)
                else %display model predictions
                    X = nanmean(binned_DV);
                    E = nanstd(binned_P_LL)./sqrt(sum(~isnan(binned_P_LL)));
                    Y = nanmean(binned_P_LL);
                    patch([X X(end:-1:1)], [Y+E Y(end:-1:1)-E(end:-1:1)], color, 'facealpha', 0.33, 'Edgecolor', 'none');
                    plot(X, Y, 'linestyle', '-', 'LineWidth', 1.25, 'color', color);
                end
            %Layout
                xlabel('Decision value'); ylabel('P(costly option)')
                xlim(X_plot_range); 
                ylim([0 1])
        end
  
return
        
%Shift bottom row subplots up
    ha = gca;
    pos = ha.Position;
    pos(2) = pos(2)+0.05;
    ha.Position = pos;

%Figure size
    set(gcf,'Position',[680 364 673 614]);
    
end
    