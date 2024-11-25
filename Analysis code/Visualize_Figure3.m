% Visualises time-resolved physiological responses to emotion inductions.
% Also applies random field theory for cluster-level significance testing.
% For this, the VBA toolbox and Fieldtrip must be installed: 
% https://github.com/MBB-team/VBA-toolbox
% https://www.fieldtriptoolbox.org/

% Prepare
    load('physiology_averages.mat')
    EmotionColors = [...
        1.0000    0.7961    0.1373
        0.0510    0.4588    0.9922
        0.8706    0.2196    0.0784
        0.6196    0.0392    0.5216
        0.5000    0.5000    0.5000
        0.2000    0.6000    0.2000];
    c_HN = [1;0;-1]; %Happy vs. neutral contrast for RFT
    c_SN = [0;1;-1]; %Sad vs. neutral contrast for RFT
    c_AN = [1;0;-1]; %Anger vs. neutral contrast for RFT
    c_FN = [0;1;-1]; %Fear vs. neutral contrast for RFT
    c_HS = [1;-1;0]; %Happy vs. sad contrast for RFT
    figure

%% Pupil
    %Happy/Sad/Neutral
        subplot(4,2,1); hold on; box on
        %Get data
            pupil = physiology_averages_HSN.pupil(:,[1,2,5]);
            pupil = pupil(~cellfun(@isempty,pupil(:,1)),:);
            X = kron(eye(3),ones(size(pupil,1),1));
            pupil = {cell2mat(pupil(:,1)),cell2mat(pupil(:,2)),cell2mat(pupil(:,3))};
        %Do RFT
            y = [pupil{1}; pupil{2}; pupil{3}];
            data = ft_preproc_smooth(y,10); y = data; %Small amount of smoothing for RFT
            [~,out_HN] = RFT_GLM_contrast(X,y,c_HN,'F',[],0); 
            [~,out_SN] = RFT_GLM_contrast(X,y,c_SN,'F',[],0); 
        %Plot pupil   
            emocolors = EmotionColors([1,2,5],:);
            for emo = 1:3
                data = cell2mat(pupil(:,emo));
                X = linspace(0,10,601);
                Y = nanmean(data);
                SEM = nanstd(data)./sqrt(sum(~isnan(data))); %#ok<*NANSTD>
                color = emocolors(emo,:);
                PrettyPatchPlot(X,Y,SEM,color)
            end
        %Plot RFT
            %Significant cluster difference of happy vs. neutral
                ii_cluster = cell2mat(out_HN.clusters.ind(out_HN.clusters.prft<0.05));
                rftH = linspace(X(ii_cluster(1)),X(ii_cluster(end)));
                plot(rftH,0.8*ones(size(rftH)),'color',emocolors(1,:),'LineWidth',4)
            %Significant cluster difference of sad vs. neutral
                ii_cluster = cell2mat(out_SN.clusters.ind(out_SN.clusters.prft<0.05));
                rftS = linspace(X(ii_cluster(1)),X(ii_cluster(end)));
                plot(rftS,0.75*ones(size(rftS)),'color',emocolors(2,:),'LineWidth',4)
        %Layout
            ylim([-0.65 0.9]); ylabel('Pupil diameter (Z)');
            xlim([0,10.05]); xticks(0:10); xlabel('Time since induction onset [s]')

    %Anger/Fear/Neutral
        subplot(4,2,2); hold on; box on
        %Get data
            pupil = physiology_averages_AFN.pupil(:,3:5);
            pupil = pupil(~cellfun(@isempty,pupil(:,1)),:);
            X = kron(eye(3),ones(size(pupil,1),1));
            pupil = {cell2mat(pupil(:,1)),cell2mat(pupil(:,2)),cell2mat(pupil(:,3))};
            pupil = cellfun(@(x)x(:,1:601),pupil,'UniformOutput',false); %crop
        %Do RFT
            y = [pupil{1}; pupil{2}; pupil{3}];
            data = ft_preproc_smooth(y,10); y = data; %Small amount of smoothing for RFT
            [~,out_AN] = RFT_GLM_contrast(X,y,c_AN,'F',[],0); 
            [~,out_FN] = RFT_GLM_contrast(X,y,c_FN,'F',[],0); 
        %Plot pupil   
            emocolors = EmotionColors(3:5,:);
            for emo = 1:3
                data = cell2mat(pupil(:,emo));
                X = linspace(0,10,601);
                Y = nanmean(data);
                SEM = nanstd(data)./sqrt(sum(~isnan(data)));
                color = emocolors(emo,:);
                PrettyPatchPlot(X,Y,SEM,color)
            end
        %Plot RFT
            %Significant cluster difference of anger vs. neutral
                ii_cluster = cell2mat(out_AN.clusters.ind(out_AN.clusters.prft<0.05));
                rftA = linspace(X(ii_cluster(1)),X(ii_cluster(end)));
                plot(rftA,0.8*ones(size(rftH)),'color',emocolors(1,:),'LineWidth',4)
            %Significant cluster difference of fear vs. neutral
                ii_cluster = cell2mat(out_FN.clusters.ind(out_FN.clusters.prft<0.05));
                rftF = linspace(X(ii_cluster(1)),X(ii_cluster(end)));
                plot(rftF,0.75*ones(size(rftF)),'color',emocolors(2,:),'LineWidth',4)
        %Layout
            ylim([-0.65 0.9]); ylabel('Pupil diameter (Z)');
            xlim([0,10.05]); xticks(0:10); xlabel('Time since induction onset [s]')

%% EDA
    %Happy/Sad/Neutral
        subplot(4,2,3); hold on; box on
        %Get data
            EDA = physiology_averages_HSN.EDA(:,[1,2,5]);
            EDA = EDA(~cellfun(@isempty,EDA(:,1)),:);
            X = kron(eye(3),ones(size(EDA,1),1));
            EDA = {cell2mat(EDA(:,1)),cell2mat(EDA(:,2)),cell2mat(EDA(:,3))};
        %Do RFT
            y = cell2mat(EDA');
            [~,out_HN] = RFT_GLM_contrast(X,y,c_HN,'F',[],0); 
            [~,out_SN] = RFT_GLM_contrast(X,y,c_SN,'F',[],0); 
        %Plot EDA
            emocolors = EmotionColors([1,2,5],:);
            for emo = 1:3
                data = cell2mat(EDA(:,emo));
                X = linspace(0,10,101);
                Y = nanmean(data);
                SEM = nanstd(data)./sqrt(sum(~isnan(data)));
                color = emocolors(emo,:);
                PrettyPatchPlot(X,Y,SEM,color)
            end
        %Plot RFT
            %Significant cluster difference of happy vs. neutral
                ii_cluster = cell2mat(out_HN.clusters.ind(out_HN.clusters.prft<0.05));
                rftH = linspace(X(ii_cluster(1)),X(ii_cluster(end)));
                plot(rftH,0.15*ones(size(rftH)),'color',emocolors(1,:),'LineWidth',4)
            %Significant cluster difference of sad vs. neutral
                ii_cluster = cell2mat(out_SN.clusters.ind(out_SN.clusters.prft<0.05));
                rftS = linspace(X(ii_cluster(1)),X(ii_cluster(end)));
                plot(rftS,0.1*ones(size(rftS)),'color',emocolors(2,:),'LineWidth',4)
        %Layout
            xlim([0,10.05]); xticks(0:10); xlabel('Time since induction onset [s]')
            ylim([-1.1 0.25]); ylabel('Skin conductance (Z)');
    %Anger/Fear/Neutral
        subplot(4,2,4); hold on; box on
        %Get data
            EDA = physiology_averages_AFN.EDA(:,3:5);
            EDA = EDA(~cellfun(@isempty,EDA(:,1)),:);
            X = kron(eye(3),ones(size(EDA,1),1));
            EDA = {cell2mat(EDA(:,1)),cell2mat(EDA(:,2)),cell2mat(EDA(:,3))};
        %Do RFT
            y = cell2mat(EDA');
            [~,out_AN] = RFT_GLM_contrast(X,y,c_AN,'F',[],0); 
            [~,out_FN] = RFT_GLM_contrast(X,y,c_FN,'F',[],0); 
        %Plot EDA
            emocolors = EmotionColors(3:5,:);
            for emo = 1:3
                data = cell2mat(EDA(:,emo));
                X = linspace(0,10,101);
                Y = nanmean(data);
                SEM = nanstd(data)./sqrt(sum(~isnan(data)));
                color = emocolors(emo,:);
                PrettyPatchPlot(X,Y,SEM,color)
            end
        %Plot RFT
            %No significant cluster difference of anger vs. neutral
                % ii_cluster = cell2mat(out_AN.clusters.ind(out_AN.clusters.prft<0.05));
                % rftA = linspace(X(ii_cluster(1)+2),X(ii_cluster(end)));
                % plot(rftA,0.15*ones(size(rftA)),'color',emocolors(1,:),'LineWidth',4)
            %No significant cluster difference of fear vs. neutral
                % ii_cluster = cell2mat(out_FN.clusters.ind(out_FN.clusters.prft<0.05));
                % rftF = linspace(X(ii_cluster(1)+2),X(ii_cluster(end)));
                % plot(rftF,0.1*ones(size(rftF)),'color',emocolors(2,:),'LineWidth',4)
            %Layout
                xlim([0,10.05]); xticks(0:10); xlabel('Time since induction onset [s]')
                ylim([-1.1 0.25]); ylabel('Skin conductance (Z)');
      
%% Zygomaticus
    %Happy/Sad/Neutral
        subplot(4,2,5); hold on; box on
        %Get data
            zygomaticus = physiology_averages_HSN.zygomaticus(:,[1,2,5]);
            zygomaticus = zygomaticus(~cellfun(@isempty,zygomaticus(:,1)),:);
            X = kron(eye(3),ones(size(zygomaticus,1),1));
            zygomaticus = {cell2mat(zygomaticus(:,1)),cell2mat(zygomaticus(:,2)),cell2mat(zygomaticus(:,3))};
        %Do RFT
            y = cell2mat(zygomaticus');
            y = ft_preproc_smooth(y,20); %extra smoothing so RFT can be performed
            [~,out_HS] = RFT_GLM_contrast(X,y,c_HS,'F',[],0); 
        %Plot zygomaticus
            emocolors = EmotionColors([1,2,5],:);
            for emo = 1:3
                data = cell2mat(zygomaticus(:,emo));
                data = ft_preproc_smooth(data,20);
                Y = nanmean(data);
                X = linspace(0,10,length(Y));
                SEM = nanstd(data)./sqrt(sum(~isnan(data)));
                color = emocolors(emo,:);
                PrettyPatchPlot(X,Y,SEM,color)
            end
        %Plot RFT: Significant cluster difference of happy vs. sad
            ii_cluster = cell2mat(out_HS.clusters.ind(out_HS.clusters.prft<0.05));
            rftH = linspace(X(ii_cluster(1)),X(ii_cluster(end)));
            plot(rftH,1*ones(size(rftH)),'color',EmotionColors(end,:),'LineWidth',4)
        %Layout
            xlim([0,10.05]); xticks(0:10); xlabel('Time since induction onset [s]')
            ylim([-0.35 1.1]); ylabel('Zygomaticus EMG (Z)');
    %Anger/Fear/Neutral
        subplot(4,2,6); hold on; box on
        %Get data
            zygomaticus = physiology_averages_AFN.zygomaticus(:,3:5);
            zygomaticus = zygomaticus(~cellfun(@isempty,zygomaticus(:,1)),:);
            zygomaticus = {cell2mat(zygomaticus(:,1)),cell2mat(zygomaticus(:,2)),cell2mat(zygomaticus(:,3))};
            emocolors = EmotionColors(3:5,:);
            for emo = 1:length(zygomaticus)
                data = cell2mat(zygomaticus(:,emo));
                data = ft_preproc_smooth(data,20);
                Y = nanmean(data);
                X = linspace(0,10,length(Y));
                SEM = nanstd(data)./sqrt(sum(~isnan(data)));
                color = emocolors(emo,:);
                PrettyPatchPlot(X,Y,SEM,color)
            end
            xlim([0,10.05]); xticks(0:10); xlabel('Time since induction onset [s]')
            ylim([-0.35 1.1]); ylabel('Zygomaticus EMG (Z)');
            
%% Corrugator
    %Happy/Sad/Neutral
        subplot(4,2,7); hold on; box on
        %Get data
            corrugator = physiology_averages_HSN.corrugator(:,[1,2,5]);
            corrugator = corrugator(~cellfun(@isempty,corrugator(:,1)),:);
            X = kron(eye(3),ones(size(corrugator,1),1));
            corrugator = {cell2mat(corrugator(:,1)),cell2mat(corrugator(:,2)),cell2mat(corrugator(:,3))};
        %Do RFT
            y = cell2mat(corrugator');
            y = ft_preproc_smooth(y,20); %extra smoothing so RFT can be performed
            [~,out_HS] = RFT_GLM_contrast(X,y,c_HS,'F',[],0); 
        %Plot corrugator
            emocolors = EmotionColors([1,2,5],:);
            for emo = 1:3
                data = cell2mat(corrugator(:,emo));
                data = ft_preproc_smooth(data,20);
                Y = nanmean(data);
                X = linspace(0,10,length(Y));
                SEM = nanstd(data)./sqrt(sum(~isnan(data)));
                color = emocolors(emo,:);
                PrettyPatchPlot(X,Y,SEM,color)
            end
        %Plot RFT: Significant cluster difference of happy vs. sad
            ii_cluster = cell2mat(out_HS.clusters.ind(out_HS.clusters.prft<0.05));
            rftH = linspace(X(ii_cluster(1)),X(ii_cluster(end)));
            plot(rftH,0.3*ones(size(rftH)),'color',EmotionColors(end,:),'LineWidth',4)
        %Layout
            xlim([0,10.05]); xticks(0:10); xlabel('Time since induction onset [s]')
            ylim([-0.7 0.4]); ylabel('Corrugator EMG (Z)');
    %Anger/Fear/Neutral
        subplot(4,2,8); hold on; box on
        %Get data
            corrugator = physiology_averages_AFN.corrugator(:,3:5);
            corrugator = corrugator(~cellfun(@isempty,corrugator(:,1)),:);
            corrugator = {cell2mat(corrugator(:,1)),cell2mat(corrugator(:,2)),cell2mat(corrugator(:,3))};
            emocolors = EmotionColors(3:5,:);
            for emo = 1:length(corrugator)
                data = cell2mat(corrugator(:,emo));
                data = ft_preproc_smooth(data,20);
                Y = nanmean(data);
                X = linspace(0,10,length(Y));
                SEM = nanstd(data)./sqrt(sum(~isnan(data)));
                color = emocolors(emo,:);
                PrettyPatchPlot(X,Y,SEM,color)
            end
            xlim([0,10.05]); xticks(0:10); xlabel('Time since induction onset [s]')
            ylim([-0.7 0.4]); ylabel('Corrugator EMG (Z)');
    
%% Patch plot
function PrettyPatchPlot(X,Y,E,color) 
    patch([X X(end:-1:1)], [Y+E Y(end:-1:1)-E(end:-1:1)], color, 'facealpha', 0.25, 'Edgecolor', color);
    plot(X, Y, 'linestyle', '-', 'LineWidth', 2, 'color', color);
end