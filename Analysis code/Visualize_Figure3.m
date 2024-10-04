% Visualises time-resolved physiological responses to emotion inductions.
% Also applies random field theory for cluster-level significance testing.
% (for this, the VBA toolbox must be installed: https://github.com/MBB-team/VBA-toolbox)

% Prepare
    load('physiology_averages.mat')
    EmotionColors = [...
        1.0000    0.7961    0.1373
        0.0510    0.4588    0.9922
        0.8706    0.2196    0.0784
        0.6196    0.0392    0.5216
        0.5000    0.5000    0.5000
        0.2000    0.6000    0.2000];
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
            c_HN = [1;0;-1]; %Happy vs. neutral 
            c_SN = [0;1;-1]; %Sad vs. neutral
            y = [pupil{1}; pupil{2}; pupil{3}];
            data = ft_preproc_smooth(y,10); y = data; %Small amount of smoothing for RFT
            [~,out_HN] = RFT_GLM_contrast(X,y,c_HN,'F',1,0); 
            [~,out_SN] = RFT_GLM_contrast(X,y,c_SN,'F',1,0); 
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
                ii_cluster = out_HN.clusters.ind{1};
                rftH = linspace(X(ii_cluster(1)),X(ii_cluster(end)-6),1000);
                scatter(rftH,0.80*ones(size(rftH)),30,emocolors(1,:),'filled')
            %Significant cluster difference of sad vs. neutral
                ii_cluster = out_SN.clusters.ind{1};
                rftS = linspace(X(ii_cluster(1)),X(ii_cluster(end)-6),1000);
                scatter(rftS,0.75*ones(size(rftS)),30,emocolors(2,:),'filled')
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
            c_AN = [1;0;-1]; %Anger vs. neutral 
            c_FN = [0;1;-1]; %Fear vs. neutral
            y = [pupil{1}; pupil{2}; pupil{3}];
            data = ft_preproc_smooth(y,10); y = data; %Small amount of smoothing for RFT
            [~,out_AN] = RFT_GLM_contrast(X,y,c_AN,'F',1,0); 
            [~,out_FN] = RFT_GLM_contrast(X,y,c_FN,'F',1,0); 
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
                ii_cluster = out_AN.clusters.ind{1};
                rftA = linspace(X(ii_cluster(1)),X(ii_cluster(end)-6),1000);
                scatter(rftA,0.80*ones(size(rftH)),30,emocolors(1,:),'filled')
            %Significant cluster difference of fear vs. neutral
                ii_cluster = out_FN.clusters.ind{2,1};
                rftF = linspace(X(ii_cluster(1)),X(ii_cluster(end)-6),1000);
                scatter(rftF,0.75*ones(size(rftF)),30,emocolors(2,:),'filled')
        %Layout
            ylim([-0.65 0.9]); ylabel('Pupil diameter (Z)');
            xlim([0,10.05]); xticks(0:10); xlabel('Time since induction onset [s]')

return

        subplot(4,2,3); hold on; %EDA
            EDA = {averages.EDA};
            EDA = vertcat(EDA{:});
            X = linspace(0,10,101);
            for emo = 1:3
                data = cell2mat(EDA(:,emo));
                Y = nanmean(data);
                SEM = nanstd(data)./sqrt(sum(~isnan(data)));
                color = emocolors(emo,:);
                PrettyPatchPlot(X,Y,SEM,color)
            end
            rftH = linspace(X(RFT.EDA.Happy(1)),X(RFT.EDA.Happy(end-1)),1000);
            rftS = linspace(X(RFT.EDA.Sad(1)),X(RFT.EDA.Sad(end-1)),1000);
            scatter(rftH,0.1*ones(size(rftH)),30,emocolors(1,:),'filled')
            scatter(rftS,0.05*ones(size(rftS)),30,emocolors(2,:),'filled')
            xlim([0,10.05]); xticks(0:10); %xlabel('Time since induction onset [s]')
            ylim([-1 0.2]); ylabel('Skin conductance (Z)');
            box on; %title('Skin conductance response');
        subplot(4,2,5); hold on; %zygo
            zygo = {averages.zygo};
            zygo = vertcat(zygo{:});
            X = linspace(0,10,501);
            for emo = 1:3
                data = cell2mat(zygo(:,emo));
                Y = nanmean(data);
                SEM = nanstd(data)./sqrt(sum(~isnan(data)));
                color = emocolors(emo,:);
                PrettyPatchPlot(X,Y,SEM,color)
            end
            rftH = linspace(X(RFT.zygo.Happy(1)+3),X(RFT.zygo.Happy(end-5)),1000);
            scatter(rftH,1*ones(size(rftH)),30,emocolors(end,:),'filled')
    %         rftS = linspace(X(RFT.zygo.Sad(1)),X(RFT.zygo.Sad(end-1)),1000);
    %         scatter(rftS,0.05*ones(size(rftS)),30,emocolors(2,:),'filled') %cluster not significant
            xlim([0,10.05]); xticks(0:10); %xlabel('Time since induction onset [s]')
            ylim([-0.35 1.1]); ylabel('Zygomaticus EMG (Z)');
            box on; %title('Zygomaticus activity');
        subplot(4,2,7); hold on; %corru
            corru = {averages.corru};
            corru = vertcat(corru{:});
            X = linspace(0,10,501);
            for emo = 1:3
                data = cell2mat(corru(:,emo));
                Y = nanmean(data);
                SEM = nanstd(data)./sqrt(sum(~isnan(data)));
                color = emocolors(emo,:);
                PrettyPatchPlot(X,Y,SEM,color)
            end
            rftH = linspace(X(RFT.zygo.Happy(1)+3),X(RFT.zygo.Happy(end-5)),1000); %Happy same as sad
            scatter(rftH,0.3*ones(size(rftH)),30,emocolors(end,:),'filled') %Happy same as sad
    %         rftH = linspace(X(RFT.corru.Happy(1)+5),X(RFT.corru.Happy(end-5)),1000);
    %         scatter(rftH,0.3*ones(size(rftH)),30,emocolors(1,:),'filled')
    %         rftS = linspace(X(RFT.corru.Sad(1)),X(RFT.corru.Sad(end-1)),1000);
    %         scatter(rftS,0.05*ones(size(rftS)),30,emocolors(2,:),'filled')
            xlim([0,10.05]); xticks(0:10); xlabel('Time since induction onset [s]')
            ylim([-0.7 0.4]); ylabel('Corrugator EMG (Z)');
            box on; %title('Corrugator activity');       
            
    %% Make figure for A, F, N
        
        %EDA
            subplot(4,2,4); hold on
            %Get data
                EDA = vertcat([all_average_EDA.induction]');
                EDA = {cell2mat(EDA(:,1)),cell2mat(EDA(:,2)),cell2mat(EDA(:,3))};
                EDA = cellfun(@(x)x(:,11:111),EDA,'UniformOutput',false); %crop
            %Do RFT
                y = cell2mat(EDA');
                X = kron(eye(3),ones(size(EDA{1},1),1));
                data = ft_preproc_smooth(y,10); y = data; %Small amount of smoothing for RFT
                [~,out_AN] = RFT_GLM_contrast(X,y,c_AN,'F',1,0); %Not significant at cluster level
                [~,out_FN] = RFT_GLM_contrast(X,y,c_FN,'F',1,0); 
            %Plot EDA
                for emo = 1:3
                    data = cell2mat(EDA(:,emo));
                    X = linspace(0,10,101);
                    Y = nanmean(data);
                    SEM = nanstd(data)./sqrt(sum(~isnan(data)));
                    color = emocolors(emo,:);
                    PrettyPatchPlot(X,Y,SEM,color)
                end
            %Plot RFT
                ii_cluster = out_FN.clusters.ind{1};
                rftF = linspace(X(ii_cluster(1)+2),X(ii_cluster(end)-1),1000);
                scatter(rftF,0.1*ones(size(rftF)),30,emocolors(2,:),'filled')
            xlim([0,10.05]); xticks(0:10); %xlabel('Time since induction onset [s]')
            ylim([-1 0.2]); %ylabel('Signal (z-scored)');
            box on; %title('Skin conductance response');
        %Zygomaticus: A,F
            subplot(4,2,6); hold on
            EMG = vertcat([all_average_EMG.zygomaticus]');
            EMG = vertcat([EMG.induction]');
            EMG = {cell2mat(EMG(:,3)),cell2mat(EMG(:,4)),cell2mat(EMG(:,5))};
            EMG = cellfun(@(x)x(:,151:651),EMG,'UniformOutput',false); %crop
            colors = settings.EmotionColors(3:5,:);
            for i = 1:length(EMG)
                data = EMG{i};
                X = linspace(0,10,501);
                Y = nanmean(data);
                SEM = nanstd(data)./sqrt(sum(~isnan(data)));
                color = colors(i,:);
                PrettyPatchPlot(X,Y,SEM,color)
            end
            xlim([0,10.05]); xticks(0:10); %xlabel('Time since induction onset [s]')
            ylim([-0.35 1.1]); %ylabel('Zygomaticus EMG activity (Z)');
            box on; %title('Zygomaticus activity');
        %Corrugator: A,F
            subplot(4,2,8); hold on
            EMG = vertcat([all_average_EMG.corrugator]');
            EMG = vertcat([EMG.induction]');
            EMG = {cell2mat(EMG(:,3)),cell2mat(EMG(:,4)),cell2mat(EMG(:,5))};
            EMG = cellfun(@(x)x(:,151:651),EMG,'UniformOutput',false); %crop
            colors = settings.EmotionColors(3:5,:);
            for i = 1:length(EMG)
                data = EMG{i};
                X = linspace(0,10,501);
                Y = nanmean(data);
                SEM = nanstd(data)./sqrt(sum(~isnan(data)));
                color = colors(i,:);
                PrettyPatchPlot(X,Y,SEM,color)
            end
            xlim([0,10.05]); xticks(0:10); xlabel('Time since induction onset [s]')
            ylim([-0.7 0.4]); %ylabel('Corrugator EMG activity (Z)');
            box on; %title('Corrugator activity');       

% Save
    saveas(gcf,'TimeResolvedPhysiology.png')
    
%% Visualization
function PrettyPatchPlot(X,Y,E,color) 
    patch([X X(end:-1:1)], [Y+E Y(end:-1:1)-E(end:-1:1)], color, 'facealpha', 0.25, ...
        'Edgecolor', color);
    plot(X, Y, 'linestyle', '-', 'LineWidth', 2, 'color', color);
end