
load('physiology_analysis.mat')

figure
subplot(2,2,1); hold on; box on
    R = nanmean(physiology_correlations,3);
    P = NaN(size(R));
    for i = 1:length(R)
        R(i,i) = NaN; 
        for j = 1:length(R)
            [~,P(i,j)] = ttest(squeeze(physiology_correlations(i,j,:)));
        end
    end
    imagesc(R);
    set(gca,'YDir','reverse')
    axis square
    xticks(1:length(R)); yticks(1:length(R))
    colorbar
subplot(4,2,2); hold on
    RH_Patchplot(nanmean(binned_physiology.rated_mood),binned_physiology.proxy_arousal);
subplot(4,2,4); hold on
    RH_Patchplot(nanmean(binned_physiology.rated_mood),binned_physiology.proxy_valence);
subplot(2,2,3); hold on
    betas_moodProxy = cell2mat(beta_mood_residuals(:,2));
    Y = betas_moodProxy(:,1) + betas_moodProxy(:,2).*linspace(-0.3,0.3);
    RH_Patchplot(linspace(-0.3,0.3),Y);
    errorbar(nanmean(binned_physiology.proxy_mood),nanmean(binned_physiology.residuals),RH_SEM(binned_physiology.residuals),'ko')
    scatter(nanmean(binned_physiology.proxy_mood),nanmean(binned_physiology.residuals),40,nanmean(binned_physiology.rated_mood),'filled')
    xlim([-0.35,0.35])
subplot(2,2,4); hold on
    betas_ratedMood = cell2mat(beta_mood_residuals(:,1));
    [~,P_res] = ttest([betas_ratedMood(:,2),betas_moodProxy(:,2)])
    scatter(betas_ratedMood(:,2),betas_moodProxy(:,2),30,'filled')
    lsline
    [R1,P1,R2,P2] = RH_Corr([betas_ratedMood(:,2),betas_moodProxy(:,2)])