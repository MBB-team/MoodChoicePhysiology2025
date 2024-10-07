
load('physiology_analysis.mat')
scatter(nanmean(binned_physiology.proxy_mood),nanmean(binned_physiology.residuals),40,nanmean(binned_physiology.rated_mood),'filled')
