function plotRADEC(idxConverge,store_zResiduals,timePlot,mainTitle)

store_zResiduals = store_zResiduals.*(180*3600)/pi; % [arcsec]

figure()
subplot(2,1,1)
scatter(timePlot,store_zResiduals(1,:),3,'filled')
raRMS = rms(store_zResiduals(1,idxConverge:end));
ylabel('RA [arcsec]')
legend(strcat('RMS = ', num2str(raRMS)))
grid on
subplot(2,1,2)
scatter(timePlot,store_zResiduals(2,:),3,'filled')
decRMS = rms(store_zResiduals(2,idxConverge:end));
ylabel('DEC [arcsec]')
legend(strcat('RMS = ', num2str(decRMS)))
sgtitle(mainTitle)
xlabel('Time Since Epoch [hrs]')
grid on