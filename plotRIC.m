function plotRIC(idxConverge,residuals,sigmaX,timePlot,mainTitle)

% Plots a figure of 6 plots (3 position and 3 velocity RIC)

figure()
subplot(3,2,1)
plot(timePlot,residuals(:,1),'r','LineStyle','-.')
posxRMS = rms(residuals(idxConverge:end,1));
hold on
plot(timePlot,3*sigmaX(:,1),'k','LineStyle','-')
hold on
plot(timePlot,-3*sigmaX(:,1),'k','LineStyle','-')
ylabel('R [km]')
legend(strcat('RMS = ',num2str(posxRMS)))
grid on
subplot(3,2,3)
plot(timePlot,residuals(:,2),'r','LineStyle','-.')
posyRMS = rms(residuals(idxConverge:end,2));
hold on
plot(timePlot,3*sigmaX(:,2),'k','LineStyle','-')
hold on
plot(timePlot,-3*sigmaX(:,2),'k','LineStyle','-')
ylabel('I [km]')
legend(strcat('RMS = ',num2str(posyRMS)))
grid on
subplot(3,2,5)
plot(timePlot,residuals(:,3),'r','LineStyle','-.')
poszRMS = rms(residuals(idxConverge:end,3));
hold on
plot(timePlot,3*sigmaX(:,3),'k','LineStyle','-')
hold on
plot(timePlot,-3*sigmaX(:,3),'k','LineStyle','-')
ylabel('C [km]')
xlabel('Time Since Epoch [hrs]')
legend(strcat('RMS = ',num2str(poszRMS)))
sgtitle(mainTitle)
grid on

subplot(3,2,2)
plot(timePlot,residuals(:,4),'r','LineStyle','-.')
posxRMS = rms(residuals(idxConverge:end,4));
hold on
plot(timePlot,3*sigmaX(:,4),'k','LineStyle','-')
hold on
plot(timePlot,-3*sigmaX(:,4),'k','LineStyle','-')
ylabel('R [km/s]')
legend(strcat('RMS = ',num2str(posxRMS)))
grid on
subplot(3,2,4)
plot(timePlot,residuals(:,5),'r','LineStyle','-.')
posyRMS = rms(residuals(idxConverge:end,5));
hold on
plot(timePlot,3*sigmaX(:,5),'k','LineStyle','-')
hold on
plot(timePlot,-3*sigmaX(:,5),'k','LineStyle','-')
ylabel('I [km/s]')
legend(strcat('RMS = ',num2str(posyRMS)))
grid on
subplot(3,2,6)
plot(timePlot,residuals(:,6),'r','LineStyle','-.')
poszRMS = rms(residuals(idxConverge:end,6));
hold on
plot(timePlot,3*sigmaX(:,6),'k','LineStyle','-')
hold on
plot(timePlot,-3*sigmaX(:,6),'k','LineStyle','-')
ylabel('C [km/s]')
xlabel('Time Since Epoch [hrs]')
legend(strcat('RMS = ',num2str(poszRMS)))
sgtitle(mainTitle)
grid on