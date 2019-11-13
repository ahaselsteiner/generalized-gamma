pdTrue = GeneralizedGamma(1.37, 0.98, 1.60);
hsPoints =   [0.5  1    1.5   2    4    6    ]; % from Ochi (1992), Fig. 4b
probPoints = [0.21 0.55 0.70  0.83 0.98 0.995]; % from Ochi (1992), Fig. 4b

% For the distribution shown in Ochi, Fig a, the MLE does not converge.
% pdTrue = GeneralizedGamma(15.77, 0.612, 8.71);
% hsPoints =   [1    2    4     6    8     10    ]; % from Ochi (1992), Fig. 4a
% probPoints = [0.10 0.49 0.90  0.97 0.997 0.9996]; % from Ochi (1992), Fig. 4a

x = [0:0.01:15];
f = pdTrue.pdf(x);
F = pdTrue.cdf(x);

fig1 = figure('position', [100 100 900 300]);
subplot(1, 2, 1);
plot(x, f);
ylabel('Density (-)');
xlabel('Significant wave height (m)');
box off

subplot(1, 2, 2);
pstariCdf = icdf('Normal',F,0,1);
plot(x, pstariCdf, '-k');
hold on
pstariPoints = icdf('Normal',probPoints,0,1);
plot(hsPoints, pstariPoints, 'xk');
set(gca, 'xtick', [1 2 4 6 8 10 15]);
set(gca, 'xscale', 'log');
pTicks = [0.05 0.10 0.20 0.5 0.8 0.9 0.95 0.99 0.999 0.9999];
set(gca, 'ytick',  icdf('Normal',pTicks, 0, 1));
set(gca, 'yticklabels', pTicks)
xlim([0.5 10])
ylim([icdf('Normal',[pTicks(1) pTicks(end)], 0, 1)])
ylabel('Probability');
xlabel('Significant wave height (m)');
grid on
suptitle(['Compare this graph with Ochi (1992), Fig. 4b ' ...
     '(doi: 10.1061/9780872629332.038).']);


n = 5000;
nOfSamples = 20;

lambdaEstimated = nan(nOfSamples, 1);
cEstimated = nan(nOfSamples, 1);
mEstimated = nan(nOfSamples, 1);
for i = 1:nOfSamples
    sample = pdTrue.drawSample(n);
    pdEstimated(i) = GeneralizedGamma();
    pdEstimated(i).fitDist(sample);
    lambdaEstimated(i) = pdEstimated(i).Lambda;
    cEstimated(i) = pdEstimated(i).C;
    mEstimated(i) = pdEstimated(i).M;
end

fig2 = figure('position', [100 100 500, 230]);
subplot(1, 3, 1)
hold on
plot([0.5 1.5], [pdTrue.Lambda pdTrue.Lambda], '-k')
boxplot(lambdaEstimated, {'$$\hat{\lambda}$$'})
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
text(1.15, mean(lambdaEstimated), [num2str(mean(lambdaEstimated), '%1.3f') '+-' ...
    num2str(std(lambdaEstimated), '%1.3f')], 'fontsize', 8, ...
    'verticalalignment', 'bottom'); 
box off

subplot(1, 3, 2)
hold on
plot([0.5 1.5], [pdTrue.C pdTrue.C], '-k')
boxplot(cEstimated, {'$$\hat{c}$$'})
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
text(1.15, mean(cEstimated), [num2str(mean(cEstimated), '%1.3f') '+-' ...
    num2str(std(cEstimated), '%1.3f')], 'fontsize', 8, ...
    'verticalalignment', 'bottom');
box off

subplot(1, 3, 3)
hold on
plot([0.5 1.5], [pdTrue.M pdTrue.M], '-k')
boxplot(mEstimated, {'$$\hat{m}$$'})
bp = gca;
bp.XAxis.TickLabelInterpreter = 'latex';
text(1.15, mean(mEstimated), [num2str(mean(mEstimated), '%1.3f') '+-' ...
    num2str(std(mEstimated), '%1.3f')], 'fontsize', 8, ...
    'verticalalignment', 'bottom'); 
box off
suptitle(['Parameter estimation, true parameters: ' ...
    num2str(pdTrue.Lambda) ', ' num2str(pdTrue.C) ', ' num2str(pdTrue.M)]);