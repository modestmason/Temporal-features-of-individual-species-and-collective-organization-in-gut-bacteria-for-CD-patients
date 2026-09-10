clear all

dataR=[0.073426573	0.079545455	0.07079646	0.07063197	0.05	0.039840637	0.055749129	0.073446328
0.034383954	0.04691358	0.034031414	0.044642857	0.037735849	0.024561404	0.036065574	0.047835991
0.057877814	0.080357143	0.085889571	0.021857923	0.054794521	0.038022814	0.057142857	0.053892216
0.042424242	0.072727273	0.063492063	0.047945205	0.05019305	0.036363636	0.078767123	0.057416268
0.050660793	0.071274298	0.064444444	0.02970297	0.050793651	0.027322404	0.055555556	0.060851927
0.038116592	0.038116592	0.038288288	0.032258065	0.037414966	0.020997375	0.030456853	0.026694045
0.076190476	0.083636364	0.105263158	0.052631579	0.074257426	0.035	0.101851852	0.063122924
0.041666667	0.049149338	0.070140281	0.069767442	0.068062827	0.027638191	0.065022422	0.061258278
];
median_R=median(dataR,1,'omitnan');
[~, sortedIndices_R] = sort(median_R, 'descend');
top25Indices_R = sortedIndices_R(1:8);
top25Data_R = dataR(:, top25Indices_R);
species_labels_R = {'R-3';'R-5';'R-14';'R-13';'R-11';'R-2';'R-10';'R-12'};  


figure('Position', [100 100 700 300]);

h = boxplot(fliplr(top25Data_R), 'Labels', flipud(species_labels_R));
boxColors = [0/255, 101/255, 248/255];  
set(findobj(h, 'Tag', 'Box'), 'Color', boxColors, 'LineWidth', 1.5);

medianColors = [39/255, 84/255, 138/255]; 
medianLines = findobj(h, 'Tag', 'Median');
for i = 1:numel(medianLines)
    set(medianLines(i), 'Color', medianColors, 'LineWidth', 2.2);
end

UpperWhiskers = findobj(h, 'Tag', 'Upper Whisker');
for i = 1:numel(UpperWhiskers)
    set(UpperWhiskers(i), 'Color', boxColors, 'LineWidth', 1,'LineStyle','-');
end
LowerWhiskers = findobj(h, 'Tag', 'Lower Whisker');
for i = 1:numel(LowerWhiskers)
    set(LowerWhiskers(i), 'Color', boxColors, 'LineWidth', 1,'LineStyle','-');
end

caps = findobj(h, 'Tag', 'Upper Adjacent Value');
set(caps, 'Visible', 'off');
caps = findobj(h, 'Tag', 'Lower Adjacent Value');
set(caps, 'Visible', 'off');
outliers = findobj(h, 'Tag', 'Outliers');
set(outliers, 'MarkerEdgeColor', medianColors);

title('Network similarties with the networks of patients in mA/R2A',...
      'FontSize', 20, 'FontWeight', 'bold');
xlabel('Patient', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Network smilarity',...
       'FontSize', 20, 'FontWeight', 'bold');

ax = gca;
ax.FontSize = 12;
ax.XGrid = 'on'; 
ax.GridLineStyle = '--';
ax.GridAlpha = 0.3;
yTicks = get(ax, 'YTick');
yLabels = get(ax, 'YTickLabel');
ax.YAxis.Color = 'w';  
ax.YAxis.Label.Color = [0 0 0];
axes('Position', get(ax, 'Position'),...
     'Color', 'none',...
     'YAxisLocation', 'left',...
     'XTick', [],...
     'YLim', get(ax, 'YLim'));

set(gca, 'YTick', yTicks,...
         'YTickLabel', yLabels,...
         'FontSize', 12);

ax.Color = [1 1 1];
box off  
set(gcf, 'Color', [1 1 1]); 
hold off;