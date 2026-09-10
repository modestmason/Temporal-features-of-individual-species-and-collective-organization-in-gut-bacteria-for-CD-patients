clear all
% 定义菌种颜色函数


% 寻找波动小而丰度水平高的菌种
%1,8,17,20:R2A;
%5,7,9,10:mA
%2,3,6,12,13,14,15,16:R
% 使用更稳健的textscan读取方式
fid = fopen('Species.csv');
raw_data = textscan(fid, '%s', 'Delimiter', '', 'Whitespace', ''); % 禁用自动分隔
fclose(fid);
species_names = raw_data{1};

%R
N=5530;
n_R=[2, 3, 5, 10, 11, 12, 13, 14];
DATA_variance_R=[];
for i = 1:length(n_R)
    
    TV_R(i,:)=meanfluc(n_R(i));
    E_TV_R(i)=mean(TV_R(i,:));
    str=['CD',num2str(n_R(i)),'.csv'];
    temp=readmatrix(str, 'NumHeaderLines', 1);
    temp(:, 1) = [];
    all_spec{n_R(i)}=temp;
    variance{n_R(i)}=zeros(1,N);
    N_0(i)=0;
    for j=1:N
        E_all_spec{n_R(i)}(j)=mean(all_spec{n_R(i)}(j,:));
        if E_all_spec{n_R(i)}(j)>=0
            variance{n_R(i)}(j)=log(TV_R(i,j));
         else
             variance{n_R(i)}(j)=NaN;
        end
    end
    DATA_variance_R=[DATA_variance_R;variance{n_R(i)}];

end

valid_columns_R = sum(~isnan(DATA_variance_R), 1) >= 6; 
filtered_data_R = DATA_variance_R(:, valid_columns_R);
median_R=median(filtered_data_R,1,'omitnan');
[~, sortedIndices_R] = sort(median_R, 'descend');
top25Indices_R = sortedIndices_R(1:20);
top25Data_R = filtered_data_R(:, top25Indices_R);

original_indices_R = find(valid_columns_R);  
species_indices_R = original_indices_R(top25Indices_R); 
species_labels_R = species_names(species_indices_R);  

figure('Position', [100 100 600 600]);

h = boxplot(fliplr(top25Data_R), 'Labels', flipud(species_labels_R),"Orientation","horizontal");
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

title('Top species by median change in R',...
      'FontSize', 20, 'FontWeight', 'bold');
xlabel('log_1_0(Averaged squared first-order difference)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Species',...
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

ax_color = axes('Position', [0 0 1 1], ...
     'Color', 'none', ...
     'Visible', 'off', ...
     'YLim', get(ax, 'YLim'), ...
     'XLim', get(ax, 'XLim'));

ax_pos = get(ax, 'Position');
x_min = ax.XLim(1);
y_min = ax.YLim(1);
y_max = ax.YLim(2);
y_range = y_max - y_min;

for i = 1:length(yTicks)
    species = yLabels{i};
    color = get_species_color(species); 
 
    y_pos = (yTicks(i) - y_min) / y_range;
    y_norm = ax_pos(2) + y_pos * ax_pos(4);
    

    x_norm = ax_pos(1) - 0.05; 
    

    text(x_norm, y_norm, species, ...
         'FontAngle', 'italic', ...
         'FontSize', 12, ...
         'Color', color, ...
         'HorizontalAlignment', 'right', ...
         'Units', 'normalized', ... 
         'Parent', ax_color);
end


set(ax_color, 'Box', 'off', 'XColor', 'none', 'YColor', 'none');
hold(ax_color, 'off');

ax.Color = [1 1 1]; 
box off  
set(gcf, 'Color', [1 1 1]);  


hold off;

%mA
n_mA=[4, 6, 8, 9];
DATA_variance_mA=[];
for i = 1:length(n_mA)
    TV_mA(i,:)=meanfluc(n_mA(i));
    E_TV_mA(i)=mean(TV_mA(i,:));
    str=['CD',num2str(n_mA(i)),'.csv'];
    temp=readmatrix(str, 'NumHeaderLines', 1);
    temp(:, 1) = [];
    all_spec{n_mA(i)}=temp;
    variance{n_mA(i)}=zeros(1,N);
    for j=1:N
        E_all_spec{n_mA(i)}(j)=mean(all_spec{n_mA(i)}(j,:));
        if E_all_spec{n_mA(i)}(j)>=0
            variance{n_mA(i)}(j)=log(TV_mA(i,j));
        else
            variance{n_mA(i)}(j)=NaN;
        end
    end
    DATA_variance_mA=[DATA_variance_mA;variance{n_mA(i)}];
end

valid_columns_mA = sum(~isnan(DATA_variance_mA), 1) >= 3; 

filtered_data_mA = DATA_variance_mA(:, valid_columns_mA);
median_mA=median(filtered_data_mA,1,'omitnan');
[~, sortedIndices_mA] = sort(median_mA, 'descend');
top25Indices_mA = sortedIndices_mA(1:20);
top25Data_mA = filtered_data_mA(:, top25Indices_mA);

original_indices_mA = find(valid_columns_mA);  
species_indices_mA = original_indices_mA(top25Indices_mA);  
species_labels_mA = species_names(species_indices_mA);  
figure('Position', [100 100 600 600]);

h = boxplot(fliplr(top25Data_mA), 'Labels', flipud(species_labels_mA),"Orientation","horizontal");
boxColors = [18/255, 153/255, 144/255];  
set(findobj(h, 'Tag', 'Box'), 'Color', boxColors, 'LineWidth', 1.5);

medianColors = [9/255, 107/255, 104/255];  
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

title('Top species by median change in mA',...
      'FontSize', 20, 'FontWeight', 'bold');
xlabel('log_1_0(Averaged squared first-order difference)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Species','FontSize', 20, 'FontWeight', 'bold');

ax = gca;
ax.FontSize = 12;
ax.XGrid = 'on';  
ax.GridLineStyle = '--';
ax.GridAlpha = 0.3;
yTicks = get(ax, 'YTick');
yLabels = get(ax, 'YTickLabel');
ax.YAxis.Color = 'w';
ax.YAxis.Label.Color = [0 0 0];


ax_color = axes('Position', [0 0 1 1], ... 
     'Color', 'none', ...
     'Visible', 'off', ... 
     'YLim', get(ax, 'YLim'), ...
     'XLim', get(ax, 'XLim'));

ax_pos = get(ax, 'Position');
x_min = ax.XLim(1);
y_min = ax.YLim(1);
y_max = ax.YLim(2);
y_range = y_max - y_min;

for i = 1:length(yTicks)
    species = yLabels{i};
    color = get_species_color(species); 
    
    y_pos = (yTicks(i) - y_min) / y_range;
    y_norm = ax_pos(2) + y_pos * ax_pos(4);
    
    x_norm = ax_pos(1) - 0.05; 
    
    text(x_norm, y_norm, species, ...
         'FontAngle', 'italic', ...
         'FontSize', 12, ...
         'Color', color, ...
         'HorizontalAlignment', 'right', ...
         'Units', 'normalized', ... 
         'Parent', ax_color);
end

set(ax_color, 'Box', 'off', 'XColor', 'none', 'YColor', 'none');
hold(ax_color, 'off');

ax.Color = [1 1 1];  
box off  
set(gcf, 'Color', [1 1 1]);  

hold off;

%R2A
n_R2A=[1, 7, 15, 16];
DATA_variance_R2A=[];
for i = 1:length(n_R2A)
    TV_R2A(i,:)=meanfluc(n_R2A(i));
    E_TV_R2A(i)=mean(TV_R2A(i,:));
    str=['CD',num2str(n_R2A(i)),'.csv'];
    temp=readmatrix(str, 'NumHeaderLines', 1);
    temp(:, 1) = [];
    all_spec{n_R2A(i)}=temp;
    variance{n_R2A(i)}=zeros(1,N);
    for j=1:N
        E_all_spec{n_R2A(i)}(j)=mean(all_spec{n_R2A(i)}(j,:));
        if E_all_spec{n_R2A(i)}(j)>=0
            variance{n_R2A(i)}(j)=log(TV_R2A(i,j));
        else
            variance{n_R2A(i)}(j)=NaN;
        end
    end
    DATA_variance_R2A=[DATA_variance_R2A;variance{n_R2A(i)}];
end

valid_columns_R2A = sum(~isnan(DATA_variance_R2A), 1) >= 3; 
filtered_data_R2A = DATA_variance_R2A(:, valid_columns_R2A);
median_R2A=median(filtered_data_R2A,1,'omitnan');
[~, sortedIndices_R2A] = sort(median_R2A, 'descend');
top25Indices_R2A = sortedIndices_R2A(1:20);

top25Data_R2A = filtered_data_R2A(:, top25Indices_R2A);
original_indices_R2A = find(valid_columns_R2A); 
species_indices_R2A = original_indices_R2A(top25Indices_R2A); 
species_labels_R2A = species_names(species_indices_R2A);  

figure('Position', [100 100 600 600]);
h = boxplot(fliplr(top25Data_R2A), 'Labels', flipud(species_labels_R2A),"Orientation","horizontal");
boxColors = [218/255, 108/255, 108/255];  
set(findobj(h, 'Tag', 'Box'), 'Color', boxColors, 'LineWidth', 1.5);

medianColors = [166/255, 60/255, 60/255];  
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

title('Top species by median change in R2A',...
      'FontSize', 20, 'FontWeight', 'bold');
xlabel('log_1_0(Averaged squared first-order difference)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Species','FontSize', 20, 'FontWeight', 'bold');

ax = gca;
ax.FontSize = 12;
ax.XGrid = 'on'; 
ax.GridLineStyle = '--';
ax.GridAlpha = 0.3;
yTicks = get(ax, 'YTick');
yLabels = get(ax, 'YTickLabel');
ax.YAxis.Color = 'w';  
ax.YAxis.Label.Color = [0 0 0];

ax_color = axes('Position', [0 0 1 1], ... 
     'Color', 'none', ...
     'Visible', 'off', ... 
     'YLim', get(ax, 'YLim'), ...
     'XLim', get(ax, 'XLim'));

ax_pos = get(ax, 'Position');
x_min = ax.XLim(1);
y_min = ax.YLim(1);
y_max = ax.YLim(2);
y_range = y_max - y_min;

for i = 1:length(yTicks)
    species = yLabels{i};
    color = get_species_color(species); 
    
    y_pos = (yTicks(i) - y_min) / y_range;
    y_norm = ax_pos(2) + y_pos * ax_pos(4);
    
    x_norm = ax_pos(1) - 0.05; 
    
    text(x_norm, y_norm, species, ...
         'FontAngle', 'italic', ...
         'FontSize', 12, ...
         'Color', color, ...
         'HorizontalAlignment', 'right', ...
         'Units', 'normalized', ... 
         'Parent', ax_color);
end

set(ax_color, 'Box', 'off', 'XColor', 'none', 'YColor', 'none');
hold(ax_color, 'off');


ax.Color = [1 1 1];  
box off  
set(gcf, 'Color', [1 1 1]);  

hold off;
