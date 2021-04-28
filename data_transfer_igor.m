%% Data transfer Excel - Igor

% Gramm
addpath('C:\Users\lena_\Projects\code_extern\Matlab\gramm-master\');

% Colormap
rgbmap = [204, 0, 4;...
    77, 0, 153;...
    0, 153, 77;...
    0, 102, 153]; % BASE; MINI; A184V; A184S
custommap = rgbmap./256;

% Load
data = readtable('C:\Users\lena_\OneDrive\Dokumente\Masterthesis\Analysis_Master\Confocal_tbls\Summary_Confocal.xlsx', 'Sheet',4);
% data = readtable('C:\Users\lena_\OneDrive\Dokumente\Masterthesis\Analysis_Master\Confocal_tbls\FOV_bounds_all.xlsx');

load('C:\Users\lena_\OneDrive\Dokumente\Masterthesis\Analysis_Master\Confocal_tbls\summary_Confocal.mat');

%% 2D Histogram
sensor ='A184V';
cond = 'MINI';
% val1 = data.Area_um_Mean(strcmp(data.Sensor, sensor) & strcmp(data.BASE_MINI, cond));
% val1 = data.ToT_Num_Events(strcmp(data.Sensor, sensor) & strcmp(data.BASE_MINI, cond));
val1 = data.Events_per_Dend_Area(strcmp(data.Sensor, sensor) & strcmp(data.BASE_MINI, cond));
val2 = data.rel_ML_pos(strcmp(data.Sensor, sensor) & strcmp(data.BASE_MINI, cond));

min1 = min(val1); max1 = max(val1);
range1 = [0 0.4];
min2 = min(val2); max2 = max(val2);
range2 = [0 100];
nbins1 = 10; nbins2 = 10;
edges1 = linspace(range1(1), range1(2), nbins1+1);
edges2 = linspace(range2(1), range2(2), nbins2+1);
% count: ML pos along y, ROI area along x
jointbincount = NaN(nbins, nbins);
for iR = 1:nbins2
    for iC = 1:nbins1
        tmpval2 =val2 >= edges2(iR) & val2 < edges2(iR+1);
        jointbincount(iR,iC) = sum(tmpval2 & val1 >= edges1(iC) & val1 < edges1(iC+1));
    end
end
jointbincount = jointbincount';

norm_jointbincount = jointbincount./max(jointbincount);
norm_jointbincount(isnan(norm_jointbincount)) = 0;
norm_jointbincount=norm_jointbincount.*100;

%% Split & summarize data
borders = [0 20 40 60 80 100];
% [0 33 66 100]
splitvar = 'rel_ML_pos'; splitpter = find(strcmp(data.Properties.VariableNames, splitvar)==1);
% rel_ML_pos
sumvar = 'dFoF_Amp'; sumpter = find(strcmp(data.Properties.VariableNames, sumvar)==1);
% Num_ROIs Area_um_Mean Amp_dFoF_Mean ToT_Num_Events
% Num_Events dFoF_Amp EvRate Ev_struct_norm_area Area
splitsum_tbl = table();
splitsum_tbl.A184V_BASE = NaN(numel(borders)-1,1);
splitsum_tbl.A184V_MINI = NaN(numel(borders)-1,1);
splitsum_tbl.A18S_BASE = NaN(numel(borders)-1,1);
splitsum_tbl.A184S_MINI = NaN(numel(borders)-1,1);
splitsum_tbl.A184V = NaN(numel(borders)-1,1);
splitsum_tbl.A184S = NaN(numel(borders)-1,1);
splitci1_tbl = splitsum_tbl; splitci2_tbl = splitci1_tbl;
splitci_tbl = splitsum_tbl;
splitn_tbl = splitsum_tbl;
anova_val = [];
anova_grp1 = [];
anova_grp2 = [];

for i = 1:numel(borders)-1
%     tmpdata = data(table2array(data(:,splitpter)) >= borders(i) & table2array(data(:,splitpter)) < borders(i+1),:);
    tmpdata = data(strcmp(data.Sensor, 'A184V') & data.Area >= 50,:);
    
    tmp = table2array(tmpdata(strcmp(tmpdata.Sensor, 'A184V') & strcmp(tmpdata.BASE_MINI, 'BASE'),sumpter));
    splitsum_tbl.A184V_BASE(i) = mean(tmp);
    splitci_tbl.A184V_BASE(i) = 1.96* std(tmp)/sqrt(numel(tmp));
    splitci1_tbl.A184V_BASE(i) = mean(tmp) - 1.96* std(tmp)/sqrt(numel(tmp));
    splitci2_tbl.A184V_BASE(i) = mean(tmp) + 1.96* std(tmp)/sqrt(numel(tmp));
%     splitn_tbl.A184V_BASE(i) = numel(tmp);
    splitn_tbl.A184V_BASE(i) = numel(unique(tmpdata.Recording(strcmp(tmpdata.Sensor, 'A184V')& strcmp(tmpdata.BASE_MINI, 'BASE'))));
    anova_val = [anova_val; tmp];
    anova_grp1 = [anova_grp1; repelem(1, size(tmp,1))'];
    anova_grp2 = [anova_grp2; repelem(i, size(tmp,1))'];
    
    tmp = table2array(tmpdata(strcmp(tmpdata.Sensor, 'A184V') & strcmp(tmpdata.BASE_MINI, 'MINI'),sumpter));
    splitsum_tbl.A184V_MINI(i) = mean(tmp);
    splitci_tbl.A184V_MINI(i) = 1.96* std(tmp)/sqrt(numel(tmp));
    splitci1_tbl.A184V_MINI(i) = mean(tmp) - 1.96* std(tmp)/sqrt(numel(tmp));
    splitci2_tbl.A184V_MINI(i) = mean(tmp) + 1.96* std(tmp)/sqrt(numel(tmp));
%     splitn_tbl.A184V_MINI(i) = numel(tmp);
    splitn_tbl.A184V_MINI(i)= numel(unique(tmpdata.Recording(strcmp(tmpdata.Sensor, 'A184V')& strcmp(tmpdata.BASE_MINI, 'MINI'))));
    anova_val = [anova_val; tmp];
    anova_grp1 = [anova_grp1; repelem(2, size(tmp,1))'];
    anova_grp2 = [anova_grp2; repelem(i, size(tmp,1))'];
    
    tmp = table2array(tmpdata(strcmp(tmpdata.Sensor, 'A184S') & strcmp(tmpdata.BASE_MINI, 'BASE'),sumpter));
    splitsum_tbl.A18S_BASE(i) = mean(tmp);
    splitci_tbl.A18S_BASE(i) = 1.96* std(tmp)/sqrt(numel(tmp)); 
    splitci1_tbl.A18S_BASE(i) = mean(tmp) - 1.96* std(tmp)/sqrt(numel(tmp)); 
    splitci2_tbl.A18S_BASE(i) = mean(tmp) + 1.96* std(tmp)/sqrt(numel(tmp));
%     splitn_tbl.A18S_BASE(i) = numel(tmp);
    splitn_tbl.A18S_BASE(i)= numel(unique(tmpdata.Recording(strcmp(tmpdata.Sensor, 'A184S')& strcmp(tmpdata.BASE_MINI, 'BASE'))));
    
    tmp = table2array(tmpdata(strcmp(tmpdata.Sensor, 'A184S') & strcmp(tmpdata.BASE_MINI, 'MINI'),sumpter));
    splitsum_tbl.A184S_MINI(i) = mean(tmp);
    splitci_tbl.A184S_MINI(i) = 1.96* std(tmp)/sqrt(numel(tmp)); 
    splitci1_tbl.A184S_MINI(i) = mean(tmp) - 1.96* std(tmp)/sqrt(numel(tmp)); 
    splitci2_tbl.A184S_MINI(i) = mean(tmp) + 1.96* std(tmp)/sqrt(numel(tmp));
%     splitn_tbl.A184S_MINI(i) = numel(tmp);
    splitn_tbl.A184S_MINI(i)= numel(unique(tmpdata.Recording(strcmp(tmpdata.Sensor, 'A184S')& strcmp(tmpdata.BASE_MINI, 'MINI'))));
    
    tmp = table2array(tmpdata(strcmp(tmpdata.Sensor, 'A184V'),sumpter));
    splitsum_tbl.A184V(i) = mean(tmp);
    splitci_tbl.A184V(i) = 1.96* std(tmp)/sqrt(numel(tmp));
    splitn_tbl.A184V(i) = numel(tmp);
    
    tmp = table2array(tmpdata(strcmp(tmpdata.Sensor, 'A184S'),sumpter));
    splitsum_tbl.A184S(i) = mean(tmp);
    splitci_tbl.A184S(i) = 1.96* std(tmp)/sqrt(numel(tmp));
    splitn_tbl.A184S(i) = numel(tmp);
end

% Graph
xticks = borders(1:end-1) + diff(borders)/2;
g(1,1) = gramm('x', xticks, 'y', splitsum_tbl.A184V_BASE,...
    'ymin', splitci1_tbl.A184V_BASE,'ymax', splitci2_tbl.A184V_BASE);
g(1,1).geom_bar('dodge',0.1,'width',0.8);
g(1,1).geom_interval('geom','black_errorbar','dodge',0.1,'width',1);
g(1,1).axe_property('XLIM', [0 100], 'XTick', borders(2:end));
g(1,1).set_color_options('map', custommap(1,:));

g(2,1) = gramm('x', xticks, 'y', splitsum_tbl.A184V_MINI,...
    'ymin', splitci1_tbl.A184V_MINI,'ymax', splitci2_tbl.A184V_MINI);
g(2,1).geom_bar('dodge',0.1,'width',0.8);
g(2,1).geom_interval('geom','black_errorbar','dodge',0.1,'width',1);
g(2,1).axe_property('XLIM', [0 100], 'XTick', borders(2:end));
g(2,1).set_color_options('map', custommap(2,:));

g.set_names('x', 'Rel. ML position', 'y', sumvar);
g.axe_property('YLIM', [0 max([splitci2_tbl.A184V_BASE; splitci2_tbl.A184V_MINI])])
g.set_title(strcat(sumvar, '(Mean & CI)')); 
g.draw();

%% Find means in certain ranges
tmp = table2array(data(strcmp(data.Sensor, 'A184V') & strcmp(data.BASE_MINI, 'MINI') & data.Area >= 50,:));
%& data.rel_ML_pos < 40   & strcmp(data.BASE_MINI, 'MINI') 
mean(tmp, 'omitnan')
numel(tmp)

%% Graph n
xticks = borders(1:end-1) + diff(borders)/2;
gn(1,1) = gramm('x', xticks, 'y', splitn_tbl.A184V_BASE);
gn(1,1).geom_bar('dodge',0.1,'width',0.8);
gn(1,1).axe_property('XLIM', [0 100], 'XTick', borders(2:end));
gn(1,1).set_color_options('map', custommap(1,:)); 

gn(2,1) = gramm('x', xticks, 'y', splitn_tbl.A184V_MINI);
gn(2,1).geom_bar('dodge',0.1,'width',0.8);
gn(2,1).axe_property('XLIM', [0 100], 'XTick', borders(2:end));
gn(2,1).set_color_options('map', custommap(2,:));

gn.set_names('x', 'Rel. ML position', 'y', '# units');
gn.axe_property('YLIM', [0 max([splitn_tbl.A184V_BASE; splitn_tbl.A184V_MINI])])
gn.set_title('Units'); 
gn.draw();


%% Graph summary statistics
% Summary stats
sum_tbl = table();
sum_tbl.sensor(1:4) = {'A184V','A184V','A184S','A184S'};
sum_tbl.cond(1:4) = {'BASE','MINI','BASE','MINI'};
anova_val = [];
anova_grp1 = [];
anova_grp2 = [];

for iC = 1:4
    tmpdata = data(strcmp(data.Sensor, sum_tbl.sensor(iC)) & strcmp(data.BASE_MINI, sum_tbl.cond(iC)),:);
    sum_tbl.dFoF_amp_mean(iC) = mean(tmpdata.Amp_dFoF_Mean);
    sum_tbl.dFoF_amp_ci_hspan(iC) = 1.96*std(tmpdata.Amp_dFoF_Mean)/sqrt(numel(tmpdata.Amp_dFoF_Mean));
    sum_tbl.roi_area_mean(iC) = mean(tmpdata.Area_um_Mean);
    sum_tbl.roi_area_ci_hspan(iC) = 1.96*std(tmpdata.Area_um_Mean)/sqrt(numel(tmpdata.Area_um_Mean));
    sum_tbl.num_rois_mean(iC) = mean(tmpdata.Num_ROIs);
    sum_tbl.num_rois_ci_hspan(iC) = 1.96*std(tmpdata.Num_ROIs)/sqrt(numel(tmpdata.Num_ROIs));
    sum_tbl.num_ev_mean(iC) = mean(tmpdata.ToT_Num_Events);
    sum_tbl.num_ev_ci_hspan(iC) = 1.96*std(tmpdata.ToT_Num_Events)/sqrt(numel(tmpdata.ToT_Num_Events));
    sum_tbl.num_evproi_mean(iC) = mean(tmpdata.Num_Events_Mean);
    sum_tbl.num_evproi_ci(iC) = 1.96*std(tmpdata.Num_Events_Mean)/sqrt(numel(tmpdata.Num_Events_Mean));
    sum_tbl.evperdend_mean(iC) = mean(tmpdata.Events_per_Dend_Area);
    sum_tbl.evperdend_ci_hspan(iC) = 1.96*std(tmpdata.Events_per_Dend_Area)/sqrt(numel(tmpdata.Events_per_Dend_Area));
    sum_tbl.dend_area_mean(iC) = mean(tmpdata.Dend_Area_um);
    sum_tbl.dend_area_ci_hspan(iC) = 1.96*std(tmpdata.Dend_Area_um)/sqrt(numel(tmpdata.Dend_Area_um));
    sum_tbl.norm_roi_area_mean(iC) = mean(tmpdata.StructNorm_Area);
    sum_tbl.norm_roi_area_ci_hspan(iC) = 1.96*std(tmpdata.StructNorm_Area)/sqrt(numel(tmpdata.StructNorm_Area));
%     sum_tbl.ev_per_norm_roi_area_mean(iC) = mean(tmpdata.ToT_Num_Events*100/tmpdata.StructNorm_Area);
%     sum_tbl.norm_roi_area_ci_hspan(iC) = 1.96*std(tmpdata.StructNorm_Area)/sqrt(numel(tmpdata.StructNorm_Area));
    anova_val = [anova_val; tmpdata.Num_Events_Mean];
    anova_grp1 = [anova_grp1; repelem(sum_tbl.sensor(iC), size(tmpdata,1))'];
    anova_grp2 = [anova_grp2; repelem(sum_tbl.cond(iC), size(tmpdata,1))'];
end

% dF/F Amplitude
g = gramm('x', sum_tbl.sensor, 'y', sum_tbl.dFoF_amp_mean,...
    'ymin', sum_tbl.dFoF_amp_mean-sum_tbl.dFoF_amp_ci_hspan,'ymax', sum_tbl.dFoF_amp_mean+sum_tbl.dFoF_amp_ci_hspan,'color', sum_tbl.cond);
g.geom_bar('dodge',0.8,'width',0.6);
g.geom_interval('geom','black_errorbar','dodge',0.8,'width',1);
g.set_color_options('map', custommap);
g.set_names('x', 'Sensor', 'y', 'dF/F amplitude', 'color', 'Condition');
g.set_title('dF/F amplitude'); 
g.draw();

% ROI Area
g = gramm('x', sum_tbl.sensor, 'y', sum_tbl.roi_area_mean,...
    'ymin', sum_tbl.roi_area_mean-sum_tbl.roi_area_ci_hspan,'ymax', sum_tbl.roi_area_mean+sum_tbl.roi_area_ci_hspan,'color', sum_tbl.cond);
g.geom_bar('dodge',0.8,'width',0.6);
g.geom_interval('geom','black_errorbar','dodge',0.8,'width',1);
g.set_color_options('map', custommap);
g.set_names('x', 'Sensor', 'y', 'Average ROI area [µm²]', 'color', 'Condition');
g.set_title('ROI Area'); 
g.draw();

% # ROIs
g = gramm('x', sum_tbl.sensor, 'y', sum_tbl.num_rois_mean,...
    'ymin', sum_tbl.num_rois_mean-sum_tbl.num_rois_ci_hspan,'ymax', sum_tbl.num_rois_mean+sum_tbl.num_rois_ci_hspan,'color', sum_tbl.cond);
g.geom_bar('dodge',0.8,'width',0.6);
g.geom_interval('geom','black_errorbar','dodge',0.8,'width',1);
g.set_color_options('map', custommap);
g.set_names('x', 'Sensor', 'y', '# ROIs', 'color', 'Condition');
g.set_title('Average Number ROIs'); 
g.draw();

% # events
g = gramm('x', sum_tbl.sensor, 'y', sum_tbl.num_ev_mean,...
    'ymin', sum_tbl.num_ev_mean-sum_tbl.num_ev_ci_hspan,'ymax', sum_tbl.num_ev_mean+sum_tbl.num_ev_ci_hspan,'color', sum_tbl.cond);
g.geom_bar('dodge',0.8,'width',0.6);
g.geom_interval('geom','black_errorbar','dodge',0.8,'width',1);
g.set_color_options('map', custommap);
g.set_names('x', 'Sensor', 'y', '# events', 'color', 'Condition');
g.set_title('Average Number events'); 
g.draw();

% Dendritic area
g = gramm('x', sum_tbl.sensor, 'y', sum_tbl.dend_area_mean,...
    'ymin', sum_tbl.dend_area_mean-sum_tbl.dend_area_ci_hspan,'ymax', sum_tbl.dend_area_mean+sum_tbl.dend_area_ci_hspan,'color', sum_tbl.cond);
g.geom_bar('dodge',0.8,'width',0.6);
g.geom_interval('geom','black_errorbar','dodge',0.8,'width',1);
g.set_color_options('map', custommap);
g.set_names('x', 'Sensor', 'y', 'Average dendritic area [µm²]', 'color', 'Condition');
g.set_title('Dendritic area'); 
g.draw();

% Structure normalized ROI area
g = gramm('x', sum_tbl.sensor, 'y', sum_tbl.norm_roi_area_mean,...
    'ymin', sum_tbl.norm_roi_area_mean-sum_tbl.norm_roi_area_ci_hspan,'ymax', sum_tbl.norm_roi_area_mean+sum_tbl.norm_roi_area_ci_hspan,'color', sum_tbl.cond);
g.geom_bar('dodge',0.8,'width',0.6);
g.geom_interval('geom','black_errorbar','dodge',0.8,'width',1);
g.set_color_options('map', custommap);
g.set_names('x', 'Sensor', 'y', 'Average ROI area normalized by dendritic area', 'color', 'Condition');
g.set_title('Structure-normalized ROI area'); 
g.draw();


% % # Events per dendritic area
% g = gramm('x', sum_tbl.sensor, 'y', sum_tbl.evperdend_mean,...
%     'ymin', sum_tbl.evperdend_mean-sum_tbl.evperdend_ci_hspan,'ymax', sum_tbl.evperdend_mean+sum_tbl.evperdend_ci_hspan,'color', sum_tbl.cond);
% g.geom_bar('dodge',0.8,'width',0.6);
% g.geom_interval('geom','black_errorbar','dodge',0.8,'width',1);
% g.set_color_options('map', custommap);
% g.set_names('x', 'Sensor', 'y', 'Average events per dendritic area', 'color', 'Condition');
% g.set_title('Events per dendritic area'); 
% g.draw();



%% FOV coverage of recordings
borders = 0:5:110;
% [0 33 66 100] [0 20 40 60 80 100]
splitvar1 = 'max_rel_dist'; splitpter1 = find(strcmp(data.Properties.VariableNames, splitvar1)==1);
splitvar2 = 'min_rel_dist'; splitpter2 = find(strcmp(data.Properties.VariableNames, splitvar2)==1);
splitn_noact_tbl = table();
splitn_noact_tbl.A184V_BASE = NaN(numel(borders)-1,1);
splitn_noact_tbl.A184V_MINI = NaN(numel(borders)-1,1);
splitn_noact_tbl.A18S_BASE = NaN(numel(borders)-1,1);
splitn_noact_tbl.A184S_MINI = NaN(numel(borders)-1,1);
splitn_noact_tbl.A184V = NaN(numel(borders)-1,1);
splitn_noact_tbl.A184S = NaN(numel(borders)-1,1);


for i = 1:numel(borders)-1
    tmpdata = data(table2array(data(:,splitpter1)) >= borders(i) & table2array(data(:,splitpter2)) <= borders(i+1),:);
    tmpdata = tmpdata(tmpdata.Activity == 1,:);
    
    tmp = table2array(tmpdata(strcmp(tmpdata.Sensor, 'A184V') & strcmp(tmpdata.BASE_MINI, 'BASE'),5));
    splitn_noact_tbl.A184V_BASE(i) = numel(tmp);
    
    tmp = table2array(tmpdata(strcmp(tmpdata.Sensor, 'A184V') & strcmp(tmpdata.BASE_MINI, 'MINI'),5));
    splitn_noact_tbl.A184V_MINI(i) = numel(tmp);
    
    tmp = table2array(tmpdata(strcmp(tmpdata.Sensor, 'A184S') & strcmp(tmpdata.BASE_MINI, 'BASE'),5));
    splitn_noact_tbl.A18S_BASE(i) = numel(tmp);
    
    tmp = table2array(tmpdata(strcmp(tmpdata.Sensor, 'A184S') & strcmp(tmpdata.BASE_MINI, 'MINI'),5));
    splitn_noact_tbl.A184S_MINI(i) = numel(tmp);
end

% Graph n
xticks = borders(1:end-1) + diff(borders)/2;
gn(1,1) = gramm('x', xticks, 'y', splitn_noact_tbl.A184V_BASE);
gn(1,1).geom_bar('dodge',0.1,'width',0.8);
gn(1,1).set_color_options('map', custommap(1,:)); 

gn(2,1) = gramm('x', xticks, 'y', splitn_noact_tbl.A184V_MINI);
gn(2,1).geom_bar('dodge',0.1,'width',0.8);
gn(2,1).set_color_options('map', custommap(2,:));

gn.set_names('x', 'Rel. ML position', 'y', '# units');
gn.axe_property('XLIM', [0 110], 'XTick', borders(2:end));
gn.axe_property('YLIM', [0 max([splitn_noact_tbl.A184V_BASE; splitn_noact_tbl.A184V_MINI])])
gn.set_title('Silent units covering ML partitions'); 
gn.draw();

% Histogram
figure, histogram(data.min_rel_dist,40);
hold on, histogram(data.max_rel_dist,40);
legend('Lower bound of FOVs','Upper bound of FOVs');
title('ML coverage of recordings');
xlabel('Rel. ML position'); ylabel('# units');


%% 2PM data
data2pm = readtable('C:\Users\lena_\OneDrive\Dokumente\Masterthesis\Analysis_Master\2PM_tbls\Summary_2PM.xlsx', 'Sheet',3);

% Summary stats
sum_tbl = table();
sum_tbl.sensor(1:2) = {'A184V','A184S'};

for iC = 1:2
    tmpdata = data2pm(strcmp(data2pm.Sensor, sum_tbl.sensor(iC)),:);
    sum_tbl.dFoF_amp_mean(iC) = mean(tmpdata.Amp_dFoF_Mean);
    sum_tbl.dFoF_amp_ci_hspan(iC) = 1.96*std(tmpdata.Amp_dFoF_Mean)/sqrt(numel(tmpdata.Amp_dFoF_Mean));
    sum_tbl.roi_area_mean(iC) = mean(tmpdata.Area_um_Mean);
    sum_tbl.roi_area_ci_hspan(iC) = 1.96*std(tmpdata.Area_um_Mean)/sqrt(numel(tmpdata.Area_um_Mean));
    sum_tbl.num_rois_mean(iC) = mean(tmpdata.Num_ROIs);
    sum_tbl.num_rois_ci_hspan(iC) = 1.96*std(tmpdata.Num_ROIs)/sqrt(numel(tmpdata.Num_ROIs));
    sum_tbl.num_ev_mean(iC) = mean(tmpdata.ToT_Num_Events);
    sum_tbl.num_ev_ci_hspan(iC) = 1.96*std(tmpdata.ToT_Num_Events)/sqrt(numel(tmpdata.ToT_Num_Events));
    sum_tbl.num_ev_mean_proi(iC) = mean(tmpdata.Num_Events_Mean);
    sum_tbl.num_ev_proi_ci_hspan(iC) = 1.96*std(tmpdata.Num_Events_Mean)/sqrt(numel(tmpdata.Num_Events_Mean));
    sum_tbl.evperdend_mean(iC) = mean(tmpdata.Events_per_Dend_Area);
    sum_tbl.evperdend_ci_hspan(iC) = 1.96*std(tmpdata.Events_per_Dend_Area)/sqrt(numel(tmpdata.Events_per_Dend_Area));
    sum_tbl.dend_area_mean(iC) = mean(tmpdata.Dend_Area_um);
    sum_tbl.dend_area_ci_hspan(iC) = 1.96*std(tmpdata.Dend_Area_um)/sqrt(numel(tmpdata.Dend_Area_um));
    sum_tbl.norm_roi_area_mean(iC) = mean(tmpdata.StructNorm_Area);
    sum_tbl.norm_roi_area_ci_hspan(iC) = 1.96*std(tmpdata.StructNorm_Area)/sqrt(numel(tmpdata.StructNorm_Area));
%     sum_tbl.ev_per_norm_roi_area_mean(iC) = mean(tmpdata.ToT_Num_Events*100/tmpdata.StructNorm_Area);
%     sum_tbl.norm_roi_area_ci_hspan(iC) = 1.96*std(tmpdata.StructNorm_Area)/sqrt(numel(tmpdata.StructNorm_Area));
    

end

% Calculate Rec mean
reclist = unique(data2pm.Rec);
nrecs = numel(reclist);

rec_tbl = table();
rec_tbl.Rec = NaN(nrecs,1);
rec_tbl.Dist_mean = NaN(nrecs,1);
rec_tbl.PeaR_mean = NaN(nrecs,1);

for iRec = 1:nrecs
    pter = strcmp(data2pm.Rec, reclist{iRec});
    rec_tbl.Rec(iRec) = iRec;
    rec_tbl.Dist_mean(iRec) = mean(data2pm.Dist_um(pter));
    rec_tbl.PeaR_mean(iRec) = mean(data2pm.PearsonR(pter));
end


%% Average Pearson R per distance
% tmpdata = data2pm(strcmp(data2pm.Sensor, 'A184V') & data2pm.PearsonR >= 0.01,:);
tmpdata = data2pm(strcmp(data2pm.Sensor, 'A184V'),:);

% distbins = 0:5:50; % max(tmpdata.Dist_um); 
distbins = 0:0.1:1; % max(tmpdata.PearsonR);
ndbins = numel(distbins)-1;
result_mtrx = table();

anova_val = [];
anova_grp = [];

for iD = 1:ndbins
%     mask = tmpdata.Dist_um >= distbins(iD) & tmpdata.Dist_um < distbins(iD+1);
    mask = tmpdata.PearsonR >= distbins(iD) & tmpdata.PearsonR < distbins(iD+1);
    result_mtrx.dist_bin(iD) = distbins(iD+1);
    result_mtrx.Cnt(iD) = sum(mask);
    result_mtrx.R_av(iD) = mean(tmpdata.PearsonR(mask));
    result_mtrx.R_CI(iD) = 1.96*std(tmpdata.PearsonR(mask))/sqrt(numel(tmpdata.PearsonR(mask)));
    result_mtrx.EvRate_av(iD) = mean(tmpdata.MeanEvRate(mask));
    result_mtrx.EvRate_CI(iD) = 1.96*std(tmpdata.MeanEvRate(mask))/sqrt(numel(tmpdata.MeanEvRate(mask)));
%     result_mtrx.Cnt_above60(iD) = numel(tmpdata.PearsonR(tmpdata.Dist_um >= distbins(iD) &...
%         tmpdata.Dist_um < distbins(iD+1) & tmpdata.PearsonR > 0.6));
    anova_val = [anova_val; tmpdata.MeanEvRate(mask)];
    anova_grp = [anova_grp; repelem(iD, sum(mask))'];
end

% dist_av_mtrx(:,2) = dist_av_mtrx(:,2)./max(dist_av_mtrx(:,2));
