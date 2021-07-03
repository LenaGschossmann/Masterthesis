figure()
ii = 2;
tmpmtrx = log(mtrx); %log(mtrx)
while ii <= size(mtrx,1)
    hold on; plot(mtrx(1,:), tmpmtrx(ii,:), 'Color', [0.8510 0.3255 0.0980], 'LineWidth',1.5); ii=ii+1;
    hold on; plot(mtrx(1,:), tmpmtrx(ii,:), 'Color', [1 0 0], 'LineWidth',1.5); ii=ii+1;
    hold on; plot(mtrx(1,:), tmpmtrx(ii,:), 'Color', [0.1765 0.2471 0.5294], 'LineWidth',1.5); ii=ii+1;
    hold on; plot(mtrx(1,:), tmpmtrx(ii,:), 'Color', [0.3922 0.8314 0.0745], 'LineWidth',1.5); ii=ii+1;
end
xlim([0 30]);ylabel('Counts log-scale'); % Counts log-scale
% ylim([0 100]);
legend('1.3', '1.5', '1.75', '1.85');
title('2PM CTRL files');


t=zeros(imsize(1), imsize(2));
testlist_cc = unique(tmp_idx_mtrx(:,1));
testn_cc = numel(testlist_cc);
for iii = 1:testn_cc
    coords = tmp_idx_mtrx(tmp_idx_mtrx(:,1)==testlist_cc(iii),2:3);
    for ii = 1:size(coords,1), t(coords(ii,1),coords(ii,2))= 1; end
end
figure, imshow(t);

%% Plot Hist of F/F values
tmptitle = split(filename,'_'); plottitle = [];
for ii = 1:numel(tmptitle), plottitle = strcat(plottitle,'-', tmptitle{ii}); end
plottitle = plottitle(2:end);
figure();
subplot(2,3,[1 4]), histogram(FoFtraces); ylabel('Counts'); xlabel('F/F values');title('Log-transformed F/F values');
subplot(2,3,[2 5]), histogram(std_noise(~exclude_px)); ylabel('Counts'); xlabel('Noise SD (bg px excluded)'); title('Noise SD');

subplot(2,3,3)
plotdata=FoFtraces./std_noise;
e = linspace(round(min(plotdata,[],'all')), round(max(plotdata,[],'all')), 50);
[n1,e1] = histcounts(plotdata(FoFtraces <= resp_thresh),e); n1 = log(n1); n1(n1<0)=0;
[n2,e2] = histcounts(plotdata(FoFtraces > resp_thresh),e);n2 = log(n2);n2(n2<0)=0;
histogram('BinEdges', e, 'BinCounts', n1, 'FaceColor', 'green');
hold on, histogram('BinEdges', e2, 'BinCounts', n2, 'FaceColor', 'red');
ylabel('Log of Counts'); xlabel('Fold SD'); title('Low bar');legend('non act.', 'act.');
% max(FoFtraces(FoFtraces > resp_thresh), [],'all')

subplot(2,3,6)
[n1,e1] = histcounts(plotdata(FoFtraces <= resp_thresh2),e); n1 = log(n1); n1(n1<0)=0;
[n2,e2] = histcounts(plotdata(FoFtraces > resp_thresh2),e);n2 = log(n2);n2(n2<0)=0;
histogram('BinEdges', e, 'BinCounts', n1, 'FaceColor', 'green');
hold on, histogram('BinEdges', e2, 'BinCounts', n2, 'FaceColor', 'red');
ylabel('Log of Counts'); xlabel('Fold SD'); title('High bar');
sgtitle(plottitle);
plottitle


%% Plottitle
tmptitle = split(filename,'_'); plottitle = [];
for ii = 1:numel(tmptitle), plottitle = strcat(plottitle,'-', tmptitle{ii}); end
plottitle = plottitle(2:end);

%% SDnoise
figure();
histogram(std_noise(~exclude_px,:)); xlabel('SDnoise'); ylabel('Count');title(plottitle);
std(std_noise(~exclude_px,:))

%% SD of FoF > 1, mirrored
test = FoFtraces-1;
test(exclude_px,:) = 0;
testfof = reshape(test,[imsize(1) imsize(2) n_frames]);
pos = testfof(testfof>0);
mirr = [-pos; pos];mirr = mirr+1;
figure();
histogram(mirr); xlabel('F/F'); ylabel('Count');title(plottitle); xlim([-9 10]);
std(mirr)

prctile(pos, 70)
prctile(pos, 80)
prctile(pos, 90)

%% FoF as SDnoise-multiples
test = FoFtraces-1;
test(exclude_px,:) = 0;
testfof = reshape(test,[imsize(1) imsize(2) n_frames]);
testsd = reshape(std_noise,[imsize(1) imsize(2)]);
testfof = (testfof)./testsd;
pos = testfof(testfof>0);
mirr = [-pos; pos];

figure();
histogram(mirr); xlabel('SD-multiples'); ylabel('Count');title(plottitle); xlim([-10 10]);
std(mirr)

prctile(pos, 70)
prctile(pos, 80)
prctile(pos, 90)


%% Component-based
figure();
histogram(cc_noise,10);
xlabel('SD of comp av. F/F<1'); ylabel('Count'); title(plottitle);

figure();
tmp = (cc_fof-1)./cc_noise;
for ii = 1:cc_n_tot
    plot(1+tmp(ii,:)); hold on;
end
xlabel('frames'); ylabel('F/F expressed as SD multiple'); title(strcat('Comp. trc - ',plottitle));

figure();
histogram(cc_pooled_pos);
xlabel('F/F expressed as SD multiple'); ylabel('Counts'); title(strcat('Comp. pooled - ',plottitle));



% ***************************
% theta1 = 2;
% theta2 = 2.75;
% cnttot = numel(pos);
% thetas = 0:0.1:4;
% results = NaN(numel(thetas),1);
% for ii=1:numel(thetas)
%     tmptheta = thetas(ii);
%     cnt=sum(pos>=tmptheta);
%     results(ii) = 100*cnt/cnttot;
% end
% figure()
% subplot(1,2,1), histogram(mirr);
% subplot(1,2,2), plot(thetas,results, 'LineWidth',1.5); ylim([0 100]);xlabel('SD multiplier');ylabel('AUC');
% yline(10);yline(20);yline(30);yline(40);
% sgtitle(plottitle);
% std(mirr)
% [n,e]=histcounts(mirr,100);
% hm = max(n)/2;
% r1=find(n>hm,1,'first');
% r2=find(n>hm,1,'last');
% fwhm = e(r2)-e(r1)
% thetas(find(results<=20,1,'first'))
% thetas(find(results<=15,1,'first'))
% thetas(find(results<=10,1,'first'))
% thetas(find(results<=5,1,'first'))
% thetas(find(results<=2.5,1,'first'))
% 
% % figure, plot(diff(results))
% sl(1) = 0.025/(thetas(find(results<=5,1,'first'))-thetas(find(results<=2.5,1,'first')));
% sl(2) = 0.075/(thetas(find(results<=10,1,'first'))-thetas(find(results<=2.5,1,'first')));
% sl(3) = 0.10/(thetas(find(results<=15,1,'first'))-thetas(find(results<=5,1,'first')));
% sl(4) = 0.05/(thetas(find(results<=10,1,'first'))-thetas(find(results<=5,1,'first')));
% sl(5) = 0.05/(thetas(find(results<=20,1,'first'))-thetas(find(results<=15,1,'first')));
% sl(6) = 0.1/(thetas(find(results<=20,1,'first'))-thetas(find(results<=10,1,'first')));
% sl(7) = 0.1/(thetas(find(results<=30,1,'first'))-thetas(find(results<=20,1,'first')));
% sl(8) = 0.2/(thetas(find(results<=30,1,'first'))-thetas(find(results<=10,1,'first')));
% sl(9) = 0.2/(thetas(find(results<=40,1,'first'))-thetas(find(results<=20,1,'first')));
% sl(10) = 0.3/(thetas(find(results<=40,1,'first'))-thetas(find(results<=10,1,'first')));
% sl(11) = 0.3/(thetas(find(results<=80,1,'first'))-thetas(find(results<=50,1,'first')));

      

%% ******* Event detection
figure, plot(traces(tmproi,1:500)), hold on, plot(trc_filtrd(tmproi,1:500)), legend()
% figure, plot(traces(tmproi,1:500)), hold on, plot(trc_filtrd_2(tmproi,1:500)), legend()

figure, plot(tmpFoF), hold on, yline(1+tmpSDnoise); yline(1-tmpSDnoise);
yline(1+tmpSDnoise*SD_factor, 'r'); yline(1+tmpSDnoise*SD_factor_safe, 'g')

figure, plot(tmpFoF), hold on;
% l = find(FoF_crossing==1);
for ii = 1:numel(onset_idx)
    xline(onset_idx(ii));
end
l = find(FoF1d_crossing==1);
for ii = 1:numel(l)
    xline(l(ii),'r');
end
hold on, yline(1+abs(mean(tmpFoF(ctrl_peak_idx))));

figure, plot(ctrl_FoF), hold on;
for ii = 1:numel(ctrl_onset_idx)
    xline(ctrl_onset_idx(ii));
end


%% Filter comparison
% Filter traces
fk1 = [1 1 1];
fk2 = [1 2 1];
fk3 = [1 4 6 4 1];

trc_filtrd1 = smooth_data(roidata.traces, fk1);
trc_filtrd2 = smooth_data(roidata.traces, fk2);
trc_filtrd3 = smooth_data(roidata.traces, fk3);

% Calculate FoF
[~,FoF, ~] = rollBase_dFoF(roidata.traces,roidata.baseline_frames,PARAMS.dFoF_baseline_shift, 'roll');
[~,FoF_filtrd1, ~] = rollBase_dFoF(trc_filtrd1,roidata.baseline_frames,PARAMS.dFoF_baseline_shift, 'roll');
[~,FoF_filtrd2, ~] = rollBase_dFoF(trc_filtrd2,roidata.baseline_frames,PARAMS.dFoF_baseline_shift, 'roll');
[~,FoF_filtrd3, ~] = rollBase_dFoF(trc_filtrd3,roidata.baseline_frames,PARAMS.dFoF_baseline_shift, 'roll');

idx = 10;
yr = [min([FoF(3,:) FoF_filtrd1(3,:) FoF_filtrd2(3,:) FoF_filtrd3(3,:)]) max([FoF(3,:) FoF_filtrd1(3,:) FoF_filtrd2(3,:) FoF_filtrd3(3,:)])]; 
range = 1:size(FoF_filtrd1,2);
xl = 1:numel(range);
figure();
subplot(4,1,1), plot(xl,FoF(idx,range)), title('No Filter'), ylim(yr);
subplot(4,1,2), plot(xl,FoF_filtrd1(idx,range)), title('[1 1 1]'), ylim(yr);
subplot(4,1,3), plot(xl,FoF_filtrd2(idx,range)), title('[1 2 1]'), ylim(yr);
subplot(4,1,4), plot(xl,FoF_filtrd3(idx,range)), title('[1 4 6 4 1]'), ylim(yr);



