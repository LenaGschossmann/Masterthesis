tmp = zeros(1200,1);
bin = zeros(1200,1);
xlidx = find(bin==1);
figure,plot(tmp, 'k','Linewidth',1.3);
title('Params3 - Roi5');
hold on, for iX = 1:numel(xlidx), xline(xlidx(iX),'r','Linewidth', 1.3); end
hold on, yline(0.2170,'k');


trc = zeros(33,1);

roi = 3; ev = 8;
figure, plot(trc, 'Linewidth',1.5);
hold on, plot(1:numel(trc), trc,'ko');
hold on, yline(eventinfo.eventdata{roi,10}(2),'k');
hold on, yline(eventinfo.eventdata{roi,10}(3),'k');
title(sprintf('Roi: %i, Event: %i', roi, ev));


xlidx = find(eventinfo.eventdata{90,1}==1);
hold on, for iX = 1:numel(xlidx), xline(xlidx(iX),'r','Linewidth', 1.3); end


figure, plot(tmpdFoF);
hold on, for iX = 1:numel(save_idx), xline(save_idx(iX),'r','Linewidth', 1.3); end
hold on, yline(mean_neg_amp*2,'k');



figure(); subplot(3,1,2), plot(a(:,1));subplot(3,1,1), plot(a(:,2));subplot(3,1,3), plot(a(:,3));
trc1=a(:,1); trc2=a(:,2);
tmpmask = trc1 > 1.5*std(trc1) & trc2 > 1.5*std(trc2);
trc_corrA = sum(trc1(tmpmask) .* trc2(tmpmask));
av1 = mean(trc1); av2 = mean(trc2);
r = sum((av1-trc1).*(av2-trc2)) / sqrt(sum((av1-trc1).^2)*sum((av2-trc2).^2));

trc1=a(:,1); trc2=a(:,3);
tmpmask = trc1 > 1.5*std(trc1) & trc2 > 1.5*std(trc2);
trc_corrA = [trc_corrA sum(trc1(tmpmask) .* trc2(tmpmask))];
av1 = mean(trc1); av2 = mean(trc2);
r = [r sum((av1-trc1).*(av2-trc2)) / sqrt(sum((av1-trc1).^2)*sum((av2-trc2).^2))];






