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




trc1=a; trc2=a2;
tmpmask = trc1 > 1.5*std(trc1) & trc2 > 1.5*std(trc2);
trc_corr = sum(trc1(tmpmask) .* trc2(tmpmask))

