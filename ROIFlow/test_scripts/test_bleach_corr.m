
% 
% rnd = poissrnd(50,[3000 1]);
% e = 20+exp(-0.2.*(1:3000));
% rnd_bleach = e-rnd;
rpx=[210 210 211 211 212 212 138 138 139 139 140 140];
cpx=[58 59 58 59 58 59 77 78 77 78 77 78];
pxlin = sub2ind([236 236], rpx,cpx);

figure();
spcnt=1;
tic;
% Y = roidata.traces(1:300,:);
% Y=dFoFtraces;
% ft = ft_s; %roidata.frametime_s
fitcnt=0;

for i=1:size(pxlin,2)
    itr = pxlin(i);
    % Y = traces(itr,:);
    % traces = Y(itr,:);
    X = 1:size(traces,2);
    % [newtrc,smtrc] = bleach_correction(traces, size(traces,2), round(4/ft_s));
    %
    
    subplot(4,3,spcnt);
    hold on,plot(traces(itr,1:600));
%     hold on, plot(sm_traces(itr,1:600),'r', 'LineWidth',1.5);
    hold on, plot(sm_traces(1:600),'r', 'LineWidth',1.5);
    %     ylim([-100 1000]);
    hold on, plot(shift_corr_traces(itr,1:600),'g');
    title(sprintf('PX (%i,%i)',rpx(i),cpx(i)));
    spcnt=spcnt+1;
end

toc

figure, plot(mean(traces,1)), hold on,plot(mean(shift_corr_traces,1),'g'),...
     plot(sm_traces+(mean(mean(traces,1),2)-mean(sm_traces)), 'r', 'LineWidth',1.5),...
     legend('original','mean', 'fit'), title('Mean of corrected traces')
