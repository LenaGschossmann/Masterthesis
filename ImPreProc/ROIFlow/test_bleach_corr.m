
figure();
spcnt=1;
tic;
% traces = roidata.traces(1:300,:);
fitcnt=0;
for itr=1:25
% Y = traces(itr,:);
traces = roidata.traces(itr,:);
X = 1:size(traces,2);
[newtrc,smtrc] = bleach_correction(traces, size(traces,2), round(1/roidata.frametime_s));
%
%     modelfun = @(b,x)b(1) + b(2).*(x(:).^b(3));
%     beta0 = [0 0 -1];
%     mdl = fitnlm(X,Y,modelfun,beta0);
%     coeff = table2array(mdl.Coefficients(:,1));
%     smtrc = coeff(1)+coeff(2).*(X.^coeff(3));
%
%     newtrc = Y./(smtrc./mean(smtrc));
%     newtrc = newtrc+mean(Y(1:5))-mean(newtrc);

%     newtrc = Y-smtrc;
%     newtrc = newtrc+mean(Y(1:5),2);

    subplot(5,5,spcnt);
    plot(traces(1:600));
    hold on, plot(smtrc(1:600),'r', 'LineWidth',1.5);
%     hold on, plot(newtrc(1:600),'g');
    spcnt=spcnt+1;
end

toc

