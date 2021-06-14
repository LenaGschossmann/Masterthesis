%% Extract offset correction values

% Add paths
wd = pwd(); p = split(wd,'\'); p = p(1:end-1); pp = []; for ii=1:numel(p), pp = fullfile(pp,p{ii}); end
addpath(genpath(pp)), clear('p','pp');

bittype = 12;

selpath = uigetdir('D:\Masterthesis stuff\2PM');
dirlist = dir(selpath); dirlist = dirlist(3:end);

for iDir = 1:size(dirlist,1)
    tmpdir = strcat(selpath,'\',dirlist(iDir).name);
    reclist = dir(tmpdir);
    reclist = reclist(3:end);
    ovw_rec = NaN(size(reclist,1),4);
    % col1: offset, col2: average value, col3: number 0-px, col4:
    % number 4096-px % col5: number total px
    for iRec = 1:size(reclist,1)
        tmpname = reclist(iRec).name;
        tmpnamesplit = split(tmpname,'.'); tmpnamesplit = tmpnamesplit{1};
        tmpnamesplit = split(tmpnamesplit,'_'); offset = tmpnamesplit{end};
        
        ovw_rec(iRec,1) = str2num(offset);
        
        % Read file
        im = read_timeseries(strcat(tmpdir,'\', tmpname), 'load');
        
        npx = numel(im);
        ovw_rec(iRec,3) = sum(im == 0, 'all');
        ovw_rec(iRec,4) = sum(im == 2^bittype-1, 'all');
        ovw_rec(iRec,5) = npx;
        ovw_rec(iRec,2) = mean(im(im > 0 & im < 2^bittype), 'all');
        % Do not regard values of 0 or maximum for bittype
    end
    [~,sortidx] = sort(ovw_rec(:,1), 'descend');
    ovw_rec = ovw_rec(sortidx,:);
    
    % Save
    settings = split(dirlist(iDir).name,'_');
    pmt = settings{1}(1:end-3); px = settings{2}(1:end-2);
   
    % Sanity check
    f = figure();
    plot(ovw_rec(:,1), ovw_rec(:,2), 'LineStyle','-', 'Color', [0,0,0], 'LineWidth',1.5)
    hold on, plot(ovw_rec(:,1), ovw_rec(:,2), '*r');
    title(sprintf('Photon count with PMT: %s and resolution: %s px', pmt, px));
    xlabel('offset'); ylabel('Photon count');
    
    tmpout = strcat(tmpdir,'\photon_count.png');
    saveas(f, tmpout, 'png');
    close(f);
    
    % Write total photon count average to spreadsheet
    tmptbl = table();
    tmptbl.offset = ovw_rec(:,1);
    tmptbl.av_photon_count = ovw_rec(:,2);
    tmptbl.count_zero_px = ovw_rec(:,3);
    tmptbl.percent_zero = ovw_rec(:,3)*100./ovw_rec(:,5);
    tmptbl.count_bitmax_px = ovw_rec(:,4);
    tmptbl.percent_bitmax = ovw_rec(:,4)*100./ovw_rec(:,5);
    tmptbl.count_total_px = ovw_rec(:,5);
    tmpout = strcat(tmpdir,'\photon_count_overview.xlsx');
    writetable(tmptbl, tmpout);
    clear('tmptbl');
    
    zeroval = ovw_rec(1,2); % After sorting 0 offset is in first line
    for iRec = 2:size(ovw_rec,1)
        tmpval = NaN(1,3);
        tmpout = strcat(tmpdir,'\diff_PMT', pmt, '_px', px, '_offset_', num2str(ovw_rec(iRec,1)), '.csv');
        tmpval(1,1) = str2num(pmt);
        tmpval(1,2) = ovw_rec(iRec,1);
        tmpval(1,3) = zeroval - ovw_rec(iRec,2);
        writematrix(tmpval, tmpout);
    end
end
