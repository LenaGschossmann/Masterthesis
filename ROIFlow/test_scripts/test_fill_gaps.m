%% Test fill-holes approach
fill_thresh = 0.2;
cond_result_cell = cell(5,1);
dim_stack = size(response);
p = 'C:\Users\lena_\Projects\iGluSnFR\Analysis\Tests\LenaMatlabCode\DD_ROI_detection\fill_holes\compare_cond';

for cc = 1:num_connc2
    %% Define rect
    cc_px = roi_idx_mtrx(roi_idx_mtrx(:,1) == cc,2:4);
    rlimit = [min(cc_px(:,1))-1 max(cc_px(:,1))+1];
    climit = [min(cc_px(:,2))-1 max(cc_px(:,2))+1];
    if rlimit(1) < 1, rlimit(1) = 1; end
    if rlimit(2) > dim_stack(1), rlimit(2) = dim_stack(1); end
    if climit(1) < 1, climit(1) = 1; end
    if climit(2) > dim_stack(2), climit(2) = dim_stack(2); end
    t = cc_px(1,3);
    
    %% Cond 1 - standard smoothing 3x1 window, horizontally&vertically
    tmpav1 = double(response(:,:,t));
    tmpav2 = double(response(:,:,t));
    % Average horizontally
    for r = rlimit(1)+1:rlimit(2)-1
        for pos = climit(1)+1:climit(2)-2
            tmpav1(r,pos) = mean(double(response(r,pos-1:pos+1)));
        end
    end
    % ...and vertically
    for c = climit(1)+1:climit(2)-1
        for pos = rlimit(1)+1:rlimit(2)-2
            tmpav2(pos, c) = mean(double(response(pos-1:pos+1,c)));
        end
    end
    tmpav = tmpav1 + tmpav2;
%     tmpav(tmpav >= fill_thresh) = 1;
    cond_result_cell{1,1} = tmpav;
    
    %% Cond 2 - standard smoothing 3x3 window
    tmpav = double(response(:,:,t));
    for pos_r = rlimit(1)+1:rlimit(2)-2
        for pos_c = climit(1)+1:climit(2)-2
            tmpav(pos_r,pos_c) = mean(double(response(pos_r-1:pos_r+1,pos_c-1:pos_c+1)),  'all');
        end
    end
%     tmpav(tmpav >= fill_thresh) = 1;
    cond_result_cell{2,1} = tmpav;
    
    %% Cond 3 - Averaging á la Lena, 3x1 window, starting
    % horizontally&vertically with original
    tmpav1 = double(response(:,:,t));
    tmpav2 = double(response(:,:,t));
    % Average horizontally
    for r = rlimit(1):rlimit(2)
        for pos = climit(1)+1:climit(2)-2
            tmpav1(r,pos) = mean(tmpav1(r,pos-1:pos+1));
        end
    end
    % ...and vertically
    for c = climit(1):climit(2)
        for pos = rlimit(1)+1:rlimit(2)-2
            tmpav2(pos, c) = mean(tmpav2(pos-1:pos+1,c));
        end
    end
    tmpav = tmpav1+tmpav2;
%     tmpav(tmpav >= 0.5) = 1;
    cond_result_cell{3,1} = tmpav;
    
    
    %% Cond 4 - Averaging á la Lena, 3x1 window, start upper left
    tmpav = double(response(:,:,t));
    % Average horizontally
    for r = rlimit(1):rlimit(2)
        for pos = climit(1)+1:climit(2)-2
            tmpav(r,pos) = mean(tmpav(r,pos-1:pos+1));
        end
    end
    % ...and vertically
    for c = climit(1):climit(2)
        for pos = rlimit(1)+1:rlimit(2)-2
            tmpav(pos, c) = mean(tmpav(pos-1:pos+1,c));
        end
    end
%     tmpav = tmpav3 >= fill_thresh;
    cond_result_cell{4,1} = tmpav;
    
    
    %% Cond 5 - Averaging á la Lena, 3x1 window, start lower right
    tmpav = double(response(:,:,t));
    % Average horizontally
    r=rlimit(2); pos=climit(2)-2;
    while r >= rlimit(1)
        while pos >= climit(1)+1
            tmpav(r,pos) = mean(tmpav(r,pos-1:pos+1));
            pos=pos-1;
        end
        r=r-1;
    end
    % ...and vertically
    c=climit(2); pos=rlimit(2)-2;
    while c >= climit(1)
        while pos >= rlimit(1)+1
            tmpav(pos, c) = mean(tmpav(pos-1:pos+1,c));
            pos=pos-1;
        end
        c=c-1;
    end
%     tmpav = tmpav3 >= fill_thresh;
    cond_result_cell{5,1} = tmpav;
    
    %% Plot
    f1=figure('Position', [30 30 1200 700]);
    rsp=3;csp=2;
    subplot(rsp,csp,1);
    colormap('gray');
    imagesc(response(:,:,t).*255);
    axis('off'); title(sprintf('ROI %i', cc));
    for iC = 1:size(cond_result_cell,1)
        subplot(rsp,csp,iC+1);
        imagesc(cond_result_cell{iC,1});
        colorbar(); axis('off');
        title(sprintf('Cond %i', iC));
    end
    saveas(f1, fullfile(p,sprintf('overview_roi_%i.png',cc))); close gcf;
    
    %% Calculate some parameters
    cond_thresh_cell = cell(size(cond_result_cell,1),1);
    mask = false(size(tmpav)); mask(rlimit(1)+1:rlimit(2)-1, climit(1)+1:climit(2)-1) = true;
    tmpav = double(response(:,:,t));
    area = zeros(size(cond_result_cell,1)+1,1);
    area(1,1) = diff(rlimit)*diff(climit);
    area(2,1) = sum(tmpav(rlimit(1)+1:rlimit(2)-1,climit(1)+1:climit(2)-1),'all');
    center = zeros(size(cond_result_cell,1)+1,2);
    linidx = find(tmpav == 1); [rrp, crp] = ind2sub([dim_stack(1) dim_stack(2)], linidx);
    center(1,1) = mean(crp); center(1,2) = mean(rrp);
    for iC = 1:size(cond_result_cell,1)
        cond_thresh_cell{iC,1} = cond_result_cell{iC,1} >= fill_thresh;
        cond_thresh_cell{iC,1}(~mask) = 0;
        % Area
        area(iC+2,1) = sum(cond_thresh_cell{iC,1}(rlimit(1)+1:rlimit(2)-1,climit(1)+1:climit(2)-1),'all');
        % Center
        linidx = find(cond_thresh_cell{iC,1} == 1); [rrp, crp] = ind2sub([dim_stack(1) dim_stack(2)], linidx);
        center(iC+1,1) = mean(crp); center(iC+1,2) = mean(rrp); 
    end
    
    f2=figure('Position', [30 30 1200 700]);
    labels = {'Original', 'Cond 1', 'Cond 2', 'Cond 3', 'Cond 4', 'Cond 5'};
    spec = {'+r','+y','+g','+b','+m','+w'};
    subplot(1,2,1);
    for iC = 1:size(center,1), hold on, plot(center(iC,1), center(iC,2),spec{iC});end
    xlim([1 dim_stack(2)]); ylim([1 dim_stack(1)]);
    set(gca,'Color',[0.7 0.7 0.7]); daspect([1 1 1]);
    legend(labels);
    title(sprintf('Center of responding px in ROI %i', cc));
    subplot(1,2,2);
    bp= bar(area); set(gca,'XTickLabel',['Encl.Rect' labels]); xtickangle(30); 
    title('ROI area for different conditions');
    saveas(f2, fullfile(p,sprintf('params_roi_%i.png',cc))); close gcf;
end


