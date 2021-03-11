%% Display - part 2 of ROI detection
backgroundim = reshape(background,[imsize(1) imsize(2)]);
% Create rgb image
rim = zeros(imsize(1), imsize(2),1);
gim = zeros(imsize(1), imsize(2),1);
bim = zeros(imsize(1), imsize(2),1);
bg_col = [0.1529 0.1882 0.6706]; % RGB
% noise_col = [0.0980 0.4118 0.2627]; % RGB
% rp_mpi = rescale(max(dFoFim,[],3)).*300;
% rp_mpi = rescale(aip).*2;

% Mark background: dark blue
rim(backgroundim==1) = bg_col(1);
gim(backgroundim==1) = bg_col(2);
bim(backgroundim==1) = bg_col(3);

rgbim = cat(3, rim, gim,bim);
% colormap(cmap');

rnd = numel(n_rois);

% Superimpose combined ROIs - all in one image
figure(); imshow(rgbim);
tmpnrois = n_rois(rnd);
bounds = bounds_cell; %bounds_cell_comb;
for cc = 1:tmpnrois
    hold on;
    B = bounds{cc,1};
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'Color',bounds{cc,2}', 'LineWidth', 1)
    end
end

% Combined ROIs - one per subplot
figure();
tmpnrois = n_rois(rnd);
csp = 4; rsp = 3; cntsp=0;
for cc = 1:tmpnrois
    cntsp = cntsp+1;
    subplot(rsp,csp,cntsp),imshow(rgbim);
    hold on;
    B = bounds_cell_comb{cc,1};
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'Color',bounds_cell_comb{cc,2}', 'LineWidth', 1)
    end
    title(sprintf('Combined ROI %i', cc));
    if cntsp == csp*rsp, figure(); cntsp=0; end
end

% Individual ROIs - one per subplot
comb_pter = [bounds_cell{:,3}]';
tmpnrois = n_rois(1);
csp = 4; rsp = 3; cntsp=0;
figure();
for cc = 1:tmpnrois
    cntsp=cntsp+1;
    subplot(rsp,csp,cntsp), imshow(rgbim);
    hold on;
    B = bounds_cell{cc,1};
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'Color',bounds_cell{cc,2}', 'LineWidth', 1)
    end
    title(sprintf('ROI %i, (comb. %i)', cc, comb_pter(cc)));
    if cntsp == csp*rsp, figure(); cntsp=0; end
end

% Individual ROIs - grouped by combination
comb_pter = [bounds_cell{:,3}]';
tmpnrois = n_rois(rnd);
figure();
csp = 4; rsp = 3;cntsp=0;
for cc_comb = 1:tmpnrois
    cntsp=cntsp+1;
    cc_ids = find(comb_pter == cc_comb);
    subplot(rsp,csp,cntsp), imshow(rgbim);
    hold on;
    for cc = 1:numel(cc_ids)
        B = bounds_cell{cc_ids(cc),1};
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'Color',bounds_cell{cc_ids(cc),2}', 'LineWidth', 1);
        end
    end
    title(sprintf('Combined ROI %i (# %i)', cc_comb, numel(cc_ids)));
    hold off;
    if cntsp == csp*rsp, figure(); cntsp=0; end
end

%% Sanity check
spr=16; spc=13;
posccim = [1/spc 13/spr 2/spc 2/spr;...
    1/spc 10/spr 2/spc 2/spr;...
    1/spc 7/spr 2/spc 2/spr;...
    1/spc 4/spr 2/spc 2/spr;...
    1/spc 1/spr 2/spc 2/spr];

poscctr = [4/spc 13/spr 8/spc 2/spr;...
    4/spc 10/spr 8/spc 2/spr;...
    4/spc 7/spr 8/spc 2/spr;...
    4/spc 4/spr 8/spc 2/spr;...
    4/spc 1/spr 8/spc 2/spr];
% figure();
% spcnter = 1;
% for cc = 1:num_connc2
%     cc_px = roi_idx_mtrx_filled(roi_idx_mtrx_filled(:,1) == cc,2:4);
%     traces = NaN(size(cc_px,1),imsize(3));
%     for ipx = 1:size(cc_px,1), traces(ipx,:) = reshape(dFoFim(cc_px(ipx,1), cc_px(ipx,2),:), [1 imsize(3)]); end
%     avtrace = mean(traces,1);
%     subplot(5,1,spcnter);
%     plot(avtrace); title(sprintf('Comp: %i',cc));
%     if spcnter == 5, spcnter = 1; figure(); 
%     else, spcnter = spcnter+1;
%     end
% end

%% Plot single ROIs wiht traces
% spcnter = 1;
% colormap(cmap');
% num_cc = n_rois(rnd);
% idx_mtrx = roi_idx_cell{1,1};
% comb_pter = roi_idx_cell{1,2};
% for comb_cc = 1:num_cc
%     tmp_idx_mtrx = idx_mtrx(comb_pter == comb_cc,:);
%     tmp_cc = unique(tmp_idx_mtx(:,1));
%     figure();
%     for cc = 1:tmp_cc
%         B = bounds_cell{cc,1};
%         cc_px = tmp_idx_mtrx(tmp_idx_mtrx(:,1) == cc,2:4);
%         t = cc_px(1,3);
%         traces = NaN(size(cc_px,1),imsize(3));
%         for ipx = 1:size(cc_px,1), traces(ipx,:) = reshape(dFoFim(cc_px(ipx,1), cc_px(ipx,2),:), [1 imsize(3)]); end
%         avtrace = mean(traces,1);
%         
%         subplot('Position', posccim(spcnter,:));
%         imshow(uint8(aip));
%         hold on
%         for k = 1:length(B)
%             boundary = B{k};
%             plot(boundary(:,2), boundary(:,1), 'Color',bounds_cell{cc,2}', 'LineWidth', 1)
%         end
%         title(sprintf('ROI: %i, t=%i',cc,t));
%         
%         subplot('Position', poscctr(spcnter,:));
%         plot(avtrace); title(sprintf('Comp: %i',cc));
%         
%         if spcnter == 5, spcnter = 1; figure();
%         else, spcnter = spcnter+1;
%         end
%     end
% end

%% Plot combined ROIs with traces
figure();
spcnter = 1;
colormap(cmap');
num_cc = n_rois(rnd);
idx_mtrx = roi_idx_cell{rnd,1};
for cc = 1:30 %num_cc
    B = bounds_cell_comb{cc,1};
    cc_px = idx_mtrx(idx_mtrx(:,1) == cc,2:3);
    traces = NaN(size(cc_px,1),imsize(3));
    for ipx = 1:size(cc_px,1), traces(ipx,:) = reshape(dFoFim(cc_px(ipx,1), cc_px(ipx,2),:), [1 imsize(3)]); end
    avtrace = mean(traces,1);
    
    subplot('Position', posccim(spcnter,:));
    imagesc(uint8(uint8_aip));colormap('gray');axis('off');daspect([1 1 1]);
    hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'Color',bounds_cell_comb{cc,2}', 'LineWidth', 1)
    end
    title(sprintf('ROI: %i',cc));
    
    subplot('Position', poscctr(spcnter,:));
    plot(avtrace); title(sprintf('Comp: %i',cc));
    
    if spcnter == 5, spcnter = 1; figure(); 
    else, spcnter = spcnter+1;
    end
end

%% Plot overlap matrix
rnd = 4;
figure();
imagesc(overlap_cell{rnd,1});
set(gca, 'XTick', 1:n_rois(rnd)); % center x-axis ticks on bins
set(gca, 'YTick', 1:n_rois(rnd)); % center y-axis ticks on bins
colormap('summer');

figure();
imagesc(overlap_cell{rnd,2});
set(gca, 'XTick', 1:n_rois(rnd)); % center x-axis ticks on bins
set(gca, 'YTick', 1:n_rois(rnd)); % center y-axis ticks on bins
colormap('gray');


%% Testing threshold for responding px classificatino
% px1 = [15,34]; px2 =[21,31]; px3 = [19,21]; px4 = [16,30]; px5= [16,27]; %set1
% px1 = [16,44]; px2 =[13,45]; px3 = [14,47]; px4 = [16,50]; px5= [10,51]; %set2
% % id=find(thresh_ratio > 1.5);% & thresh_ratio > 1.1);
% % [idr,idc]=ind2sub(imsize, id(1:5));
% % px1 = [idr(1),idc(1)]; px2 =[idr(2),idc(2)]; px3 = [idr(3),idc(3)]; px4 = [idr(4),idc(4)]; px5= [idr(5),idc(5)];
% figure();
% subplot(5,1,1); plot(reshape(dFoFim(px1(1),px1(2),:),[imsize(3) 1])); hold on, yline(threshold(px1(1),px1(2)), 'b');
% tmprp = find(reshape(response(px1(1),px1(2),:),[1 imsize(3) 1]) == 1); for ii=1:numel(tmprp), xline(tmprp(ii),'r'); end, title(sprintf('Px 1, thresh: %i', threshold(px1(1),px1(2))));
% subplot(5,1,2); plot(reshape(dFoFim(px2(1),px2(2),:),[imsize(3) 1])); hold on, yline(threshold(px2(1),px2(2)), 'b');
% tmprp = find(reshape(response(px2(1),px2(2),:),[1 imsize(3) 1]) == 1); for ii=1:numel(tmprp), xline(tmprp(ii),'r'); end, title(sprintf('Px 2, thresh: %i', threshold(px2(1),px2(2))));
% subplot(5,1,3); plot(reshape(dFoFim(px3(1),px3(2),:),[imsize(3) 1])); hold on, yline(threshold(px3(1),px3(2)), 'b');
% tmprp = find(reshape(response(px3(1),px3(2),:),[1 imsize(3) 1]) == 1); for ii=1:numel(tmprp), xline(tmprp(ii),'r'); end, title(sprintf('Px 3, thresh: %i', threshold(px3(1),px3(2))));
% subplot(5,1,4); plot(reshape(dFoFim(px4(1),px4(2),:),[imsize(3) 1])); hold on, yline(threshold(px4(1),px4(2)), 'b');
% tmprp = find(reshape(response(px4(1),px4(2),:),[1 imsize(3) 1]) == 1); for ii=1:numel(tmprp), xline(tmprp(ii),'r'); end, title(sprintf('Px 4, thresh: %i', threshold(px4(1),px4(2))));
% subplot(5,1,5); plot(reshape(dFoFim(px5(1),px5(2),:),[imsize(3) 1])); hold on, yline(threshold(px5(1),px5(2)), 'b');
% tmprp = find(reshape(response(px5(1),px5(2),:),[1 imsize(3) 1]) == 1); for ii=1:numel(tmprp), xline(tmprp(ii),'r'); end, title(sprintf('Px 5, thresh: %i', threshold(px5(1),px5(2))));
