function threshold = find_prctl_threshold(dFoFtraces, background, perc_thresh_range,...
    perc_cutoff, perc_fac,connc_px_thresh, imsize, savepath, filename)

disp('Determine the appropriate threshold...');

test_vals = perc_thresh_range(1):perc_thresh_range(2);
detect_roi_list = zeros(numel(test_vals),1);
mask1 = repmat(reshape(background, [imsize(1) imsize(2)]), [1 1 imsize(3)]);

for iTh = 1:numel(test_vals)
    tmp_thresh = test_vals(iTh);
    num_cc_px = 0;
    % Detect responding px
    resp_thresh = prctile(dFoFtraces,tmp_thresh,2);
    response = dFoFtraces > perc_fac.*resp_thresh;
    response = reshape(response,[imsize(1) imsize(2) imsize(3)]);
    response(mask1) = false;
    
    % Component detection in 2D
    % Detect connected pixel (connectivity: 8) in each frame and save as linear index
    % list (2D array indices!)
    all_connc_cell=cell(imsize(3),1);
    tic;
    parfor t =1:imsize(3)
        bw=response(:,:,t);
        connc = bwconncomp(bw);
        all_connc_cell{t,1} = connc.PixelIdxList;
    end
    toc
    tic;
    for t =1:imsize(3)
        conn_px = cellfun(@numel, all_connc_cell{t,1}, 'UniformOutput', false);
        conn_px = cellfun(@(x) x, conn_px);
        num_cc_px = num_cc_px+sum(conn_px(conn_px > connc_px_thresh));
    end
    toc
    detect_roi_list(iTh) = num_cc_px;
end

c_cum = cumsum(detect_roi_list)./max(cumsum(detect_roi_list));
[~,idx] = min(abs(c_cum - perc_cutoff));
threshold = test_vals(idx);

f = figure();
plot(test_vals, c_cum); xlabel('Percentile Threshold'); ylabel('Cummulative sum of connected px');
xline(threshold);
outputname = strcat(savepath,'\prctile_threshold_',filename,'.png');
if isa(outputname,'cell'), outputname= outputname{1}; end
saveas(f, outputname);
close gcf;

end