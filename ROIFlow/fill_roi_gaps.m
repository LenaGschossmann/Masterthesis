function [roi_idx_mtrx_filled] = fill_roi_gaps(imsize, roi_idx_mtrx, cid, fill_thresh)

% response_filled = zeros(imsize);

roi_idx_mtrx_filled = [];
list_cc = unique(roi_idx_mtrx(:,1));
num_cc = numel(list_cc);

for cc = 1:num_cc
    cc_px = roi_idx_mtrx(roi_idx_mtrx(:,1) == list_cc(cc),cid(1):cid(2));
    if ~isempty(cc_px)
        response = zeros(imsize(1), imsize(2));
        for ipx = 1:size(cc_px,1), response(cc_px(ipx,1), cc_px(ipx,2)) = 1; end
        
        % Define enclosing rectangle
        rlimit = [min(cc_px(:,1))-1 max(cc_px(:,1))+1];
        climit = [min(cc_px(:,2))-1 max(cc_px(:,2))+1];
        if rlimit(1) < 1, rlimit(1) = 1; end
        if rlimit(2) > imsize(1), rlimit(2) = imsize(1); end
        if climit(1) < 1, climit(1) = 1; end
        if climit(2) > imsize(2), climit(2) = imsize(2); end
        
        tmpav = double(response);
        %     tmpav1 = double(response);
        %     tmpav2 = double(response);
        
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
        %     tmpav = tmpav1+tmpav2;
        tmpav = tmpav >= fill_thresh;
        
        mask = false(size(tmpav)); mask(rlimit(1)+1:rlimit(2)-1, climit(1)+1:climit(2)-1) = true;
        tmpav(~mask) = 0;
        
        % Make sure the resulting px are connected
        test_connc = bwconncomp(tmpav);
        if numel(test_connc.PixelIdxList) > 1
            [~,m] = max(cellfun(@numel, test_connc.PixelIdxList));
            new_px = test_connc.PixelIdxList{m};
            tmpav(:) = 0; tmpav(new_px) = 1;
        end
        
        if ~all(~tmpav,'all')
            [r_f, c_f] = ind2sub([imsize(1) imsize(2)], find(tmpav == 1));
            r_f = reshape(r_f, [numel(r_f) 1]);
            c_f = reshape(c_f, [numel(c_f) 1]);
            if diff(cid)>1
                t = cc_px(1,3);
                roi_idx_mtrx_filled = [roi_idx_mtrx_filled; repelem(list_cc(cc),numel(r_f))' r_f c_f repelem(t,numel(r_f))'];
            else
                roi_idx_mtrx_filled = [roi_idx_mtrx_filled; repelem(list_cc(cc),numel(r_f))' r_f c_f];
            end
        else, cc=cc-1;
        end
    end
end
end