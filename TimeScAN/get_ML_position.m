function [roidistances] = get_ML_position(traceindices)

global SOMAINFO LOCAREAINFO


somares = SOMAINFO(end,:); somalocs = SOMAINFO(1:end-1,:);
roilocs = LOCAREAINFO(traceindices,:);

% Calculate location of ROIs in z-stack projection
roilocs(:,1) = 0.5* somares(1) + roilocs(:,1);
roilocs(:,2) = 0.5* somares(2) + roilocs(:,2);

% Fit line through somata
tbl.x = somalocs(:,1);
tbl.y =  somalocs(:,2);
lnfit = fitlm(tbl.x, tbl.y, 'linear');
coeffs = lnfit.Coefficients.Estimate;
numcoeffs = numel(coeffs);
c1 = coeffs(2); % X Coeff.
c2 = 1; % Y Coeff.
c3 = coeffs(1); % Intercept

% Calculate distance of ROI centroids in ML of Cb from somata
numroictr = size(roilocs,1);
roidistances = zeros(numroictr,1);
for iR = 1:numroictr
    tmpctr = roilocs(iR,1:2);
    roidistances(iR) = abs(c1*tmpctr(1) - c2*tmpctr(2) + c3)/ sqrt(c1^2 + c2^2);
end

end
