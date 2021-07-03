function yval = lut_create(xval, lut_params)

mv = lut_params.lut_log_maxval-lut_params.lut_log_minval;
k = lut_params.lut_log_growthrate;
mp = lut_params.lut_log_midpoint;

yval = lut_params.lut_log_minval + (mv ./ (1 + exp(-k .*(xval - mp))));

end


% mp=25;
% xval = 1:50;
% mv = 1.5; gr= 0.2; sty = 5; plx = 1 + xval./100;
% % mv = 2; gr = 0.3; sty = 6; plx = 1.25 + xval./100;
% % function yval = lut_create(xval, mv, gr, mp)
% yval = sty + (mv./(1+exp(-gr.*(xval-mp))));
% % end
% figure, plot(plx, yval), ylim([1 10]);
% hold on, plot(xval, yval)