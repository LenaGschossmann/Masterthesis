function gen_new_fig()
% Generates overview of raw, averaged, and dFoF image;
% the actual plotting of these subplots happens by calling update_fig()

% Declare globally shared variables
global POSITIONFIG figINFO FIGCNTER

FIGCNTER = FIGCNTER+1;
figure('Position', POSITIONFIG)
set(gcf, 'UserData', FIGCNTER);
set(gcf, 'Name', sprintf('Figure #%i', FIGCNTER));

% Save figure-specific parameters
figINFO(FIGCNTER).IDs = FIGCNTER;
figINFO(FIGCNTER).name = sprintf('Figure #%i', FIGCNTER);
figINFO(FIGCNTER).saved = false;
update_fig();

end