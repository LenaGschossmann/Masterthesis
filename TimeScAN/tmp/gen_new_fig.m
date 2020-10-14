function gen_new_fig()

% Declare globally shared variables
global POSITIONFIG figINFO FIGCOUNTER

FIGCOUNTER = FIGCOUNTER+1;
figure('Position', POSITIONFIG)
set(gcf, 'UserData', FIGCOUNTER);
set(gcf, 'Name', sprintf('Figure #%i', FIGCOUNTER));

figINFO(FIGCOUNTER).IDs = FIGCOUNTER;
figINFO(FIGCOUNTER).name = sprintf('Figure #%i', FIGCOUNTER);
figINFO(FIGCOUNTER).saved = false;
update_fig();

end