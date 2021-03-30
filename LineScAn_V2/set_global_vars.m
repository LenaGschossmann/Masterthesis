function set_global_vars()

% Declare globally shared variables
global SCRSZ WHFIG WHTRCFIG POSITIONFIG POSITIONTRCFIG POSITIONROISELECT FONTSIZE SCMAPDDITEMS...
    SCMAP WINSZDDITEMS WINSZ PLOTRANGE CLIMRAW CLIMUI TRCYLABEL TRCYLABELDFOF TRCXLABEL TRCPLOTCOL1...
    TRCPLOTCOL2 TRCALPHA THRESHLW THRESHALPHA EVPLOTCOL PLOTPEAKS...
    figINFO roiINFO traceINFO FIGCNTER ROICNTER CURRFILE...
    SAVEPARENT dFWIN FNAME PARAMS SMTHWIN EVDATA

FIGCNTER = 0;
ROICNTER = 0;
CURRFILE = 1;
SCRSZ = get(0, 'Screensize'); SCRSZ = SCRSZ(3:4);
WHFIG = [800 600];
WHTRCFIG = [1000 600];
POSITIONFIG = [100 SCRSZ(2)-WHFIG(2)-100 WHFIG];
POSITIONTRCFIG = [100 SCRSZ(2)-WHFIG(2)-100 WHTRCFIG];
POSITIONROISELECT = [0.05 0.1 0.9 0.85];
FONTSIZE = 8;
SCMAPDDITEMS = {'gray', 'hot', 'summer', 'winter', 'parula'};
SCMAP = SCMAPDDITEMS{1};
WINSZDDITEMS = {'4', '8', '16', '32'};
WINSZ = str2double(WINSZDDITEMS{1}); % Line average
PLOTRANGE = [1 1000];
CLIMRAW = [0 255];
CLIMUI = [-3 3];
TRCYLABEL = 'Av. Intensity [a.u.]';
TRCYLABELDFOF = 'DeltaF/F';
TRCXLABEL = 'time [s]';
TRCPLOTCOL1 = [0.6510 0.0667 0.2784];
TRCPLOTCOL2 = [0.0471 0.3137 0.4588];
TRCALPHA = 0.5;
THRESHLW = 1;
THRESHALPHA = 0.7;
dFWIN = 0.5; % [s]
EVPLOTCOL = [0 0.1882 0.1059];
FNAME = '';
PLOTPEAKS = true; % otherwise plot THRESHOLD crossings
SMTHWIN = 0.1; % smoothing for peak detection [s]
% roiINFO = struct(); roiINFO(1).name = []; roiINFO(1).mask = []; roiINFO(1).position = []; roiINFO(1).ID = []; roiINFO(1).selected = []; roiINFO(1).saved = []; roiINFO(1).mode = []; roiINFO(1).PLOTRANGE = [];
roiINFO = struct('name', [], 'mask', [], 'position', [], 'ID', [], 'selected', [], 'saved', [], 'mode', [], 'plotrange', []);
figINFO = struct('IDs', [], 'name', [], 'PLOTRANGE',[], 'cSCMAP',[], 'csclimits', [], 'avwinsize',[], 'saved',[]); % Saves parameters associated with a certain overview figure
traceINFO = struct('figID',[], 'fig_params',[],'roiID',[], 'binned_roi_av',[],'dFoF_roi_av',[], 'timestamp',[], 'save',[], 'currmode',[], 'showtot', []);
SAVEPARENT = [];
PARAMS =  get_predefined_params('ev_perc_thresh',97, 'critical_fac',1.7, 'safe_fac', 1.3);
EVDATA = cell(1,12);
end