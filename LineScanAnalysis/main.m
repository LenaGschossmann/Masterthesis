%% ** Script for line scan analysis **

% Declare globally shared variables
global SCRSZ WHFIG WHTRCFIG POSITIONFIG POSITIONTRCFIG POSITIONROISELECT FONTSIZE SCMAPDDITEMS...
    SCMAP WINSZDDITEMS WINSZ PLOTRANGE CLIMRAW CLIMUI TRCYLABEL TRCYLABELDFOF TRCXLABEL TRCPLOTCOL1...
    TRCPLOTCOL2 TRCALPHA THRESHLW THRESHALPHA EVPLOTCOL PLOTPEAKS RISETIME THRESHOLD PEAKTHRESHOLD SMTHWIN...
    figINFO roiINFO traceINFO ROICNTID FIGCOUNTER CURRFILE...
    SAVEPARENT IMHW FTIME COMPOSITE2D SAVEPATH FULLFILENAMES NUMFILES FNAME IMMETA

% Initiate variables
FIGCOUNTER = 0;
ROICNTID = 1;
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
TRCALPHA = 0.25;
THRESHLW = 1;
THRESHALPHA = 0.7;
EVPLOTCOL = [0 0.1882 0.1059];
PLOTPEAKS = true; % otherwise plot THRESHOLD crossings
RISETIME = inputdlg('Upper limit sensor rise time [s]:','Sensor info', [1 60], {'0.2'}); RISETIME = str2double(RISETIME{1});
THRESHOLD = 2;
PEAKTHRESHOLD = THRESHOLD+0.5;
SMTHWIN = 0.1; % smoothing for peak detection [s]
roiINFO = struct(); roiINFO(1).name = []; roiINFO(1).mask = []; roiINFO(1).position = []; roiINFO(1).ID = []; roiINFO(1).selected = []; roiINFO(1).saved = []; roiINFO(1).mode = []; roiINFO(1).PLOTRANGE = [];
figINFO = struct('IDs', [], 'name', [], 'PLOTRANGE',[], 'cSCMAP',[], 'csclimits', [], 'avwinsize',[], 'saved',[]);
traceINFO = struct('figID',[], 'fig_params',[],'roiID',[], 'binned_roi_av',[],'dFoF_roi_av',[], 'timestamp',[], 'save',[], 'currmode',[], 'showtot', []);
SAVEPARENT = [];

close all;
ini_ctrl_box();
