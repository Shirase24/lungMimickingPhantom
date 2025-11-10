%% =======================================================================
%  generateLungPhantom_pro.m
%  Developer-grade 2D Lung Phantom Generator (toolbox-free, error-handled)
%  Author: OpenAI ChatGPT (auto-generated)
%  Target: Simulation pipelines (openEMS/CST/MEEP/gprMax) + teaching demos
% ========================================================================

function generateLungPhantom_pro
    %GENERATELUNGPHANTOM_PRO Build and export a heterogeneous 2D lung phantom.
    %   This script creates a label map with simple analytic shapes that mimic
    %   thoracic anatomy, assigns nominal electrical properties, produces
    %   visualizations, exports data products, and can synthesize a breathing
    %   animation. All dependencies are available in base MATLAB.

    clc; close all;
    t0 = tic; logmsg('--- Lung Phantom: start ---');

    %% ---------------- CONFIG (edit here) ----------------
    cfg = struct();
    cfg.H = 256;                  % image height  (pixels)
    cfg.W = 256;                  % image width   (pixels)
    cfg.freqGHz = 1.0;            % (metadata only)
    cfg.seed = 123;               % RNG seed for reproducibility

    % Electrical props (~1 GHz placeholders; replace with literature when needed)
    cfg.eps = struct('AIR',1.0006,'FAT',5.0,'MUSCLE',52.0,'LUNG',20.0,'HEART',60.0,'BONE',12.0,'PARENCHYMA',25.0);
    cfg.sig = struct('AIR',0.00,  'FAT',0.05,'MUSCLE',0.80,'LUNG',0.40,'HEART',1.00,'BONE',0.20,'PARENCHYMA',0.60);

    % Heterogeneity (lung texture) parameters
    cfg.alpha = struct('min',0.35,'max',0.70,'sigma_px',6,'boost',1.00);

    % Visualization ranges (keep consistent across runs)
    cfg.range = struct('epsr',[1 65],'sigma',[0 1.2]);

    % Lesion example (set radius<=0 to disable)
    cfg.lesion = struct('enable',true,'relRadius',0.08,'centerOffset',[0.18 -0.05], ...
                        'epsrFactor',1.35,'sigmaFactor',1.45);

    % Output
    cfg.outDir  = "output_phantom";
    cfg.pngName = "phantom_epsr_sigma.png";
    cfg.matName = "phantom_data.mat";
    cfg.csvPrefix = "phantom_";
    cfg.makeMP4 = true;           % breathing animation
    cfg.mp4Name = "phantom_breathing.mp4";
    cfg.mp4_fps = 20; cfg.mp4_seconds = 6; cfg.mp4_bpm = 15; cfg.mp4_amp = 0.10;

    rng(cfg.seed);

    % Labels
    LABEL = struct('AIR',0,'FAT',1,'MUSCLE',2,'LUNG',3,'HEART',4,'BONE',5);

    % Validate config
    mustBeScalarPosInt(cfg.H,'cfg.H');
    mustBeScalarPosInt(cfg.W,'cfg.W');
    mustBeRange(cfg.range.epsr,'cfg.range.epsr');
    mustBeRange(cfg.range.sigma,'cfg.range.sigma');

    % Ensure output dir
    ensureDir(cfg.outDir);

    %% ---------------- BUILD PHANTOM ----------------
    [labels, lungsMask] = buildLabelMap(cfg, LABEL);
    [epsr, sigma, alpha] = assignProps(cfg, labels, lungsMask, LABEL);
    [epsr, sigma] = addDefaultLesion(cfg, epsr, sigma, lungsMask);

    %% ---------------- VISUALIZE + SAVE ----------------
    fig = plotClean(cfg, labels, epsr, sigma, LABEL);
    safeExportPNG(fig, fullfile(cfg.outDir,cfg.pngName), 220);

    safeSave(fullfile(cfg.outDir,cfg.matName), labels, epsr, sigma, alpha, cfg);
    safeWriteCSV(fullfile(cfg.outDir, cfg.csvPrefix+"labels.csv"), labels);
    safeWriteCSV(fullfile(cfg.outDir, cfg.csvPrefix+"epsr.csv"),   epsr);
    safeWriteCSV(fullfile(cfg.outDir, cfg.csvPrefix+"sigma.csv"),  sigma);

    %% ---------------- ANIMATION (optional) ----------------
    if cfg.makeMP4
        makeBreathingMP4(cfg, epsr, alpha, lungsMask);
    end

    logmsg('--- Lung Phantom: done (%.3f s) ---', toc(t0));
end

%% ======================= CORE FUNCTIONS =======================

function [labels, lungs] = buildLabelMap(cfg, LABEL)
    H = cfg.H; W = cfg.W;

    % Normalized coordinate system (-1 .. 1) centered in the thorax
    [xIdx, yIdx] = meshgrid(1:W, 1:H);
    x = (xIdx - (W+1)/2) / (W/2);
    y = (yIdx - (H+1)/2) / (H/2);

    labels = uint8(zeros(H, W) + LABEL.AIR);

    % Body contour and subcutaneous fat layer (elliptical torso)
    bodyOuter = ((x/0.95).^2 + ((y+0.05)/1.05).^2) <= 1;
    labels(bodyOuter) = LABEL.FAT;

    % Muscle band between outer torso and inner thoracic cavity
    thoraxOuter = ((x/0.82).^2 + ((y+0.02)/0.98).^2) <= 1;
    thoraxInner = ((x/0.64).^2 + ((y+0.02)/0.86).^2) <= 1;
    muscleShell = thoraxOuter & ~thoraxInner;
    labels(muscleShell) = LABEL.MUSCLE;

    % Mediastinal soft tissue in the centre between the lungs
    mediastinum = ((x/0.30).^2 + ((y+0.02)/0.95).^2) <= 1 & (x > -0.08);
    labels(mediastinum) = LABEL.MUSCLE;

    % Spine (posterior column)
    spine = (abs(x) <= 0.08) & (((y+0.10)/0.88).^2 <= 1) & (y > -0.75);
    labels(spine) = LABEL.BONE;

    % Clavicles (superior bone bridges)
    clavLeft  = ((x+0.45)/0.18).^2 + ((y+0.72)/0.06).^2 <= 1 & y < -0.58;
    clavRight = ((x-0.45)/0.18).^2 + ((y+0.72)/0.06).^2 <= 1 & y < -0.58;
    labels(clavLeft | clavRight) = LABEL.BONE;

    % Ribs approximated as curved bone arcs wrapping the thorax
    ribCenters = linspace(-0.20, 0.60, 5);
    for k = 1:numel(ribCenters)
        y0 = ribCenters(k);
        outer = ((x/0.92).^2 + ((y - y0)/0.09).^2) <= 1;
        inner = ((x/0.70).^2 + ((y - y0)/0.07).^2) <= 1;
        rib = outer & ~inner & (y > y0-0.10) & (y < y0+0.18) & bodyOuter;
        rib = rib & ~spine; % keep mediastinum clear
        labels(rib) = LABEL.BONE;
    end

    % Superior airway (trachea) carved out of mediastinum
    trachea = (abs(x+0.02) <= 0.05) & (y < -0.55) & (y > -0.95);
    labels(trachea) = LABEL.AIR;

    % Diaphragm curvature limiting inferior extent of the lungs
    diaphragm = 0.60 + 0.06*cos(pi*x/1.7);

    % Left lung (patient left / image right) with cardiac notch
    lungL = (((x+0.30)/0.36).^2 + ((y+0.05)/0.78).^2) <= 1;
    lungL = lungL & (x < 0.10) & (y <= diaphragm);
    cardiacNotch = (((x+0.02)/0.20).^2 + ((y-0.02)/0.28).^2) <= 1 & (x > -0.15);
    lungL = lungL & ~cardiacNotch;

    % Right lung larger, spanning more inferiorly
    lungR = (((x-0.22)/0.40).^2 + ((y+0.08)/0.82).^2) <= 1;
    lungR = lungR & (x > -0.32) & (y <= diaphragm + 0.02);

    % Hilum clearance near mediastinum
    hilum = (abs(x) < 0.06) & (y > -0.25) & (y < 0.25);
    lungL = lungL & ~hilum;
    lungR = lungR & ~hilum;

    lungs = (lungL | lungR) & bodyOuter;
    labels(lungs) = LABEL.LUNG;

    % Heart occupying mediastinum and intruding into left lung
    heartUpper = (((x-0.02)/0.26).^2 + ((y+0.10)/0.28).^2) <= 1 & (y < 0.20);
    heartLower = (abs(x-0.05) + 0.55*(y-0.05)) <= 0.38 & (y >= -0.10) & (y <= 0.45);
    heart = (heartUpper | heartLower) & bodyOuter;
    labels(heart) = LABEL.HEART;

    % Ensure lungs do not include the heart or spine regions
    lungs = (labels == LABEL.LUNG);
end

function [epsr, sigma, alpha] = assignProps(cfg, labels, lungs, LABEL)
    % Base property assignment
    H=cfg.H; W=cfg.W; epsr=zeros(H,W); sigma=zeros(H,W);
    tissueOrder = {'AIR','FAT','MUSCLE','LUNG','HEART','BONE'};
    for k = 1:numel(tissueOrder)
        key = tissueOrder{k};
        valEps = cfg.eps.(key);
        valSig = cfg.sig.(key);
        mask = labels == LABEL.(key);
        epsr(mask) = valEps;
        sigma(mask) = valSig;
    end

    % Heterogeneity field alpha (0..1) only in lungs
    n = smooth2(randn(H,W), cfg.alpha.sigma_px);
    n = (n - min(n(:)));
    if max(n(:)) > 0
        n = n / max(n(:));
    end
    alpha = (cfg.alpha.min + (cfg.alpha.max-cfg.alpha.min)*n) * cfg.alpha.boost;
    alpha = min(max(alpha,0),1);
    alpha(~lungs) = 0;

    % Blend lungs between air and parenchyma
    eps_air=cfg.eps.AIR; eps_par=cfg.eps.PARENCHYMA;
    sig_air=cfg.sig.AIR; sig_par=cfg.sig.PARENCHYMA;
    epsr(lungs)  = alpha(lungs).*eps_par + (1-alpha(lungs)).*eps_air;
    sigma(lungs) = alpha(lungs).*sig_par + (1-alpha(lungs)).*sig_air;
end

function [epsr,sigma] = addDefaultLesion(cfg, epsr, sigma, lungs)
    if ~isfield(cfg,'lesion') || isempty(cfg.lesion) || ~cfg.lesion.enable
        return;
    end

    radius = max(1, round(cfg.lesion.relRadius * min(cfg.H,cfg.W)));
    if radius <= 0
        return;
    end

    [X,Y] = meshgrid(1:cfg.W,1:cfg.H);
    cx = cfg.W/2 + cfg.lesion.centerOffset(1)*cfg.W;
    cy = cfg.H/2 + cfg.lesion.centerOffset(2)*cfg.H;
    lesionMask = ((X-cx).^2 + (Y-cy).^2) <= radius^2;
    lesionMask = lesionMask & lungs;
    if ~any(lesionMask(:))
        logmsg('Lesion mask empty; skipping default lesion.');
        return;
    end

    lungEps = cfg.eps.PARENCHYMA;
    lungSig = cfg.sig.PARENCHYMA;
    epsr(lesionMask) = lungEps * cfg.lesion.epsrFactor;
    sigma(lesionMask) = lungSig * cfg.lesion.sigmaFactor;
    logmsg('Lesion injected: radius=%d px, epsr x%.2f, sigma x%.2f', radius, cfg.lesion.epsrFactor, cfg.lesion.sigmaFactor);
end

function fig = plotClean(cfg, labels, epsr, sigma, LABEL)
    if exist('turbo','file') == 2
        turboMap = turbo(256);
    else
        turboMap = getTurboFallback(256);
    end

    labNames = {'Air','Fat','Muscle','Lung','Heart','Bone/Ribs'};
    labColors = [0.90 0.90 0.90;  % air
                 1.00 0.89 0.65;  % fat
                 0.80 0.10 0.10;  % muscle
                 0.50 0.80 0.85;  % lung
                 0.80 0.20 0.50;  % heart
                 0.65 0.65 0.65]; % bone/spine
    labCmap = interp1(1:size(labColors,1), labColors, linspace(1,size(labColors,1),256));

    fig = figure('Name','Lung Phantom','Position',[120 120 1220 420]);
    tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

    % 1) Labels
    nexttile; imagesc(labels); axis image off;
    colormap(gca, labCmap);
    c=colorbar; c.Ticks=0:5; c.TickLabels=labNames; c.Box='off';
    title('Tissue Map','FontWeight','bold');

    % 2) Permittivity
    nexttile; imagesc(epsr, cfg.range.epsr); axis image off;
    colormap(gca, turboMap); colorbar; title('\epsilon_r @ 1 GHz','FontWeight','bold');

    % 3) Conductivity
    nexttile; imagesc(sigma, cfg.range.sigma); axis image off;
    colormap(gca, hot(256)); colorbar; title('\sigma (S/m)','FontWeight','bold');

    set(gcf,'Color','w');
end

function makeBreathingMP4(cfg, epsr, alpha, lungs)
    try
        v = VideoWriter(fullfile(cfg.outDir,cfg.mp4Name),'MPEG-4');
    catch
        logmsg('(!) VideoWriter MPEG-4 not available on this MATLAB/OS. Skipping MP4.');
        return;
    end
    v.FrameRate = cfg.mp4_fps; open(v);
    t = linspace(0, cfg.mp4_seconds, cfg.mp4_seconds*cfg.mp4_fps);
    bpm = cfg.mp4_bpm; amp = cfg.mp4_amp;

    eps_air=cfg.eps.AIR; eps_par=cfg.eps.PARENCHYMA;

    logmsg('MP4: rendering breathing animation ...');
    for k=1:numel(t)
        a = alpha .* (1 + amp*sin(2*pi*(bpm/60)*t(k)));
        a = min(max(a,0),1);
        e = epsr;
        e(lungs) = a(lungs).*eps_par + (1-a(lungs)).*eps_air;

        im = mat2gray(e, cfg.range.epsr);
        fr = ind2rgb(gray2ind(im,256), turboOrFallback(256));
        writeVideo(v, im2frame(im2uint8(fr)));
    end
    close(v);
    logmsg('MP4 saved -> %s', fullfile(cfg.outDir,cfg.mp4Name));
end

%% ======================= UTILITIES =======================

function ensureDir(p)
    if ~exist(p,'dir')
        ok = mkdir(p);
        if ~ok, error('Failed to create output dir: %s', p); end
    end
end

function mustBeScalarPosInt(x, name)
    if ~(isscalar(x) && isnumeric(x) && x>0 && mod(x,1)==0)
        error('%s must be a positive integer scalar.', name);
    end
end

function mustBeRange(r, name)
    if ~(isnumeric(r) && numel(r)==2 && r(1)<r(2))
        error('%s must be a 1x2 numeric [min max] with min<max.', name);
    end
end

function safeExportPNG(fig, filename, dpi)
    try
        exportgraphics(fig, filename, 'Resolution', dpi);
        logmsg('PNG saved -> %s', filename);
    catch ME
        warning('exportgraphics failed (%s). Using print() fallback.', ME.message);
        try
            set(fig,'PaperPositionMode','auto'); print(fig, filename, '-dpng', sprintf('-r%d',dpi));
            logmsg('PNG (print fallback) -> %s', filename);
        catch ME2
            warning('print fallback failed: %s', ME2.message);
        end
    end
end

function safeSave(matFile, labels, epsr, sigma, alpha, cfg)
    try
        save(matFile,'labels','epsr','sigma','alpha','cfg','-v7.3');
        logmsg('MAT saved -> %s', matFile);
    catch ME
        warning('save -v7.3 failed (%s). Retrying without -v7.3.', ME.message);
        try
            save(matFile,'labels','epsr','sigma','alpha','cfg');
            logmsg('MAT saved (fallback) -> %s', matFile);
        catch ME2
            warning('MAT save failed: %s', ME2.message);
        end
    end
end

function safeWriteCSV(path, A)
    try
        if exist('writematrix','file') == 2
            writematrix(A, path);
        else
            csvwrite(path, A); %#ok<CSVWR>
        end
        logmsg('CSV saved -> %s', path);
    catch ME
        warning('CSV write failed for %s: %s', path, ME.message);
    end
end

function S = smooth2(A, sigma_px)
    % Gaussian-like separable blur (toolbox-free)
    if sigma_px<=0, S=A; return; end
    r = max(1, ceil(3*sigma_px));
    x = -r:r; g = exp(-(x.^2)/(2*sigma_px^2)); g = g/sum(g);
    S = conv2(conv2(A, g, 'same'), g', 'same');
end

function logmsg(fmt, varargin)
    fprintf('%s  %s\n', datestr(now,'HH:MM:SS'), sprintf(fmt, varargin{:}));
end

function cmap = turboOrFallback(n)
    if exist('turbo','file') == 2
        cmap = turbo(n);
    else
        cmap = getTurboFallback(n);
    end
end

function cmap = getTurboFallback(n)
    % Approximate turbo colormap coefficients (MATLAB File Exchange variant)
    % Reference: Anton Semechko, 2020 (public domain approximation)
    t = linspace(0,1,n)';
    r = 0.135 - 0.13*t + 0.645*t.^2 + 0.365*t.^3 - 1.17*t.^4 + 0.725*t.^5;
    g = 0.12  + 0.885*t - 0.305*t.^2 - 1.8*t.^3 + 2.77*t.^4 - 1.3*t.^5;
    b = 0.45  + 0.17*t - 1.6*t.^2 + 2.31*t.^3 - 1.5*t.^4 + 0.34*t.^5;
    cmap = max(min([r g b],1),0);
end

