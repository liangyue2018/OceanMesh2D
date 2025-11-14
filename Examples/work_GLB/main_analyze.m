function main_analyze()
% set parameters
yList = 2018;
mList = 10;
oceanSegLen = 10.24; % [km]
depth_threshold = -10; % [m]
distance_threshold = 3; % [km]
version = '007';
%weight_threshold = 0.2; % for v006
weight_threshold = 0.03; % for v007
bufWin = 20; % [m]
interval = 5; % [m]
method = 'linear';
lackwidth = 20; % [m]
missing_threshold = 0.01;
psd_flag = false;
save_flag = true;

% Parallel processing flag
useParallel = true;

for year = yList
    for month = mList
        tic;
        fprintf('+-----------------------------------------------------------------------------+\n');
        fprintf('Start processing ICESat-2 ATL03 granules for %04d-%02d\n', year, month);
        
        if useParallel
            nDays = eomday(year, month);
            poolobj = gcp('nocreate'); % If no pool, do not create new one.
            if isempty(poolobj)
                c = parcluster('Processes');
                c.NumWorkers = min(nDays, c.NumWorkers);
                parpool(c);
            end
            cleanupObj = onCleanup(@() delete(gcp('nocreate')));

            % use parfor
            parfor day = 1:nDays
                datestr = sprintf('%04d%02d%02d', year, month, day);
                analyze_icesat2_ocean(datestr,  'oceanSegLen', oceanSegLen, ...
                                                'depth_threshold', depth_threshold, ...
                                                'distance_threshold', distance_threshold, ...
                                                'version', version, ...
                                                'weight_threshold', weight_threshold, ...
                                                'bufWin', bufWin, ...
                                                'interval', interval, ...
                                                'method', method, ...
                                                'lackwidth', lackwidth, ...
                                                'missing_threshold', missing_threshold, ...
                                                'psd_flag', psd_flag, ...
                                                'save_flag', save_flag);
            end

            % concatenate results
            concat_results(year, month, version);
        else
            datestr = sprintf('%04d%02d**', year, month);
            analyze_icesat2_ocean(datestr,  'oceanSegLen', oceanSegLen, ...
                                            'depth_threshold', depth_threshold, ...
                                            'distance_threshold', distance_threshold, ...
                                            'version', version, ...
                                            'weight_threshold', weight_threshold, ...
                                            'bufWin', bufWin, ...
                                            'interval', interval, ...
                                            'method', method, ...
                                            'lackwidth', lackwidth, ...
                                            'missing_threshold', missing_threshold, ...
                                            'psd_flag', psd_flag, ...
                                            'save_flag', save_flag);
        end

        elapsedTime = toc;
        fprintf('Finished processing ICESat-2 ATL03 granules for %04d-%02d in %s.\n', year, month, formatElapsedTime(elapsedTime));
        fprintf('+-----------------------------------------------------------------------------+\n');
    end
end
end

function elapsedTimeStr = formatElapsedTime(elapsedTime)
% Format elapsed time (seconds) for display
if elapsedTime < 60
    elapsedTimeStr = sprintf('%.2f seconds', elapsedTime);
elseif elapsedTime < 3600
    elapsedTimeStr = sprintf('%.2f minutes', elapsedTime / 60);
else
    elapsedTimeStr = sprintf('%.2f hours', elapsedTime / 3600);
end
end

function concat_results(year, month, version)
% Concatenate results for a given year and month with specified version
save_dir = fullfile('results', "ATL03_v" + version);
pattern = sprintf('icesat2_atl03_%04d%02d*_v%s.csv', year, month, version);
fstruct = dir(fullfile(save_dir, pattern));
if ~isempty(fstruct)
    names = {fstruct.name}';
    rx = sprintf('^icesat2_atl03_%04d%02d\\d{2}_v%s\\.csv$', year, month, version);
    ok = ~cellfun(@isempty, regexp(names, rx));
    fstruct = fstruct(ok);
end
cfile = fullfile({fstruct.folder}', {fstruct.name}');
TT_month = cell(numel(cfile),1);
for k = 1:numel(cfile)
    opts = detectImportOptions(cfile{k});
    opts = setvartype(opts, {'beam' 'rcs'}, {'string' 'string'});
    opts = setvaropts(opts, 'Time', InputFormat='uuuu-MM-dd HH:mm:ss', TimeZone='UTC');
    TT = readtimetable(cfile{k}, opts);
    TT_month{k} = TT;
    delete(cfile{k});
end
TT_month = vertcat(TT_month{:});

% save
outfnm = fullfile(save_dir, sprintf('icesat2_atl03_%04d%02d_v%s.csv', year, month, version));
writetimetable(TT_month, outfnm, DateLocale='en_US');
end