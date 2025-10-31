function main_analyze()
% set parameters
yList = 2018;
mList = 10;
oceanSegLen = 10.24;
depth = -10;
distance = [];
save_flag = 1;

for year = yList
    for month = mList
        tic;
        fprintf('+-----------------------------------------------------------------------------+\n');
        fprintf('Start processing ICESat-2 ATL03 granules for %04d-%02d...\n', year, month);
        datestr = sprintf('%04d%02d**', year, month);
        S = analyze_icesat2_ocean(datestr,  'oceanSegLen', oceanSegLen, ...
                                            'depth', depth, ...
                                            'distance', distance, ...
                                            'save_flag', save_flag);
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