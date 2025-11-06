%[text] # 1 ICESat-2 ATL03 Data Access
%[text] ## 1.1 Get coastal polygons
% define parameters
reso = 'c';
bbox = [-180 180; -60 66];
bufwidth = {3 .25}; % {outer inner} degrees
tol = 0.1; % degrees
min_area = 100; % km2

% specify query points manually to remove land polygons
% latq = [];
% lonq = [];
latq = [65.95 65.95 65.95 65.98];
lonq = [-174 -162 44 72];

% clip into 3 parts
% clip_lat = [];
clip_lat = [-27 27];
save_flag = false;

shapeoutC = make_coast_polygon(reso, bbox, bufwidth, tol, 'min_area', min_area, ... %[output:group:309f7b4b] %[output:592f385c]
                                                          'latq', latq, ... %[output:592f385c]
                                                          'lonq', lonq, ... %[output:592f385c]
                                                          'clip_lat', clip_lat, ... %[output:592f385c]
                                                          'save_flag', save_flag); %[output:group:309f7b4b] %[output:592f385c]
%%
%[text] ## 1.2 Download ICESat-2 ATL03 data
% using python script `download_icesat2.py`
%%
%[text] ## 1.3 Analyze ICESat-2 ATL03 data
main_analyze()
%%
%[text] ## 1.4 Visualization
datestr = '201810';
version = '006';
save_flag = true;

plot_is2_swh(datestr, version, save_flag);

%[appendix]{"version":"1.0"}
%---
%[metadata:styles]
%   data: {"code":{"fontSize":"12"},"heading1":{"fontSize":"18"},"heading2":{"fontSize":"14"},"heading3":{"fontSize":"12"},"normal":{"fontFamily":"Arial","fontSize":"12"},"title":{"fontSize":"24"}}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:592f385c]
%   data: {"dataType":"text","outputData":{"text":"Simplify coastal polygon (Npts=3956) with tol=0.10 deg, err=1.1%.\n","truncated":false}}
%---
