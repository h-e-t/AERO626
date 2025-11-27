function fname = datedFilename(prefix, ext, folder, fmt)
    arguments
        prefix 
        ext 
        folder 
        fmt 
    end
    % datedFilename  Generate a filename based on current date/time.
    % 
    % fname = datedFilename()                      -> 'file_YYYYMMDD_HHMMSS.txt' in pwd
    % fname = datedFilename(prefix, ext)           -> custom prefix and extension
    % fname = datedFilename(prefix, ext, folder)   -> place file in folder
    % fname = datedFilename(prefix, ext, folder, fmt) -> custom datetime format
    %
    % Defaults:
    %   prefix = 'file'
    %   ext    = '.txt'
    %   folder = pwd
    %   fmt    = 'yyyyMMdd_HHmmss'   (examples: 20251126_153045)
    
    if nargin < 1 || isempty(prefix), prefix = 'file'; end
    if nargin < 2 || isempty(ext), ext = '.txt'; end
    if nargin < 3 || isempty(folder), folder = pwd; end
    if nargin < 4 || isempty(fmt), fmt = 'yyyyMMdd_HHmmss'; end
    
    % Ensure extension starts with a dot
    if ext(1) ~= '.'
        ext = '.' + ext
    end

    
    % Build base name
    ts = char(datetime('now','Format',fmt));   % timestamp string
    base = prefix + "_" + ts + ext;
    fullpath = folder + base; 
    
    % If file exists, append an index: _1, _2, ...
    if exist(fullpath, 'file')
        i = 1;
        while true
            candidate = fullfile(folder, sprintf('%s_%d%s', base, i, ext));
            if ~exist(candidate, 'file')
                fullpath = candidate;
                break
            end
            i = i + 1;
        end
    end
    
    fname = fullpath;
end