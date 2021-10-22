function [data,delimiter,formatstring] = readtxtfile(filename,varargin)
% DATA = READTXTFILE(FILENAME,PARAMNAMES)
% DATA = READTXTFILE(..,NAME,VALUE..)
% Attempt to read a text file of the form:
%
%     # PARAMETER: VALUE
%     # ..
%     # COMMENTS
%     [#] HEADER
%     DATA
%
% If cell-string PARAMNAMES is provided (default is {}), all COMMENT/HEADER lines will be scanned 
% for propperty-value pairs of the form X: VALUE, where X matches a given PARAMNAMES{j,1}. If
% PARAMNAMES has a second column, matches will be saved with field-names PARAMNAMES{j,2}.
%
% FILENAME (path & name of file to import) can be omitted for UI import via PICKFILE()
%
% Output:
% DATA.PARAMS: structure with fields in PARAMNAMES for which matches were found in commented lines.
%   All values will be RETURNED AS STRINGS, no parsing effort whatsoever.
% DATA.COMMENTS: cell-array of strings with any lines (at the beginning of the file) that start 
%   with '#'. Spaces at the beginning, and spaces/delimiters at the end will be removed.
%   Any lines matching a propperty-value pair (see above) will also be removed.
% DATA.HEADER: 1·N cell-array of column headers obtained from the last row of COMMENTS or the 
%   first row of DATA. The header line must use the same delimiters and have the same number of 
%   columns as DATA, if the row is not commented, it also must contain less numeric fields 
%   than DATA.
% DATA.DATA: 1·N cell-array of data-columns. Columns with numeric data (as tested by STR2DOUBLE)
%   will be numeric vectors, all others will be a cell-vector of strings.
%
% If the file format doesn't match that specified above, the result will be an empty structure.
%
% [DATA,DELIM,FORMATSTR] = READTXTFILE(...) returns the used delimiter and format-string.
%
% DATA = READTXTFILE(..,'comment',c) use custom comment-character, instead of #
% DATA = READTXTFILE(..,'delim',c) use custom delimiter-character. By default, the first appearance
%        of {tab, semicolon, comma, or space} after commented lines will be used as delimiter.
%        NOTE: both HEADER and DATA should be separater by the same delimiter.
% DATA = READTXTFILE(..,'formatstring',s), use custom format-string s in TEXTSCAN, instead of one
%        automatically detected from first non-header line. Only used in Data, not for Header.
% DATA = READTXTFILE(..,'ignorecase',true) - case-insensitive match for PARAMNAMES, resulting
%        structure names will be lowercased.
% DATA = READTXTFILE(..,'regexpmatch',true) - use regexp(PAR,STR) to match PARAMNAMES.
% DATA = READTXTFILE(..,'colon',c) - use character c instead of default ':' to separate PARAMETER-
%   VALUE pairs on the file header.
% DATA = READTXTFILE(..,'data',false) - read header only
%
% DATA = READTXTFILE(..,param,value) should also work with almost(*) all parameter-value pairs in
%        TEXTSCAN (e.g. CollectOutput, TreatAsEmpty, DateLocale, etc.). Exceptions are:
%           Delimiter, HeaderLines: ignored, handled by delim, and (indirectly) comment options.
%           CommentStyle: passed, don't use unless there are comments mixed with data
%           MultipleDelimsAsOne: passed, set to TRUE by default
%
% See also: DLMREAD, IMPORTDATA, TEXTSCAN

    % Delimiter and HeaderLines not included
    TS_PARAMS = {'CollectOutput', 'CommentStyle', 'DateLocale', 'EmptyValue', 'EndOfLine',...
        'ExpChars','ReturnOnError','TreatAsEmpty','Whitespace','TextType','MultipleDelimsAsOne'};
    TS_DEFAULTS(1:numel(TS_PARAMS)) = {[]};
    TS_DEFAULTS{end} = true; % MultipleDelimsAsOne defaults true

    % Custom (READTXTFILE-specific) parameters
    RTF_PARAMS = {'delim','comment','ignorecase','formatstring','regexpmatch','colon','data'};
    RTF_DEFAULTS = {'','#',false,'',false,':',true};

    [opt,remargs] = getpairedoptions(varargin,[TS_PARAMS,RTF_PARAMS],[TS_DEFAULTS,RTF_DEFAULTS]);
    assert(numel(remargs)<=1, shortliststr(remargs(2:end),'Unrecognized argument','colon',':'));
    
    if ~isempty(remargs)
        assert(iscellstr(remargs{1}),'Expecting cell-string of parameter names');
        paramnames = remargs{1};
    else, paramnames = {}; 
    end
    if size(paramnames,2) > 2, paramnames = paramnames'; end
    if ~isempty(paramnames)
        assert(any(size(paramnames,2) == [1,2]),'Expecting vector or 2-column cell-array of parameters');
    end
        
    if nargin < 1, filename = pickfile(); end
    
    % Handle TEXTSCAN parameters separately
    tsargs = rmfield(opt,RTF_PARAMS);
    tsargs = [fieldnames(tsargs),struct2cell(tsargs)]';
    tsargs(:,cellfun(@isempty,tsargs(2,:))) = []; % remove empty...
        
    % Attempt to read a generic file with #comments followed by delimited data
    fileID = fopen(filename);
	if fileID < 0, error('getmeteodata:open','Could not open file'); end
    comments = textscan(fileID,[opt.comment '%[^\n\r]'],tsargs{:}); comments = comments{1};
        
    eoc = ftell(fileID); % remember position at end of comments
    firstdatalines = textscan(fileID,'%[^\n\r]',2,tsargs{:});
    
    % Try to guess delimiter, if not set explicitly
    if isempty(opt.delim) && ~isempty(firstdatalines{1})
        opt.delim = [char(9),'; ,']; % opt.delimiter priority: (tab) ; (space) ,
        for j = 1:numel(opt.delim)
            % Guess delimiters from second data line (if possible), in case first is header
            if any(firstdatalines{1}{end}==opt.delim(j)), opt.delim = opt.delim(j); break; end
        end
    end

    % Clean-up comments and try to find parameter-value pairs
    [data,comments] = parsecomments(comments,paramnames,opt);
    delimiter = opt.delim;
    
    % If there's no data that's it (data may be non-empty, if there are comments)
    if isempty(firstdatalines{1}) || isempty(delimiter), fclose(fileID); return; end
    
    % Try to split data lines in tokens
    fields = cellfun(@(s) strsplit(s,delimiter,'CollapseDelimiter',false),firstdatalines{1}(:),'unif',0);
    for j = 1:numel(fields), fields{j} = fields{j}(~cellfun(@isempty,fields{j})); end   % remove empty
    %for j = 1:numel(fields), fields{j}(cellfun(@isempty,fields{j})) = {NaN}; end   % empties -> nans

    % If both lines have different number of tokens, we're lost
    if numel(fields) > 1 && numel(fields{2}) ~= numel(fields{1}), fclose(fileID); return; end

    fields = cat(1,fields{:});
    numfields = ~isnan(str2double(fields)); 
    if ~isempty(opt.TreatAsEmpty)
        % Assume any TreatAsEmpty characters are also numbers
        numfields = numfields | strcmp(fields,opt.TreatAsEmpty);
    end
    fseek(fileID,eoc,'bof'); % move back to end of comments

    if ~isempty(comments) && numel(strsplit(comments{end},delimiter)) == size(fields,2)
    % Header-line is last line of comments
        headers = strsplit(comments{end},delimiter);
        data.comments = data.comments(1:end-1);
    elseif size(numfields,1) > 1 && nnz(numfields(2,:)) > nnz(numfields(1,:)) 
    % Header-line is first line of non-comments
        headers = fields(1,:);
        textscan(fileID,'%[^\n\r]',1,tsargs{:}); % move one row after end of comments
    elseif size(numfields,1) == 1 || nnz(numfields(2,:)) == nnz(numfields(1,:))
    % No header-liine
        headers = {};
    else
        fclose(fileID); return; % No clue, return empty data
    end

    data(1).headers = headers;
    
    % Read delimited data and save everything to structure
    if isempty(opt.formatstring)
        formatstring = repmat({'%f'},1,size(fields,2));
        formatstring(~numfields(end,:)) = {'%s'};
        formatstring = strjoin(formatstring);
    else
        formatstring = opt.formatstring;
    end
    if opt.data
        data(1).data = textscan(fileID,formatstring,'Delimiter',delimiter,tsargs{:});
    end
    fclose(fileID);
end

function [data,comments] = parsecomments(comments,paramnames,opt)
% Remove empty comment lines, remove extra delimiters, and try to match parameter-value pairs.
% Return an empty structure if there are no comments left.

    data = struct('params',{},'comments',{},'headers',{},'data',{});
        
    % Remove trailing delimiters from comments (happens with csv files)
    if isempty(opt.delim), opt.delim = '\s;,'; end % anything
    for j = 1:numel(comments)
       comments{j} = strtrim(strjoin(regexp(comments{j},['([' opt.delim ']*)$'],'split'),''));
    end
    % Clear empty rows
    comments = comments(~cellfun(@isempty,comments));
    
    if isempty(comments), return; end
    
    if isempty(paramnames)
        data(1).comments = comments;
    else
        % Try to match propperty-value pairs, dump the rest to comments.
        [data(1).params,data(1).comments] = findproperties(comments,paramnames,opt);
    end
end

function [S,C] = findproperties(C,paramnames,opt)
% Find property-value pairs in the format 'property_name: value' within the cell-array of strings
% C, and arrange them in a structure with fields S.property_name = value.
% Recognized property names are only PARAMNAMES, any remaining lines will be returned as a
% cell-array of strings in C.
%
% ALL VALUES RETURNED AS STRINGS, no parsing effort whatsoever.

    if isempty(paramnames) || isempty(C), S = struct(); return; end
    
    % Find anything that remotely has the form foo:bar
    S = regexp(C,['^(.+?)' opt.colon '(.+)$'],'tokens');
    validpairs = ~cellfun(@isempty,S);
    
    if ~any(validpairs), S = struct(); return; end
    
    % Arrange into a cell-array {foo,bar;...}
    S = cat(1,S{validpairs}); 
    S = cat(1,S{:});
    S = cellfun(@strtrim,S,'unif',0);
    
    % Filter fields which are in fact recognized PARAMNAMES
    if opt.ignorecase
       paramnames(:,1) = lower(paramnames(:,1));
       S(:,1) = lower(S(:,1));
    end
        
    if opt.regexpmatch
        for j = size(paramnames,1):-1:1
            c(j,:) = ~cellfun(@isempty,regexp(S(:,1),paramnames{j,1})); 
        end
        [p,~] = find(c);                % propperty matched
        c = any(c,1);                   % logical index on S for valid pairs
        validpairs(validpairs) = c; 
    else
        [~,p,c] = intersect(paramnames(:,1),S(:,1),'stable'); % integer index on S for valid pairs
        validpairs(validpairs) = sparse(c,1,true,nnz(validpairs),1);
    end
    
    if size(paramnames,2) == 1
        fieldnames = matlab.lang.makeValidName(S(c,1));
    else
        fieldnames = matlab.lang.makeValidName(paramnames(p,2));
    end
        
    % Build a structure with those fields, dump everything else to S.info
    S = cell2struct(S(c,2),fieldnames);
    C = C(~validpairs);
end
