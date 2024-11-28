function varargout = errorstats(x,varargin)
% [MSG,MBE,RMS,MAD] = ERRORSTATS(X) - return statistics for X, along with a summary message
% [MSG,nMBE,nRMS,nMAD] = ERRORSTATS(X,Y,NORM) - return statistics for (X-Y)/NORM(X)
% ERRORSTATS(X,..) - print summary message only
% 
%  ..,'delim',',' - set message delimiter
%  ..,'units','' - set message units. Using 'units','%' will premultiply stats by 100
%  ..,'format','%0.2f' 
%
% EXAMPLES:
%    x = rand(100,1); y = rand(100,1);
%    errorstats(x)
%    errorstats(x,y,'units',' m/s')
%    errorstats(x,y,@mean,'units','%')
%    errorstats(x,[],'units','%','delim',newline())

    opt.y = 0;
    opt.norm = 1;
    opt.delim = ', ';
    opt.tags = {'RMSE','MBE','MAD'};
    opt.units = '';
    opt.format = '%0.2f';
    
    [opt,~,isdef] = parseoptions(varargin,{'omitnan'},opt,'dealrest',2);
    validateattributes(x,'numeric',{'vector','real'});
    if isempty(opt.y), y = 0;
    else
        y = opt.y;
        if isscalar(opt.y)
            validateattributes(y,'numeric',{'scalar','real'});
        else
            validateattributes(y,'numeric',{'vector','real','size',size(x)});
        end
    end
    e = x-y;
    if opt.omitnan
        f = ~isnan(e);
        if ~all(f)
            x = x(f);
            e = e(f);
        end
    end
    
    if isa(opt.norm,'function_handle')
        n = opt.norm(x);
        if isdef.tags, opt.tags = cellfun(@(x) ['n' x],opt.tags,'unif',0); end
    else, n = opt.norm;
    end
    if isscalar(n)
        validateattributes(n,'numeric',{'scalar','real'});
    else
        validateattributes(n,'numeric',{'vector','real','size',size(x)});
    end

    if nargin < 3, n = 1; end
    
    MBE = mean(e)./n;
    RMS = hypot(MBE,std(e)./n);
    MAD = mad(e)./n;
    
    x = [RMS,MBE,MAD];
    if isequal(opt.units,'%'), x = x*100; end
    msg = arrayfun(@(j) sprintf(['%s = ' opt.format '%s'],opt.tags{j},x(j),opt.units),1:3,'unif',0);
    msg = strjoin(msg,opt.delim);
    if nargout == 0
        fprintf('%s\n',msg);
    else
        varargout = {msg,MBE,RMS,MAD};
    end
end