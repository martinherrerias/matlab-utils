classdef table_ish
% A wrapper for a table/timetable, that returns isfield/fieldnames as a structure would, and most
% importantly, doesn't complain that "Class 'table' is Sealed and may not be used as a superclass"
%
% SUBSREF and SUBSASGN overloads are meant to allow variables inside the table (obj.data) to 
% be accessed as if they were properties, e.g. obj.data.a == obj.a .Notice, however, that the 
% overloaded methods are NOT called inside sub-class definitions (e.g. MeteoData).

properties (Hidden, GetAccess=public, SetAccess=protected)
    data
end
methods
    function F = table_ish(X)
        if nargin == 0, F.data = table.empty; return; end

        if isa(X,'table_ish'), F = X;
        elseif istable(X) || istimetable(X)
            F.data = X;
            
            if isa(F,'MeteoData')
                if istable(X)
                    F.data = table2timetable(F.data,'RowTimes',NaT(size(F.data,1),1));
                end
                F.data.Properties.DimensionNames{1} = 't';
            end
            
            S = F.data.Properties.UserData;
            if ~isempty(S) && isstruct(S)
                
                props = properties(F);
                P = rmfield(S,setdiff(fieldnames(S),props));
                F.data.Properties.UserData = rmfield(S,intersect(fieldnames(S),props));
                
                for j = 1:numel(props)
                   if isfield(P,props{j})
                       F.(props{j}) = P.(props{j});
                   end
                end
            end
        elseif isstruct(X)
            
            fld = fieldnames(X);
            
            props = [{'Properties'};properties(F)];
            P = rmfield(X,setdiff(fld,props));
            X = rmfield(X,intersect(fld,props));
           
            % handle fields with different numbers of rows
            fld = fieldnames(X);
            n = cellfun(@(f) size(X.(f),1),fld);
            if ~all(n == n(1))

                if isa(F,'MeteoData')
                    ia = ismember(fld,MeteoData.varnames('all'));
                    if ~any(ia), ia = n > 1; end
                else
                    ia = n > 1;
                end
                assert(any(ia),'Failed to find common field size');
                N = mode(n(ia));
                R = rmfield(X,fld(n==N));
                X = rmfield(X,fld(n~=N));
            else
                R = struct(); 
            end
            F.data = struct2table(X);
            if ~isempty(fieldnames(R))
               F.data.Properties.UserData = R;
            end
            
            if isa(F,'MeteoData')
                F.data = table2timetable(F.data,'RowTimes',NaT(size(F.data,1),1));
                F.data.Properties.DimensionNames{1} = 't';
            end
        
            for j = 2:numel(props)
               if isfield(P,props{j})
                   F.(props{j}) = P.(props{j});
               end
            end
            
            if isfield(P,'Properties') && ...
               ( isa(P.Properties,'matlab.tabular.TableProperties') || ...
                 isa(P.Properties,'matlab.tabular.TimetableProperties') )

                if isa(P.Properties,'matlab.tabular.TimetableProperties') && istable(F.data)
                    F.data = table2timetable(F.data,'RowTimes',P.Properties.RowTimes);
                end

                [ia,ib] = ismember(fieldnames(F),P.Properties.VariableNames);
                ib = ib(ia);
                if ~any(ia)
                    warning('Unable to use table Properties: non matching names');
                    return;
                end
                
                if ~isempty(P.Properties.VariableUnits)
                    F.data.Properties.VariableUnits(ia) = P.Properties.VariableUnits(ib);
                end
                if ~isempty(P.Properties.VariableDescriptions)
                    F.data.Properties.VariableDescriptions(ia) = P.Properties.VariableDescriptions(ib);
                end

                if ~isempty(P.Properties.CustomProperties)
                    warning_resetter = naptime('MATLAB:structOnObject'); %#ok<NASGU>
                    p = struct(P.Properties.CustomProperties);
                    s = fieldnames(p.perTableProps);
                    for j = 1:numel(s)
                        F.data = addprop(F.data,s{j},'table');
                        F.data.Properties.CustomProperties.(s{j}) = p.perTableProps.(s{j});
                    end
                    s = fieldnames(p.perVarProps);
                    for j = 1:numel(s)
                        F.data = addprop(F.data,s{j},'variable');
                        v = p.perVarProps.(s{j});
                        if isempty(v)
                            F.data.Properties.CustomProperties.(s{j}) = v;
                        elseif all(ia)
                            F.data.Properties.CustomProperties.(s{j}) = v(ib);
                        else
                            F.data.Properties.CustomProperties.(s{j}) = revertfilter(v(ib),ib,2);
                        end
                    end
                end
            end
        end
    end

    function varargout = subsref(obj,s)
        if strcmp(s(1).type,'.') 
            
            if isprop(obj,s(1).subs) % || ismethod(obj,s(1).subs)
            
                if numel(s) > 1 && isprop(obj,s(1).subs) && isa(obj.(s(1).subs),'table_ish')
                    [varargout{1:nargout}] = subsref(obj.(s(1).subs),s(2:end));
                    return;
                end

                % Allow column names as property indices, e.g. x = obj.prop{'key'}
                if isprop(obj,s(1).subs) && ~strcmp(s(1).subs,'data') && numel(s) > 1 &&  ...
                   contains(s(2).type,{'{}','()'}) && hascharkeys(s(2).subs)

                    if isequal(s(2).subs,{':'})
                        x = builtin('subsref',obj,s(1));
                        varargout = x;
                        return;
                    else
                        keys = cellstr(s(2).subs);
                        [ia,ib] = ismember(keys,obj.data.Properties.VariableNames);
                        if ~all(ia)
                            ERR = MException('table_ish:subsasgnkey',...
                                shortliststr(keys(~ia),'Unknown key','quotes',''''));
                            throwAsCaller(ERR);
                        end
                        s(2).subs = {ib};
                    end
                end
            elseif isfield(obj,s(1).subs)
            % Make OBJ.field... calls equivalent to OBJ.data.field...
            
                [varargout{1:nargout}] = subsref(obj.data,s);
                return;
            end
            
            [varargout{1:nargout}] = builtin('subsref',obj,s);
        else
            [varargout{1:nargout}] = builtin('subsref',obj.data,s);
        end
    end
    
    function n = numArgumentsFromSubscript(obj,s,indexingContext)
        if isprop(obj,s(1).subs) && numel(s) > 1 && contains(s(2).type,'{}') ...
             && ~strcmp(s(1).subs,'data') 
            
            if isnumeric(s(end).subs{1}), n = numel(s(2).subs{1});
            elseif isequal(s(end).subs,{':'}), n = numel(builtin('subsref',obj,s(1)));
            else, n = numel(cellstr(s(end).subs));
            end
        elseif any(isfield(obj,s(1).subs))
            n = builtin('numArgumentsFromSubscript',obj.data,s,indexingContext);
        else
            n = builtin('numArgumentsFromSubscript',obj,s,indexingContext);
        end
    end
    
    function obj = subsasgn(obj,s,varargin)

        % Allow subscripted assignment to uninitialized variable
        if isequal(obj,[]), obj = table_ish.empty; end

        if strcmp(s(1).type,'.') && ( isprop(obj,s(1).subs) || ismethod(obj,s(1).subs) )
            
            % Allow column names as property indices, e.g. obj.prop{'key'} = ...
            if numel(s) > 1 && contains(s(2).type,{'{}','()'}) && hascharkeys(s(2).subs)
                
                keys = cellstr(s(2).subs);
                [ia,ib] = ismember(keys,obj.data.Properties.VariableNames);
                if ~all(ia)
                    ERR = MException('table_ish:subsasgnkey',...
                        shortliststr(keys(~ia),'Unknown key','quotes',''''));
                    throwAsCaller(ERR);
                end
                s(2).subs = {ib};
            end
            obj = builtin('subsasgn',obj,s,varargin{:});
        else
            obj.data = builtin('subsasgn',obj.data,s,varargin{:});
        end
    end
    
    function T = table(obj)
        T = obj.data;
        obj.data = [];
        T.Properties.UserData = struct(obj);
    end
    
    function TT = timetable(obj,varargin)
        if istimetable(obj.data), TT = obj.data;
        else
            TT = table2timetable(obj.data,varargin{:});
        end
        obj.data = [];
        TT.Properties.UserData = struct(obj);
    end
    
    function S = struct(obj)
        
        persistent propertylists
        if isempty(propertylists), propertylists = struct(); end
        c = class(obj);
        if ~isfield(propertylists,c)
            MC = metaclass(obj);
            S = MC.PropertyList;
            S( ~strcmp({S.GetAccess},'public') | [S.NonCopyable] | ...
               [S.Dependent] | [S.Constant] | [S.Abstract] | [S.Hidden]) = [];
            propertylists.(c) = {S.Name};
        end

        if ~isempty(obj.data)
            if isa(obj.data,'table') 
                S = table2struct(obj.data,'toscalar',true);
            else
                S = timetable2struct(obj.data);
            end
            S.Properties = obj.data.Properties;
        end

        prop = propertylists.(c);
        for j = 1:numel(prop)
           S.(prop{j}) = obj.(prop{j});
        end
    end
    
    function s = size(obj,varargin), s = size(obj.data,varargin{:}); end
    function x = isempty(obj), x = builtin('isempty',obj) || isempty(obj.data); end
%     function disp(obj)
%         if isempty(obj) || isempty(obj.data), builtin('disp',obj);
%         else
%             disp(obj.data);
%         end
%     end
    function disp(obj)
        if numel(obj) == 0, builtin('disp',obj); return; end
        % name = class(obj);
        % header = [ regexprep(mat2str(obj.size),{'[',']',' '},{'','','x'}),...
        %     ' <a href = "matlab:help ' name '">' name '</a>' newline()];
        % disp(header);
        header = evalc('builtin(''disp'',obj)');
        header = strrep(header,'array with properties','object');
        disp(header);
        if ~isempty(obj) > 0 
            disp(head(obj.data))
        end
    end
    
    function x = head(obj), x = head(obj.data); end
    function x = fieldnames(obj), x = obj.data.Properties.VariableNames'; end
        % if istimetable(obj.data)
        %     x = [obj.data.Properties.DimensionNames(1);x];
        % end  
    function x = isfield(obj,f)
        x = ismember(f, obj.data.Properties.VariableNames); 
    end
    
    function x = horzcat(varargin), x = cat(2,varargin{:}); end
    function x = vertcat(varargin), x = cat(1,varargin{:}); end
    function cat(varargin)
        error('table_ish object arrays are currently forbidden');
    end
    
    function X = renamevars(X,varargin), X.data = renamevars(X.data,varargin{:}); end
end
% methods (Access=protected)
%    function propgrp = getPropertyGroups(obj)
%         propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
%         propList = struct('Department',obj.Department,...
%             'JobTitle',obj.JobTitle,...
%             'Name',obj.Name,...
%             'Salary','Not available',...
%             'Password',pd);
%         propgrp = matlab.mixin.util.PropertyGroup(propList);
%     end
%     
%     function s = getHeader(obj)
%         s = getHeader@matlab.mixin.CustomDisplay(obj);
%     end 
% end
end

function yn = hascharkeys(x)
    if iscell(x), t = x{1};
    else, t = x;
    end
    yn = ~isnumeric(t) && ~islogical(t) && ~isequal(t,':');
end