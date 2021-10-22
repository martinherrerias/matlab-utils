function varargout = dispnested(S)
% DISPNESTED(S) - display contents of nested structure, indenting and aligning semicolons at all
%   hierarchical levels.

    TABSIZE = 4;

    f = fieldnames(S);
    substructures = cellfun(@(x) isstruct(x) && isscalar(x),struct2cell(S));
    
    txt = evalc('disp(S)');
    txt = strsplit(txt,newline());
    txt(end) = []; % empty end line
    
    pivot = strfind(txt{1},':');
    pivot = pivot(1);
    
    if any(substructures)
        for j = find(substructures)'
            txt{j}(pivot+1:end) = []; % remove ' [1×1 struct]' after field:
            [subtxt,subpivot] = dispnested(S.(f{j}));
            if (pivot - subpivot)+4 < 0
                txt = pad(txt,(subpivot-pivot)-TABSIZE);
                pivot = subpivot-TABSIZE;
            else
                subtxt = pad(subtxt,(pivot - subpivot)+TABSIZE);
            end
            txt{j} = [txt{j},newline(),subtxt];
        end
    end
    txt = strjoin(txt,newline());

    if nargout > 0
        varargout = {txt,pivot};
    else
        fprintf('%s\n',txt);
    end
    
    function txt = pad(txt,n)
        padding = repmat(' ',1,n);
        if iscell(txt)
            txt = cellfun(@(l) [padding,l],txt,'unif',0);
        else
            txt = [padding,strrep(txt,newline(),[newline(),padding])];
        end
    end
end