function struct2csv(S,N,filename)
% STRUCT2CSV(S,N,filename)
% Take all fields of structure (or cell array of structures) S which have size N 
% (scalar or two-vector) and copy them as columns in a csv file called filename.

	B = cell(numel(S),1);
	if isscalar(N), N(2) = 1; end
    
    % if S is a single structure, turn it into a cell-array with 1 structure
    if isstruct(S), trash = S; S = cell(1); S{1} = trash; clear trash; end
    
    fileID = fopen(filename,'w');
    
    % for all struc
    for k = 1:numel(S)
		fields = fieldnames(S{k});
		A = zeros(N(1),N(2)*numel(fields));
		used = 0;
		for j = 1:numel(fields)
            % if the field is numeric and matches size
			if ismatrix(S{k}.(fields{j})) && size(S{k}.(fields{j}),1)==N(1) &&...
                    size(S{k}.(fields{j}),2)==N(2)&& isnumeric(S{k}.(fields{j})(1))
                for l = 1:N(2), fprintf(fileID,'%s',[fields{j} ';']); end
				A(:,used+1:used+N(2)) = S{k}.(fields{j});
				used = used+N(2);
			end
		end
		B{k} = A(:,1:used);
    end
    fprintf(fileID,'\n\r');
    fclose(fileID);
    
    if isempty(cat(2,B{:}))
        warning('struct2xls:empty','No numeric fields seem to match specified size');
    end
	dlmwrite(filename,cat(2,B{:}),'-append','Delimiter',';');
end
