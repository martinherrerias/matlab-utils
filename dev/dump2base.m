function dump2base()
% Copies the caller workspace to 'base', optionally clearing it beforehand.
% Intended to be used for debugging purposes

    r = evalin('caller','whos();');
    if ~isempty(r)
        switch questdlg(['The base-workspace is not empty. Dumping this workspace''s variables '...
                         'might change things in unsuspected — possibly interesting, yet '...
                         'probably unwanted — ways. What do you want to do?'],'dump2base',...
                         'See what happens','Clear base before','Stop!','Stop!')
            case 'See what happens' % Don't do a thing
            case 'Clear base before', evalin('base','clear variables');
            case {'Stop!',''}, return;
        end
    end
    filename = ['~ws_dump_' datestr(now,'yymmddhhMMss') '.mat'];
    evalin('caller',['save(''' filename ''');']);
    evalin('base',['load ' filename]);
    delete(filename);
end
