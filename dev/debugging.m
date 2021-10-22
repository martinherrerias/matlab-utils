function b = debugging()
    b = ~isempty(dbstatus) || feature('IsDebugMode');
end
