
cd '/home/martin/solar-simulation/simulation'
codefiles = arrayfun(@(d) fullfile(d.folder,d.name),dir('**/*.m'),'unif',0);

reqfiles = {'splitscript.m','splitGUI.m','samplesystem.m','GUIsolver.m','GUIterrain.m'...
    'GUIarrdef.m','GUIfigure.m','GUIirrtrans.m','GUIlayout.m','GUImeteo.m',...
    'GUImodels.m','GUIphysmod.m','GUIphystrck.m','GUIreduce.m','GUIsetup.m','GUIshading.m'}';

reqfiles = cellfun(@(f) matlab.codetools.requiredFilesAndProducts(f),reqfiles,'unif',0);
reqfiles = uniquecell(cat(2,reqfiles{:})');
    
relativepath(setdiff(codefiles,reqfiles))
