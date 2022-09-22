function startup = startExperiment(experimentName, isPractice)

disp('working from directory')
disp([pwd ' ...'])
disp(' ')

if IsWin
    dataDir = [pwd '\data\'];
else
    dataDir = [pwd '/data/'];
end

startupTime = clock;

if isPractice
    currentFile = [experimentName ' PRACTICE' datestr(startupTime,' yyyy-mm-dd HH-MM-SS') '.mat'];
else
    currentFile = [experimentName datestr(startupTime,' yyyy-mm-dd HH-MM-SS') '.mat'];
end
dataFile    = [dataDir currentFile];

startup.dataDir        = dataDir;
startup.dataFile       = dataFile;
startup.startupTime    = startupTime;
startup.startupTimeStr = datestr(startupTime);