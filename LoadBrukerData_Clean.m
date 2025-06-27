function LoadBrukerData_Clean(path)

%%%%% Temp hardcode 
path = "\\rds6.cchmc.org\PulMed-43\CPIR_Share\Carter\01_For Samal_19Jun2025\PreclinicalRadialRecon_TestNewScanner\NewScanner\New_Ungated";
%%%%%
cd(path)
fid_File = dir('*.job0');

%% IDK if this is necessary or not for the other bits of the code library yet --CBM
% Read in Configuration File
%If it doesn't exist, create one
ConfigPresent = exist(fullfile(path,'ConfigFile_V4.txt'),'file');

if ConfigPresent == 0
    gen_config_file_v4(path)
end


Config_Params = read_config_file_v4(path);

%% Read Method File
disp('Reading Method File')

[traj,Method_Params] = read_method(path);


%If we are doing a phantom, want to pretend that we are triggered so that
%we don't try to retrogate (which would be pointless)
if strcmpi(Config_Params.Subject,'Phantom') %non case sensitive string compare
    Method_Params.Trigger = 'On';
end
%% Load FIDs
disp('Reading FID File')
fid_File = dir('*.job0'); 
FIDs = Bruker_Load(fid_File.name);
%% Divide up the raw data vector by nPro and coil
FIDs = reshape(FIDs,[],Method_Params.NPro*Method_Params.NumTEs*Method_Params.Repetitions);
FIDs_coil1 = FIDs(1:70, :);
FIDs_coil2 = FIDs(71:140, :);
if size(traj,2) > size(FIDs_coil1,1)
    traj = traj(:,1:size(FIDs_coil1,1),:);
end
