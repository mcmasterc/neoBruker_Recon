% %Split FIDs into two arrays of double precision floats with 70 rows
% below was attempt 1, after running radial_recon.m
% 42 FIDs = reshape(FIDs,[],Method_Params.NPro*Method_Params.NumTEs*Method_Params.Repetitions);%,Method_Params.Repetitions);
FIDs_coil1 = FIDs(1:70, :);
FIDs_coil2 = FIDs(71:140, :);
% Triggered example mouse data was garbage, leading my suspicions to
% believe the FIDs from each coil channel were mixed up maybe by the above
% 3 lines of code.

%% attempt 2: FIDs reshape after Bruker_load
% 
% % make two raw FID vectors, one for each channel.
% FIDs_coil1 = FIDs(1:length(FIDs)/2);
% FIDs_coil2 = FIDs((length(FIDs)/2+1):end);
% 
% % now reshape using reshape function and NPro
% FIDs_coil1 = reshape(FIDs_coil1,[],Method_Params.NPro);
% FIDs_coil2 = reshape(FIDs_coil2,[],Method_Params.NPro);

%{
plot(abs(FIDs_coil1(:,1:30)),'LineWidth',1)
title("FIDs coil 1")
legend("proj1","proj2","...")
xlabel("Readout")
ylabel("Intensity [a.u.]")
xlabel("Projection")
%}