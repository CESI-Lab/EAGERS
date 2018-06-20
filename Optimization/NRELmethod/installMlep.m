% INSTALLMLEP code to install MLE+
%      Run this script before using MLE+.
%
%      Use: installMlep
%
%      In installMlep you need to specify whether you want to use the GUI
%      mode or the Manual mode. Set manualInstall = 0 if you do not want to
%      use the GUI for installation.
%      GUI (manualInstall = 1): A installation screen will pop up according
%      to your operating system (PC or UNIX)
%      MANUAL (manualInstall = 0): You need to specify the E+ and Java
%      directory for Windows machines and only the E+ directory for UNIX
%      systems.
%
%      To save your configurations permanently in Matlab, you need to have
%      admin privileges. If you get the following: "Warning: Unable to save
%      apth to file ..." You either need to run this script everytime you
%      open Matlab or pass a path to the savepath function at the end of
%      the script.
%
% Last Modified by Willy Bernal - Willy.BernalHeredia@nrel.gov 30-Jul-2015

%% MODIFY
% true = Install Manually
% fase = Install through GUI
manualInstall = true;

% Paths
if ispc
    % Windows
    eplusPath = 'C:\EnergyPlusV8-8-0';
%     javaPath = 'C:\Program Files\Java\jre1.8.0_171\bin';
    javaPath = 'C:\Program Files (x86)\Java\jre1.8.0_171\bin';
else
    % Unix
    eplusPath = '/Applications/EnergyPlus-8-4-0';
end

%% DO NOT MODIFY
% Extract MLE+ Main Path
filename = mfilename('fullpath');
[mlepPath, ~, ~] = fileparts(filename);

% Add Path
% addpath(mlepPath);
addpath(genpath(fullfile(mlepPath,'MLE+','core')));
% addpath(fullfile(mlepPath,'gui'));
addpath(genpath(fullfile(mlepPath,'MLE+','install')));
%addpath(fullfile(mlepPath,'settings'));
addpath(genpath(fullfile(mlepPath,'MLE+','version')));
% addpath(fullfile(mlepPath,'mlepHelp'));
addpath(genpath(fullfile(mlepPath,'MLE+','library')));
bcvtbPath = fullfile(mlepPath,'MLE+','bcvtb');
addpath(genpath(bcvtbPath));

% Installation Dialog
if ~manualInstall %========================================================
    if ispc
        % WIN
        installationWin;
    else
        % UNIX
        installationUnix;
    end
else %=====================================================================
    % Manual Install
    if ispc
        % No Changes
    else
        % Java Path
        javaPath = '';
    end
    
    % Save Settings =======================================================
    mlepSaveSettings(mlepPath, eplusPath, javaPath, bcvtbPath);
    %%
    disp('================================');
    disp('        MLE+ COSIMULATION       ');
    disp('     INSTALLATION COMPLETED     ');
    disp(' Install everytime Matlab opens.');
    disp('================================');
end

% Save Paths
% savepath;
savepath ./pathdef.m


