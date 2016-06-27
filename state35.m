% [] = state35(obj)    Method that gets called at the end of every trial
%
% This method assumes that it has been given read/write access to a SPH
% called n_done_trials (which will be updated by one upon entry to
% state35.m), and read-only access to a SPH called
% trial_finished_actions. This last should be a cell vector of strings,
% each of which will be eval'd in sequence (type "help eval" at the
% Matlab prompt if you don't know what that means).
%
% If you put everything into trial_finished_actions, the code for this
% method should be universal for all protocols, and there should be no
% need for you to modify this file.
%

% CDB Feb 06


function [] = state35(obj)

GetSoloFunctionArgs;
% SoloFunctionAddVars('state35', 'rw_args', 'n_done_trials', ...
%         'ro_args', 'trial_finished_actions');


SoloParamhandle(obj, 'iti', 'value', 0); % for debugging motor move time
tic

n_done_trials.value = n_done_trials + 1;

for i=1:length(trial_finished_actions),
    eval(trial_finished_actions{i});
end;

n_started_trials.value = n_started_trials + 1;



if ismember(SessionType, {'Discrim_DHO','2port-Discrim','Manual-Training'})
        disp(['Starting move...'])     
%         tic
        MotorsSection(obj,'move_next_side'); % chooses and *moves* to next side
%         error('TESTING!')
        disp(['Move finished. Ready for trial ' int2str(value(n_started_trials))]) 
%         toc        
else
    % nothing
end

iti.value = toc;




