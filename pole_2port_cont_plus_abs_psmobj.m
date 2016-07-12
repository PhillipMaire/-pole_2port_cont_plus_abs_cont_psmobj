function [obj] = pole_2port_cont_plus_abs_psmobj(name)
    
% -------- BEGIN Magic code that all protocol objects must have ---
%
% No need to alter this code: will work for every protocol. Jump to 
% "END Magic code" to start your own code part.
%
   
% Default object:
obj = class(struct, mfilename);

% If creating an empty object, return without further ado:
if nargin==0 || (nargin==1 && ischar(name) && strcmp(name, 'empty')), 
   return; 
end;

% Delete previous vars owned by this object:
delete_sphandle('owner', ['^@', mfilename '$']); 

% Non-empty: proceed with regular init of this object
if nargin~=1 || ~isstr(name), 
   error(['To initialize this protocol, you need to call the constructor ' ...
          'with the name of the protocol as its only arg.']); 
end;

% Make default figure. We remember to make it non-saveable; on next run
% the handle to this figure might be different, and we don't want to
% overwrite it when someone does load_data and some old value of the
% fig handle was stored as SoloParamHandle "myfig"
SoloParamHandle(obj, 'myfig', 'saveable', 0); myfig.value = figure;

% Give close.m access to myfig, so that upon closure the figure may be
% deleted: 
SoloFunctionAddVars('close', 'ro_args', 'myfig');

% Give MotorsSection.m, SMControlSection.m access to myfig, so that they can switch active fig back
% upon opening/closing their own windows:
SoloFunctionAddVars('MotorsSection', 'ro_args', 'myfig');
SoloFunctionAddVars('SMControlSection', 'ro_args', 'myfig');
SoloFunctionAddVars('AnalysisSection', 'ro_args', 'myfig');
SoloFunctionAddVars('NotesSection', 'ro_args', 'myfig');

% Make the title of the figure be the protocol name, and if someone tries
% to close this figure, call Exper's ModuleClose function, so it'll know
% to take it off Exper's list of open protocols.
name = lower(name);
set(value(myfig), 'Name', name, 'Tag', name, ...
                  'closerequestfcn', ['ModuleClose(''' name ''')'], ...
                  'NumberTitle', 'off', 'MenuBar', 'none');

%
% -------- END Magic code that all protocol objects must have ---
%



% Ok, at this point we have one SoloParamHandle, myfig 

% Let's put the figure where we want it and give it a reasonable size:
set(value(myfig), 'Position', [485   100   450   780]);

% Let's declare some globals that everybody is likely to want to know about.
% Number of finished trials:
SoloParamHandle(obj, 'n_done_trials', 'value', 0);

% Number of trials started (*during* a trial, n_started_trials will be
% one greater than n_done_trials; during an intertrial interval, the
% two will be the same):
SoloParamHandle(obj, 'n_started_trials', 'value', 0);

% History of hit/miss:
SoloParamHandle(obj, 'hit_history', 'value', []);

% History of EXACT zaber pole positions:
SoloParamHandle(obj, 'pole_position_history', 'value', []);


% History of EXACT lickport positions:
SoloParamHandle(obj, 'lickport_position_history', 'value', []);

% History of left port probabilities, AFTER bias correction
SoloParamHandle(obj, 'left_port_prob_history', 'value', []);

% History of right port probabilities, not yet adjusted to work with bias correction-psm
SoloParamHandle(obj, 'right_port_prob_history', 'value', []);

% Position of poles for each trial:%%% I DONT USE THIS-SP
npos=5; SoloParamHandle(obj, 'pole_position_list', 'value', cell(npos,1));



% Every function will be able to read these, but only those explicitly
% given r/w access will be able to modify them:
DeclareGlobals(obj, 'ro_args', {'n_done_trials', 'n_started_trials', ...
                    'hit_history', 'pole_position_history', 'pole_position_list',...
                    'lickport_position_history','left_port_prob_history'});

% Let RewardsSection, the part that parses what happened at the end of
% a trial, write to hit_history:
SoloFunctionAddVars('RewardsSection', 'rw_args', 'hit_history');

% Let MotorsSection access pole_position_history
SoloFunctionAddVars('MotorsSection', 'rw_args', 'pole_position_history');
SoloFunctionAddVars('MotorsSection', 'rw_args', 'lickport_position_history');

% Sides section gets left ort prob
SoloFunctionAddVars('MotorsSection', 'rw_args', 'left_port_prob_history');
% ----------

x = 1; y = 1;                     % Initial position on main GUI window

% Section of the code that takes care of saving/loading data and
% settings. This section's GUIs will be placed at x, y, and upon returning
% from its initialization, it will have updated x and y so that the next
% section's GUIs doesn't go on top.
[x, y] = SavingSection(obj, 'init', x, y); 

% next_row(y);
[x, y] = SessionTypeSection(obj, 'init', x, y);

[x, y] = NotesSection(obj, 'init', x, y);

% Section that allows interaction with the RTLinux state machine.
[x, y] = SMControlSection(obj, 'init', x, y);

% Online analysis:
[x, y] = AnalysisSection(obj, 'init', x, y);

% Control of motors:
[x, y] = MotorsSection(obj, 'init', x, y);

next_column(x); y = 1;
% Section that sets ITI times, etc. 
[x, y] = TimesSection(obj, 'init', x, y);

% next_column(x); y = 1;
% % % Section that sets the sound for each side. 
% [x, y] = SoundSection(obj, 'init', x, y);

% Section of the code that takes care of how long to keep the water and
% pneumatic slider valves open:
[x, y] = ValvesSection(obj, 'init', x, y);  

% Section that chooses the next correct side and keeps track of previous
% side history. 
[x, y] = SidesSection(obj, 'init', x, y);



figpos = get(gcf, 'Position');
HeaderParam(obj, 'prot_title', 'Whisker Discrimination', ...
            x, y, 'position', [10 figpos(4)-25, 400 20]);


% Initialize section that takes care of preparing and uploading the state
% matrix, and pload up the first trial's state matrix.
make_and_upload_state_matrix(obj, 'init');
% The first trial starts upon uploading the first state matrix:
n_started_trials.value = 1;


% Initialize section that parses the results of each trial, and keeps
% track of hits and errors. 
RewardsSection(obj, 'init');



% ------------------------------------------------------------------
% List of functions to call, in sequence, when a trial is finished:
% If adding a function to this list,
%    (a) Declare its args with a SoloFunctionAddVars() call
%    (b) Add your function as a method of the current object
%    (c) As the first action of your method, call GetSoloFunctionArgs;
%
% Below, push_history(class(obj)) will push the history of all GUI
% SoloParamHandles; 
SoloParamHandle(obj, 'trial_finished_actions', 'value', { ...
  'RewardsSection(obj, ''trial_finished'');'             ; ...
  'ComputeAutoSet'                                       ; ...
  'SidesSection(obj, ''choose_next_side'');'             ; ...
  'SidesSection(obj, ''update_plot'');'                  ; ...
  'AnalysisSection(obj, ''update'');'                    ; ... 
%   'MotorsSection(obj, ''move_next_side'');'                ; ... 
  'make_and_upload_state_matrix(obj, ''next_matrix'');'  ; ... 
  'push_history(class(obj));'                            ; ... % no args
});
% Here's a breakdown of what we did above:
% 1) parse the result of the latest trial, update hit_history
% 2) go through any autoset strings defined for GUI parameters
% 3) Choose which side is going to be the next correct side:
% 4) update the plot of sides.
% 5) upload new sounds if necessary
% 6) prepare and send the next state matrix

% state 35 will add one to n_done_trials, then go through the
% trial_finished actions, then add one to n_started_trials
SoloFunctionAddVars('state35', ...
    'rw_args', {'n_done_trials', 'n_started_trials'}, ...
    'ro_args', 'trial_finished_actions');

% List of functions to call, in sequence, when an update call is made.
% (Every 350 ms or so, the RPbox.m module polls the RT machine to ask for
% the latest events, and then calls the update.m method on the
% current protocol.)
% This has a similar structure to trial_finished_actions. Here we call 
% update on RewardsSection to ask it to accumulate the events for the
% trial.
SoloParamHandle(obj, 'within_trial_update_actions', 'value', { ...
    'RewardsSection(obj, ''update'');'                    ; ...
    'MotorsSection(obj, ''update'');'                    ; ...
    });
SoloFunctionAddVars('update', 'ro_args', 'within_trial_update_actions');

% ------------------------------------------------------------------

return;


