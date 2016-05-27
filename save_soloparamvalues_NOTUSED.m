% [] = sp_autosave(ratname, varargin)
%
% Loosely moeled on save_solouiparamvalues; this will 
%  1) obtain or create a session_id stored in global autosave_session_id
%  2) use user-specified animal name to, if needbe make a directory
%  3) in that directory, it will
%     a) create a file {animalid}_{sessionid}.mat, saving ALL trials to it
%     b) OR copy above file to {animalid}_{sessionid}.mat_last and then a)
%
%
% PARAMETERS:
% ----------
%
% ratname     This will determine which directory the file goes into.
%


function [] = sp_autosave_NOTUSED(ratname, varargin)

pairs = { ...
    'child_protocol', [] ; ...
    'asv', 0; ...
    'interactive'      1 ; ...
    'commit'           0 ; ...
    'remove_asvs'      1 ; ...
    };
parse_knownargs(varargin, pairs);

if isempty(child_protocol),
    owner = determine_owner;
else
    owner = class(value(child_protocol));   % the child protocol owns all vars
   % owner = c(2:end);
end;

% Don't do anything if user is retarted
if (strcmp(ratname, 'JF000'))
    disp(['sp_autosave::you must set your mousename to non-default JF000 for autosave to work.']);
    return;
end

% Create session ID
if (length(whos('global','autosave_session_id')) == 0)
    global autosave_session_id;
    autosave_session_id = datestr(now,30);
else
    global autosave_session_id;
end
    
% Did I mention I hate most solo code because there is almost NO comments?

% --- prelims -- access workspace, figure datapath otu
   owner = owner;
   global Solo_datadir;
   if isempty(Solo_datadir), 
      Solo_datadir=[pwd filesep '..' filesep 'SoloData'];
   end;
   data_path = [Solo_datadir filesep 'Data'];
   if ~exist(data_path, 'dir'),
      success = mkdir(Solo_datadir, 'Data');
      if ~success, error(['Couldn''t make directory ' data_path]); end;
   end;
   
   handles = get_sphandle('owner', owner);
   k = zeros(size(handles));
   for i=1:length(handles),
      if get_saveable(handles{i}), k(i) = 1; end;
   end;
   handles = handles(find(k));
   
   % --- start building the written arrays
   saved = struct; saved_history = struct; saved_autoset = struct;
   for i=1:length(handles),
      saved.(get_fullname(handles{i}))        = value(handles{i});
      saved_history.(get_fullname(handles{i}))= get_history(handles{i});
      saved_autoset.(get_fullname(handles{i}))= get_autoset_string(handles{i});
   end;

   
   protocol_name = get_sphandle('owner', owner, 'name', 'protocol_name');
   if ~isempty(protocol_name),
      protocol_name = protocol_name{1};   
      fig_position = get(...
        findobj(get(0, 'Children'), 'Name', value(protocol_name)),'Position');
   else
      fig_position = [];
   end;
   
   % Now set the path for the data file
   
   if data_path(end)~=filesep, data_path=[data_path filesep]; end;
   
   rat_dir = [data_path ratname];
    if ~exist(rat_dir)
    success = mkdir(data_path, ratname);
    if ~success, error(['Couldn''t make directory ' rat_dir]); end;
end;
if rat_dir(end)~=filesep, rat_dir=[rat_dir filesep]; end;
    
% --- ACTUAL NAMES:
    pname = rat_dir;
    fname = [ 'data_' ratname '_' ...
        autosave_session_id '.mat'];

    % --- copy to _last if needbe -- in the unlikely event that you 
    %     experience a crash, one of these is guaranteed to survive
    if (exist([pname fname], 'file') ~= 0)
        copyfile([pname fname], [pname fname '_last']);
    end
    
    % --- The WRITE:
    save([pname fname], 'saved', 'saved_history', 'saved_autoset', ...
        'fig_position');

    disp([datestr(now) '::Autosaved to: ' pname fname]);
