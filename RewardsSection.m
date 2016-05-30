 function [] = RewardsSection(obj, action)
   global next_side %made 'next_side' global!!!-psm
   GetSoloFunctionArgs;
   
   switch action
    case 'init',
      SoloParamHandle(obj, 'pstruct');
      SoloParamHandle(obj, 'LastTrialEvents', 'value', []);
      
    case 'trial_finished',
      % Make sure we've collected up to the latest events from the RT machine:  
      RewardsSection(obj, 'update');
      
      % Parse the events from the last trial:
      pstruct.value = parse_trial(value(LastTrialEvents), RealTimeStates);
      
%       % Parse the events from the last trial:
%       pstruct.value = parse_trial(value(LastTrialEvents), RealTimeStates,...
%           {'statename_list', 'all'}, {'parse_pokes', false});
      
      % Take the current raw events and push them into the history:
      push_history(LastTrialEvents); LastTrialEvents.value = [];
      
      if rows(pstruct.punish)>0 
          hit = 0;
          %this will work fine for the trial 'a' condition, just have to plot based on 
          %lick port and trial type, trial type 'a' then lick of wither port is an error
          %thus plot red and also distingish L or R port.
      elseif rows(pstruct.miss)>0 %miss just means mouse didn't lick, thats it. 
          %So for abscent correct rejection, must define this based on the trial type. 
          if next_side=='a'
              hit = 2; %inicating a correct rejection-psm
          else
          hit = -1;
          end
      else
          hit = 1;
      end
      
      
      hh = hit_history(:);
      hh(n_started_trials) = hit;
      if (size(hh,1) == 1) ; hh = hh'; end
      hit_history.value = [hh];
      
      
    case 'update',
      Event = GetParam('rpbox', 'event', 'user');
      LastTrialEvents.value = [value(LastTrialEvents) ; Event];
        
    case 'reinit',

      % Delete all SoloParamHandles who belong to this object and whose
      % fullname starts with the name of this mfile:
      delete_sphandle('owner', ['^@' class(obj) '$'], ...
                      'fullname', ['^' mfilename]);

      % Reinitialise 
      feval(mfilename, obj, 'init');
   end;
   
   
      