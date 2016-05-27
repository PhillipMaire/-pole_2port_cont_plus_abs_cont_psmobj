function [x, y] = TimesSection(obj, action, x, y)

GetSoloFunctionArgs;

switch action
    case 'init',
        % Save the figure and the position in the figure where we are
        % going to start adding GUI elements:
        SoloParamHandle(obj, 'my_gui_info', 'value', [x y gcf]);
        
        %       NumeditParam(obj, 'ExtraITIOnError', 2, x, y);
        NumeditParam(obj, 'PostTrialPauseTime', 2.5, x, y, 'TooltipString', ...
            'Time after everything -- ensures minimal desired trial length');
        next_row(y);
        NumeditParam(obj, 'RewardCollectTime', 1, x, y, 'TooltipString', ...
            'Time after animal gets reward where lickoprt not withdrawn');
        next_row(y);        MenuParam(obj, 'ExtraITIOnError', {'0', '2','3','4','5','6','7','8','9','10','11','12','13','14','15','20','25','30'},'2', x, y);
        next_row(y);
        NumeditParam(obj, 'AnswerPeriodTime', 1, x,y); % Can be modified by make_and_upload_state_matrix
        next_row(y, 1);
        EditParam(obj, 'LickportTravelTime', '.001', x, y, 'TooltipString', ...
            'How long for lickport to move into reach? NOT part of trial duration');
        next_row(y);
        MenuParam(obj, 'RestartPreAnsOnLick', {'on','off'},'off', x, y);      
        next_row(y);          
        EditParam(obj, 'PreAnswerTime', '.5', x, y, 'TooltipString', ...
            'Time animal waits after pole is retracted but before answer counts');
        next_row(y);
        MenuParam(obj, 'SamplingPeriodTime', {'0.2','0.5','0.75','1','1.25','1.5','1.75','2','2.5','3'},'1', x, y);
        next_row(y);
        NumeditParam(obj, 'PoleRetractTime', 0.2, x, y, 'TooltipString', ...
            'Time allowed for removing the pole before start answering period');
        next_row(y);
        NumeditParam(obj, 'PreTrialPauseTime', 1, x, y, 'TooltipString', ...
            'Time before anything happens -- basically to allow F_0 sampling');
        next_row(y);

%         if value(SamplingPeriodTime)>1
%             aptm = 1;
%         else
%             aptm = 2 - value(SamplingPeriodTime);
%         end

        SoloFunctionAddVars('make_and_upload_state_matrix', 'ro_args', ...
            {'ExtraITIOnError','SamplingPeriodTime','PreTrialPauseTime', ...
            'PreAnswerTime','PostTrialPauseTime','PoleRetractTime','RewardCollectTime', ...
            'LickportTravelTime','RestartPreAnsOnLick'});
        SoloFunctionAddVars('make_and_upload_state_matrix', 'rw_args', ...
            {'AnswerPeriodTime'});
            
        SubheaderParam(obj, 'title', 'Times', x, y);
        next_row(y, 1.5);
        
    case 'reinit',
        currfig = gcf;
        
        % Get the original GUI position and figure:
        x = my_gui_info(1); y = my_gui_info(2); figure(my_gui_info(3));
        
        % Delete all SoloParamHandles who belong to this object and whose
        % fullname starts with the name of this mfile:
        delete_sphandle('owner', ['^@' class(obj) '$'], ...
            'fullname', ['^' mfilename]);
        
        % Reinitialise at the original GUI position and figure:
        [x, y] = feval(mfilename, obj, 'init', x, y);
        
        % Restore the current figure:
        figure(currfig);
end;

   
      