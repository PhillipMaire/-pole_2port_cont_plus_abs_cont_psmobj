% [x, y] = NotesSection(obj, action, x, y)
%
% Text box for taking notes during the experiment.
%
%
% PARAMETERS:
% -----------
%
% obj      Default object argument.
%
% action   One of:
%            'init'      To initialise the section and set up the GUI
%                        for it;
%
%            'reinit'    Delete all of this section's GUIs and data,
%                        and reinit, at the same position on the same
%                        figure as the original section GUI was placed.
%           
%            Several other actions are available (see code of this file).
%
% x, y     Relevant to action = 'init'; they indicate the initial
%          position to place the GUI at, in the current figure window
%
% RETURNS:
% --------
%
% [x, y]   When action == 'init', returns x and y, pixel positions on
%          the current figure, updated after placing of this section's GUI.
%
% x        When action = 'get_next_side', x will be either 'l' for
%          left or 'r' for right.
%

function [x, y] = NotesSection(obj, action, x, y)

GetSoloFunctionArgs;


switch action

    case 'init',   % ------------ CASE INIT ----------------
        
        % Save the figure and the position in the figure where we are
        % going to start adding GUI elements:
        SoloParamHandle(obj, 'my_gui_info', 'value', [x y gcf]); 
        next_row(y);

        MenuParam(obj, 'notes_show', {'view', 'hide'}, 'hide', x, y, 'label', 'Notes', 'TooltipString', 'Experiment notes');
        set_callback(notes_show, {mfilename,'hide_show'});

        next_row(y,1);
        SubheaderParam(obj, 'sectiontitle', 'Notes', x, y);

        parentfig_x = x; parentfig_y = y;

        SoloParamHandle(obj, 'notesfig', 'saveable', 0);
        notesfig.value = figure('Position', [3 750 200 400], 'Menubar', 'none',...
            'Toolbar', 'none','Name','Notes','NumberTitle','off');
        
        SoloParamHandle(obj, 'noteshandle');
        noteshandle.value = annotation(value(notesfig),'textbox',[.1 .1 .8 .8]);
        plotedit(value(notesfig));
        
        SoloParamHandle(obj, 'notes');
                
        SoloFunctionAddVars('SavingSection', 'rw_args', {'notes', 'noteshandle'});
        
        x = 1; y = 1;
        
        NotesSection(obj,'hide_show');
        
        x = parentfig_x; y = parentfig_y;
        set(0,'CurrentFigure',value(myfig));
        return;


    case 'hide_show'
        if strcmpi(value(notes_show), 'hide')
            set(value(notesfig), 'Visible', 'off');
        elseif strcmpi(value(notes_show),'view')
            set(value(notesfig),'Visible','on');
        end;
        return;


    case 'reinit',   % ------- CASE REINIT -------------
        currfig = gcf;

        % Get the original GUI position and figure:
        x = my_gui_info(1); y = my_gui_info(2); figure(my_gui_info(3));

        delete(value(myaxes));

        % Delete all SoloParamHandles who belong to this object and whose
        % fullname starts with the name of this mfile:
        delete_sphandle('owner', ['^@' class(obj) '$'], ...
            'fullname', ['^' mfilename]);

        % Reinitialise at the original GUI position and figure:
        [x, y] = feval(mfilename, obj, 'init', x, y);

        % Restore the current figure:
        figure(currfig);
        return;
end


