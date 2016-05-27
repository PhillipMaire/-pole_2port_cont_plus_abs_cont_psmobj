% [x, y] = SavingSection(obj, action, x, y)
%
% Section that takes care of saving/loading, etc.
%
% PARAMETERS:
% -----------
%
% obj      Default object argument.
%
% action   One of:
%            'init'      To initialise the section and set up the GUI
%                        for it
%
%            'reinit'    Delete all of this section's GUIs and data,
%                        and reinit, at the same position on the same
%                        figure as the original section GUI was placed.
%
%            'savesets'  Save GUI settings to a file
%
%            'loadsets'  Load GUI settings to a file
%
%            'savedata'  Save all SoloParamHandles to a file
%
%            'loaddata'  Load all SoloParamHandles from a file
%
%
% x, y     Only relevant to action = 'init'; they indicate the initial
%          position to place the GUI at, in the current figure window
%
% RETURNS:
% --------
%
% [x, y]   When action == 'init', returns x and y, pixel positions on
%          the current figure, updated after placing of this section's GUI. 
%

function [x, y] = SavingSection(obj, action, x, y)
   
   GetSoloFunctionArgs;
   
   switch action
    case 'init',      % ------------ CASE INIT --------------------
      % Save the figure and the position in the figure where we are
      % going to start adding GUI elements:
      SoloParamHandle(obj, 'my_gui_info', 'value', [x y gcf]);
      
%       PushButtonParam(obj, 'loaddata', x, y, 'label', 'Load Data');      
%       set_callback(loaddata, {mfilename, 'loaddata'}); next_row(y);     
%       PushButtonParam(obj, 'savedata', x, y, 'label', 'Save Data');
%       set_callback(savedata, {mfilename, 'savedata'}); next_row(y, 1.5);   
      EditParam(obj, 'MouseName', 'JF000', x, y); next_row(y, 1.5);
            
      PushButtonParam(obj, 'loadsets', x, y, 'label', 'Load Settings');
      set_callback(loadsets, {mfilename, 'loadsets'});
      next_row(y);     
      PushButtonParam(obj, 'savesets', x, y, 'label', 'Save Settings');
      set_callback(savesets, {mfilename, 'savesets'});
      next_row(y, 1.5);     
      
%       PushButtonParam(obj, 'loaddata', x, y, 'label', 'Load Data');      
%       set_callback(loaddata, {mfilename, 'loaddata'});
%       next_row(y);     
      PushButtonParam(obj, 'savedata', x, y, 'label', 'Save Data');
      set_callback(savedata, {mfilename, 'savedata'});
      next_row(y, 1.5);     

%       SubheaderParam(obj, 'sectiontitle', 'Saving/Loading', x, y);
%       next_row(y, 1.5);
% 
      SoloParamHandle(obj, 'SaveTime');       
      SoloParamHandle(obj, 'hostname', 'value', get_hostname);
      
      
      
      

      EditParam(obj, 'Weight', 'Not recorded', x, y); next_row(y);
      EditParam(obj, 'WeightAfterExp', 'Not recorded', x, y); next_row(y);
      MenuParam(obj, 'HeadFixed', {'Fixed w/o anesthesia','Fixed after anesthesia','Not fixed'},'Fixed after anesthesia', x, y); next_row(y);
      MenuParam(obj, 'Lighting', {'IR Security Cam','Diffuse 940nm for Security Cam','940nm High-speed Imaging','White light','Totally Dark'},'Diffuse 940nm for Security Cam', x, y); next_row(y);
      MenuParam(obj, 'Punishment', {'None','ExtraITI','60 psi',...
          '55 psi','50 psi','45 psi','40 psi','35 psi','30 psi','25 psi',...
          '20 psi','15 psi','10 psi','Quinine'},'ExtraITI', x, y); next_row(y);  
      MenuParam(obj, 'WhiskerSet', {'Full','C2 only','None','Row C only'},'Full', x, y); 
      next_row(y, 1);
      
      % Record which rig experiment was run on:
      [status,host] = system('hostname');
      SoloParamHandle(obj, 'RigComputer', 'value', host);
      SoloParamHandle(obj, 'SaveTime');            
%       PushButtonParam(obj, 'loadsets', x, y, 'label', 'Load Settings');
%       set_callback(loadsets, {mfilename, 'loadsets'});
%       next_row(y);     
%       PushButtonParam(obj, 'savesets', x, y, 'label', 'Save Settings');
%       set_callback(savesets, {mfilename, 'savesets'});
%       next_row(y, 1.5);     
      
  

      SubheaderParam(obj, 'sectiontitle', 'Documentation', x, y);
      next_row(y, 1.5);
     
      return;
      
    case 'savesets',       % ------------ CASE SAVESETS --------------------
        SaveTime.value = datestr(now);
        save_solouiparamvalues(value(MouseName));
        disp(['Settings saved at ' value(SaveTime) ' for ' value(MouseName)]);
      return;
      
    case 'loadsets',       % ------------ CASE LOADSETS --------------------
      load_solouiparamvalues(value(MouseName));
      
      return;
      
    case 'autosave',       % ------------ CASE AUTOSAVE --------------------
      sp_autosave(value(MouseName));
      
      return;

      
    case 'savedata',       % ------------ CASE SAVEDATA --------------------
        SaveTime.value = datestr(now);
        % Get string from annotation object, put in a SoloParamHandle (notes) which will be saved:       
%      notes.value = get(value(noteshandle),'String');

      save_soloparamvalues(value(MouseName));
        disp(['Data saved at ' value(SaveTime) '  for  ' value(MouseName)]);
      return;
      
%     case 'loaddata',       % ------------ CASE LOADDATA --------------------
%       load_soloparamvalues(value(MouseName));
%       
%       return;
    
    
    case 'reinit',       % ------------ CASE REINIT --------------------
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
   
   
      