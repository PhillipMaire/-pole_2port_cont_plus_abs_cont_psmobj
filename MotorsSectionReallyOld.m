% [x, y] = MotorsSection(obj, action, x, y)
%
% Section that takes care of controlling the stepper motors.
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

function [x, y] = MotorsSection(obj, action, x, y)

GetSoloFunctionArgs;

global Solo_Try_Catch_Flag
global pole_motors_properties;
global lickport_motors_properties;
global pole_motors;
global lickport_motors; 

switch action

    case 'init',   % ------------ CASE INIT ----------------
        
        if strcmp(pole_motors_properties.type,'@FakeZaberTCD1000')
            pole_motors = FakeZaberTCD1000;
        else
            pole_motors = ZaberTCD1000(pole_motors_properties.port);
        end
        serial_open(pole_motors);

        if strcmp(lickport_motors_properties.type,'@FakeZaberTCD1000')
            lickport_motors = FakeZaberTCD1000;
        else
            lickport_motors = ZaberTCD1000(lickport_motors_properties.port);
        end
        serial_open(lickport_motors);

        % Save the figure and the position in the figure where we are
        % going to start adding GUI elements:
        SoloParamHandle(obj, 'my_gui_info', 'value', [x y gcf]); next_row(y,1.5);

        SoloParamHandle(obj, 'trial_ready_times', 'value', 0);

        MenuParam(obj, 'motor_show', {'view', 'hide'}, 'view', x, y, 'label', 'Motor Control', 'TooltipString', 'Control motors');
        set_callback(motor_show, {mfilename,'hide_show'});

        next_row(y);
        SubheaderParam(obj, 'sectiontitle', 'Motor Control', x, y);

        parentfig_x = x; parentfig_y = y;
       
        
        % ---  Make new window for motor configuration
        SoloParamHandle(obj, 'motorfig', 'saveable', 0);
        motorfig.value = figure('Position', [3 500 400 450], 'Menubar', 'none',...
            'Toolbar', 'none','Name','Motor Control','NumberTitle','off');

        x = 1; y = 1;

        %       PushButtonParam(obj, 'serial_open', x, y, 'label', 'Open serial port');
        %       set_callback(serial_open, {mfilename, 'serial_open'});
        %       next_row(y);

        PushButtonParam(obj, 'pole_serial_reset', x, y, 'label', 'Reset serial port connection');
        set_callback(pole_serial_reset, {mfilename, 'pole_serial_reset'});
        next_row(y);

%         PushButtonParam(obj, 'reset_motors_firmware', x, y, 'label', 'Reset Zaber firmware parameters',...
%             'TooltipString','Target acceleration, target speed, and microsteps/step');
%         set_callback(reset_motors_firmware, {mfilename, 'reset_motors_firmware'});
%         next_row(y);

        PushButtonParam(obj, 'pole_motors_home', x, y, 'label', 'Home Pole motors');
        set_callback(pole_motors_home, {mfilename, 'pole_motors_home'});
        next_row(y);

        PushButtonParam(obj, 'pole_motors_stop', x, y, 'label', 'Stop Pole motors');
        set_callback(pole_motors_stop, {mfilename, 'pole_motors_stop'});
        next_row(y);

        PushButtonParam(obj, 'pole_motors_reset', x, y, 'label', 'Reset Pole motors');
        set_callback(pole_motors_reset, {mfilename, 'pole_motors_reset'}); 
        
        next_row(y,6);       
        PushButtonParam(obj, 'lickport_serial_reset', x, y, 'label', 'Reset serial port connection');
        set_callback(lickport_serial_reset, {mfilename, 'lickport_serial_reset'});
        next_row(y);

        PushButtonParam(obj, 'lickport_motors_home', x, y, 'label', 'Home lickport motors');
        set_callback(lickport_motors_home, {mfilename, 'lickport_motors_home'});
        next_row(y);

        PushButtonParam(obj, 'lickport_motors_stop', x, y, 'label', 'Stop lickport motors');
        set_callback(lickport_motors_stop, {mfilename, 'lickport_motors_stop'});
        next_row(y);

        PushButtonParam(obj, 'lickport_motors_reset', x, y, 'label', 'Reset lickport motors');
        set_callback(lickport_motors_reset, {mfilename, 'lickport_motors_reset'}); 
        next_row(y);        
        
        PushButtonParam(obj, 'lickport_motor_set_ml_center', x, y, 'label', 'Set center');
        set_callback(lickport_motor_set_ml_center, {mfilename, 'lickport_motor_set_ml_center'}); 
         
        next_row(y);        
        next_row(y);        
        next_row(y);        
        PushButtonParam(obj, 'lickport_motor_set_ap_inreach', x, y, 'label', 'Set In Reach Position');
        set_callback(lickport_motor_set_ap_inreach, {mfilename, 'lickport_motor_set_ap_inreach'}); 
                
        
        next_column(x); y = 1;
        
        PushButtonParam(obj, 'read_pole_positions', x, y, 'label', 'Read pole position');
        set_callback(read_pole_positions, {mfilename, 'read_pole_positions'});
        
        next_row(y);
        NumeditParam(obj, 'ap_pole_motor_position', 0, x, y, 'label', ...
            'AP motor position','TooltipString','Absolute position in microsteps of motor.');
        set_callback(ap_pole_motor_position, {mfilename, 'ap_pole_motor_position'});
         
        next_row(y);
        NumeditParam(obj, 'ml_pole_motor_position', 0, x, y, 'label', ...
            'ML motor position','TooltipString','Absolute position in microsteps of motor.');
        set_callback(ml_pole_motor_position, {mfilename, 'ml_pole_motor_position'});
               
        next_row(y);
        SubheaderParam(obj, 'title', 'Read/set pole position', x, y);

        next_row(y);
        EditParam(obj, 'left_position', '180000', x, y, 'label', ...
            'Left position','TooltipString','Left port rewarded trial position in microsteps.');
        
        next_row(y);
       % NumeditParam(obj, 'right_position', 0, x, y, 'label', ...
       %     'Go position','TooltipString','Go trial position in microsteps.');
        EditParam(obj, 'right_position', '(3e4,9e4)', x, y, 'label', ...
            'Right position','TooltipString','Right port rewarded trial position in microsteps.');

        % bump on top
        next_row(y);
        NumeditParam(obj, 'inreach_moveby', 0, x, y, 'label', ...
            'In-reach move by','TooltipString','How much to move once pole is in reach in microsteps.');
        SoloParamHandle(obj, 'motor_inreach_moved', 'value', 0);
        
        next_row(y);
        SubheaderParam(obj, 'title', 'Trial AP pole position', x, y);
        
        %%% Lickport motor

        next_row(y, 2);        
        PushButtonParam(obj, 'read_lickport_positions', x, y, 'label', 'Read lickport position');
        set_callback(read_lickport_positions, {mfilename, 'read_lickport_positions'});
        
        next_row(y);
        NumeditParam(obj, 'ap_lickport_motor_position', 0, x, y, 'label', ...
            'AP motor position','TooltipString','Absolute position in microsteps of motor.');
        set_callback(ap_lickport_motor_position, {mfilename, 'ap_lickport_motor_position'});
         
        next_row(y);
        NumeditParam(obj, 'ml_lickport_motor_position', 0, x, y, 'label', ...
            'ML motor position','TooltipString','Absolute position in microsteps of motor.');
        set_callback(ml_lickport_motor_position, {mfilename, 'ml_lickport_motor_position'});
               
        next_row(y);
        SubheaderParam(obj, 'title', 'Read/set lickport position', x, y);

        % Positional raneg for auto positioning of motor for lickport
        next_row(y);
        NumeditParam(obj, 'ml_lickport_motor_cenpos', 90000, x, y, 'label', ...
            'ML position center','TooltipString','Center position for ML lickport motor in automode.');
        next_row(y);
        NumeditParam(obj, 'ml_lickport_motor_rangepos', 20000, x, y, 'label', ...
            'ML position range','TooltipString','Position range -- spans center +/- range/2.');
        next_row(y);

        % Modefor ML motor
        MenuParam(obj, 'MLMotorMode', {'Manual','Auto'}, ...
                    'Manual', x, y);
        
        
        % Settings for lickport withdrawal
        next_row(y);
        NumeditParam(obj, 'ap_lickport_motor_inreach_pos', 90000, x, y, 'label', ...
            'Inreach Pos','TooltipString','AP lick motor position animal can reach.');
        next_row(y);
        NumeditParam(obj, 'ap_lickport_wdraw_moveby', 10000, x, y, 'label', ...
            'Withdraw Distance','TooltipString','How many zaber units to move lickport by when withdrawn.');
        next_row(y);
        NumeditParam(obj, 'ap_lickport_wdraw_prob', 1, x, y, 'label', ...
            'Withdraw Probability','TooltipString','Fraction of trials to do withdrawal on');
        next_row(y);

        % Withdraw on?
        MenuParam(obj, 'LPWithdrawMode', {'off','postOnly','on','postAndPremature'}, ...
                    'off', x, y);
        next_row(y);
        SoloParamHandle(obj, 'motor_withdraw_moved', 'value', 0);
        
        %give make and upload smatix axxess to this
        SoloFunctionAddVars('make_and_upload_state_matrix', 'ro_args', 'LPWithdrawMode');   
        
        SubheaderParam(obj, 'title', 'Automatic Lickport Controls', x, y);
        
        % For debugging motor
        SoloParamhandle(obj, 'motor_move_time', 'value', 0);


        MotorsSection(obj,'hide_show');
        MotorsSection(obj,'read_pole_positions');
        MotorsSection(obj,'read_lickport_positions');
        
        x = parentfig_x; y = parentfig_y;
        set(0,'CurrentFigure',value(myfig));
        return;

    case 'move_next_side', % --------- CASE MOVE_NEXT_SIDE -----
        
        next_side = SidesSection(obj,'get_next_side');
        
        % --- the lickport -- positon 0 is left, 18e4 is right
        mlmm = value (MLMotorMode);
        lickportPos= value(ml_lickport_motor_position);

        if (strcmp(mlmm, 'Auto') & length(previous_sides) > 1)
            % get #trials used
            ntbc=value(NumTrialsBiasCalc);
            cenpos = value(ml_lickport_motor_cenpos);
            rangepos = value(ml_lickport_motor_rangepos);
            minpos = cenpos - rangepos/2; % LEFMOST
            maxpos = cenpos + rangepos/2; % RIGHTMOST
           
            rT = find(previous_sides(1:end-1) == 'r');
            lT = find(previous_sides(1:end-1) == 'l');            
              
            % determince frac correct for last # trials each side
            if (length(rT) >= ntbc & length(lT) >= ntbc)
                valL = find(hit_history(lT) >= 0);
                valR = find(hit_history(rT) >= 0);
                  
                % enuff VALID trials per side? i.e. with lick
                if (length(valL) >= ntbc & length(valR) >= ntbc)
                    % restrict to nTrials
                    valL = valL(end-ntbc+1:end);
                    valR = valR(end-ntbc+1:end);
                      
                    % compute frac correct
                    fL = sum(hit_history(lT(valL)))/length(valL);
                    fR = sum(hit_history(rT(valR)))/length(valR);
                      
                    % assign motor position
                    lpp = 0.5 - ((fL - fR)/2); % ranges from [0 1]
                    if (lpp < 0) ; lpp = 0; elseif (lpp > 1) ; lpp = 1 ; end
                    
                    % positino
                    finpos = minpos + (lpp)*(maxpos-minpos); % minpos is LEFT MOST
                    disp(['Using position (0=R 1=L): ' num2str(lpp) ' finpos ' num2str(finpos)]);

                    move_absolute(lickport_motors,finpos,lickport_motors_properties.ml_motornum);       
                    MotorsSection(obj,'read_lickport_positions');
                    lickportPos=finpos;
                end
            else % center
                move_absolute(lickport_motors,cenpos,lickport_motors_properties.ml_motornum);       
                MotorsSection(obj,'read_lickport_positions');
                lickportPos =cenpos;
            end
            
        end
         
     
        
        % --- the pole
        % pull position and process if range
        rightPoss = value(right_position);
        leftPos = value(left_position);
        
        % cam be a range specified by (a,b) where a is min b is max
        if (~isnumeric(rightPoss)) % assume format is correct!
          comIdx = find(rightPoss == ',');
          if (length(comIdx) > 0)
              minValueR = str2num(rightPoss(2:comIdx-1));
              maxValueR = str2num(rightPoss(comIdx+1:end-1));

              rightPoss = minValueR + round((maxValueR-minValueR)*rand(1));
          else
            disp('WARNING: right value format wrong');
            minValueR = rightPoss;
            maxValueR = leftPos;
          end
          
        else
          minValueR = rightPoss;
          maxValueR = leftPos;
        end
     
        if (~isnumeric(leftPos)) % assume format is correct!
          comIdx = find(leftPos == ',');
          minValueL = str2num(leftPos(2:comIdx-1));
          maxValueL = str2num(leftPos(comIdx+1:end-1));
          
          leftPos = minValueL + round((maxValueL-minValueL)*rand(1));
        else
          minValueL = leftPos;
          maxValueL = leftPos;
        end        
        % Manually start pedestal at mid-point (90000).
        if next_side == 'r'
              position = rightPoss;
        elseif next_side == 'l'
              position = leftPos;
        else
            error('Invalid next side.')
        end
        
        % OLD : halfway point is pre-position point
        halfpoint = abs(round((leftPos-rightPoss)/2)) + min(leftPos,rightPoss);
        prePoint = halfpoint;
        
        % NEW : preposition point is simply drawn RANDOMLY from net
        % position range
        minPos = min(minValueL,minValueR);
        maxPos = max(maxValueL,maxValueR);
        rangePos = maxPos-minPos;
        prePoint = minPos + round(rand(1)*rangePos);
       
        
%         if Solo_Try_Catch_Flag == 1
%             try
%                 move_absolute_sequence(motors,{90000,position},value(right_motor_num));
%             catch
%                 pause(1)
%                 warning('Error in move_absolute_sequence, forcing state 35...');
%                 SMControlSection(obj,'force_state_35');
%                 return
%             end
%         else
             tic
             move_absolute_sequence(pole_motors,{prePoint,position},pole_motors_properties.ap_motornum);
             movetime = toc;
             motor_move_time.value = movetime;
             if movetime<1.9 % Should make this min-ITI a SoloParamHandle
                 pause(1.9-movetime); %4
             end
%         end
        
         
        previous_sides(n_started_trials+1) = next_side;
        MotorsSection(obj,'read_pole_positions');
        
%         trial_ready_times.value = datestr(now);
        trial_ready_times.value = clock;
        
        pph = pole_position_history(:);
        pph(n_started_trials) = position;
        if (size(pph,1) == 1) ; pph = pph'; end
        
        % Store commanded position (don't rely on MotorsSection.position)
        pole_position_history.value = [pph];
     
        lph = lickport_position_history(:);
        lph(n_started_trials) = lickportPos;
        if (size(lph,1) == 1) ; lph = lph'; end
        
        % Store commanded position (don't rely on MotorsSection.position)
        lickport_position_history.value = [lph];
        return;
        
    case 'pole_motors_home',
        move_home(pole_motors);
        return;

    case 'pole_serial_open',
        serial_open(pole_motors);
        return;

    case 'pole_serial_reset',     
        close_and_cleanup(pole_motors);
        
        global pole_motors_properties;
        global pole_motors; 
        
        if strcmp(pole_motors_properties.type,'@FakeZaberTCD1000')
            pole_motors = FakeZaberTCD1000;
        else
            pole_motors = ZaberTCD1000;
        end
        
        serial_open(pole_motors);
        return;

    case 'pole_motors_stop',
        stop(pole_motors);
        return;

    case 'pole_motors_reset',
        reset(pole_motors);
        return;

    case 'reset_pole_motors_firmware',
        set_initial_parameters(pole_motors)
        display('Reset speed, acceleration, and motor bus ID numbers.')
        return;

    case 'ap_pole_motor_position',
        position = value(ap_pole_motor_position);
        if position > pole_motors_properties.ap_max_position | position < 0
            p = get_position(pole_motors,pole_motors_properties.ap_motornum);
            ap_pole_motor_position.value = p;
        else
            move_absolute(pole_motors,position,pole_motors_properties.ap_motornum);
        end
        return;


    case 'ml_pole_motor_position',
        position = value(ml_pole_motor_position);
        if position > pole_motors_properties.ml_max_position | position < 0
            p = get_position(pole_motors,pole_motors_properties.ml_motornum);
            ml_pole_motor_position.value = p;
        else
            move_absolute(pole_motors,position,pole_motors_properties.ml_motornum);
        end
        return;
    
    case 'read_pole_positions'
        p = get_position(pole_motors,pole_motors_properties.ap_motornum);
        ap_pole_motor_position.value = p;
        p = get_position(pole_motors,pole_motors_properties.ml_motornum);
        ml_pole_motor_position.value = p;
        return;
        
    case 'ap_lickport_motor_position',
        position = value(ap_lickport_motor_position);
        if position > lickport_motors_properties.ap_max_position | position < 0
            p = get_position(lickport_motors,lickport_motors_properties.ap_motornum);
            ap_lickport_motor_position.value = p;
        else
            move_absolute(lickport_motors,position,lickport_motors_properties.ap_motornum);
        end
        return;


    case 'ml_lickport_motor_position',
        position = value(ml_lickport_motor_position);
        if position > lickport_motors_properties.ml_max_position | position < 0
            p = get_position(lickport_motors,lickport_motors_properties.ml_motornum);
            ml_lickport_motor_position.value = p;
        else
            move_absolute(lickport_motors,position,lickport_motors_properties.ml_motornum);
        end
        return;
    
    case 'read_lickport_positions'
        p = get_position(lickport_motors,lickport_motors_properties.ap_motornum);
        ap_lickport_motor_position.value = p;
        p = get_position(lickport_motors,lickport_motors_properties.ml_motornum);
        ml_lickport_motor_position.value = p;
        return;

    
    case 'lickport_motor_set_ml_center'
        % read ML position
        p = get_position(lickport_motors,lickport_motors_properties.ml_motornum);
        
        % set value
        ml_lickport_motor_cenpos.value = p;
        return;

    
    case 'lickport_motor_set_ap_inreach'
        % read ML position
        p = get_position(lickport_motors,lickport_motors_properties.ap_motornum);
        
        % set value
        ap_lickport_motor_inreach_pos.value = p;
        return;
        
    case 'lickport_motors_home',
        move_home(lickport_motors);
        return;

    case 'lickport_serial_open',
        serial_open(lickport_motors);
        return;

    case 'lickport_serial_reset',     
        close_and_cleanup(lickport_motors);
        
        global lickport_motors_properties;
        global lickport_motors; 
        
        if strcmp(lickport_motors_properties.type,'@FakeZaberTCD1000')
            lickport_motors = FakeZaberTCD1000;
        else
            lickport_motors = ZaberTCD1000;
        end
        
        serial_open(lickport_motors);
        return;

    case 'lickport_motors_stop',
        stop(lickport_motors);
        return;

    case 'lickport_motors_reset',
        reset(lickport_motors);
        return;
      
    % --------- CASE HIDE_SHOW ---------------------------------

    case 'hide_show'
        if strcmpi(value(motor_show), 'hide')
            set(value(motorfig), 'Visible', 'off');
        elseif strcmpi(value(motor_show),'view')
            set(value(motorfig),'Visible','on');
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
        
    case 'update' 
        
        % check for sampling period -- if in it and user wants to move pole, do so.
        state = GetParam('rpbox','state','value');
        if (state == RealTimeStates.sampling_period)
            if (value(motor_inreach_moved) < 3)  
                motor_inreach_moved.value = value(motor_inreach_moved) + 1;
                if (value(motor_inreach_moved) == 3)
                    moveBy = value(inreach_moveby);
                    if (moveBy ~= 0)
                        pph = value(pole_position_history);
                        finalPos = pph(end) + moveBy;
                        if (finalPos > 180000) ; finalPos = 180000; end
                        if (finalPos < 0) ; finalPos = 0; end
                    
                        move_absolute(pole_motors,finalPos,pole_motors_properties.ap_motornum);
                    end               
                end
            end
        else
            motor_inreach_moved.value = 0;
        end
        
        % check for lickport withdrawal
        wdrawMode = value(LPWithdrawMode);
        wdrawProb = value(ap_lickport_wdraw_prob);
        if (strcmp(wdrawMode, 'on')) % moev in only during response period
             if (state == RealTimeStates.present_lickport)     % move IN
                if(value(motor_withdraw_moved) == 0 )
                    motor_withdraw_moved.value = 1;
                    moveBy = value(ap_lickport_wdraw_moveby);
                    if (moveBy ~= 0)
                        sPos = value(ap_lickport_motor_inreach_pos);
                        move_absolute(lickport_motors,sPos,lickport_motors_properties.ap_motornum);                    
                    end
                end
            elseif (state == RealTimeStates.posttrial_pause) % move OUT at end of trial
                 if(value(motor_withdraw_moved) == 1 )
                    motor_withdraw_moved.value = 0;
                    moveBy = value(ap_lickport_wdraw_moveby);
                    if (moveBy ~= 0)
                        if (rand <= wdrawProb)
                            sPos = value(ap_lickport_motor_inreach_pos);
                            move_absolute(lickport_motors,sPos + moveBy,lickport_motors_properties.ap_motornum);                    
                        end
                    end
                 end
             end
        elseif (strcmp(wdrawMode,'postOnly')) % lickport usually in, but @end of trial popout
             if (state == RealTimeStates.pretrial_pause)     % move IN
                if(value(motor_withdraw_moved) == 0 )
                    motor_withdraw_moved.value = 1;
                    moveBy = value(ap_lickport_wdraw_moveby);
                    if (moveBy ~= 0)
                        sPos = value(ap_lickport_motor_inreach_pos);
                        move_absolute(lickport_motors,sPos,lickport_motors_properties.ap_motornum);                    
                    end
                end
            elseif (state == RealTimeStates.posttrial_pause) % move OUT at end of trial
                 if(value(motor_withdraw_moved) == 1 )
                    motor_withdraw_moved.value = 0;
                    moveBy = value(ap_lickport_wdraw_moveby);
                    if (moveBy ~= 0)                        
                        if (rand <= wdrawProb)
                            sPos = value(ap_lickport_motor_inreach_pos);
                            move_absolute(lickport_motors,sPos + moveBy,lickport_motors_properties.ap_motornum);                    
                        end                  
                    end
                 end
             end
        elseif (strcmp(wdrawMode,'postAndPremature')) % withdraw at end / on premature lick
             if (state == RealTimeStates.pretrial_pause)     % move IN
                if(value(motor_withdraw_moved) == 0 )
                    motor_withdraw_moved.value = 1;
                    moveBy = value(ap_lickport_wdraw_moveby);
                    if (moveBy ~= 0)
                        sPos = value(ap_lickport_motor_inreach_pos);
                        move_absolute(lickport_motors,sPos,lickport_motors_properties.ap_motornum);                    
                    end
                end
            elseif (state == RealTimeStates.posttrial_pause | ...
                    state == RealTimeStates.withdraw_lickport_preans | ...
                    state == RealTimeStates.withdraw_lickport_sample ) % move OUT at end of trial, on premie lick
                 if(value(motor_withdraw_moved) == 1 )
                    motor_withdraw_moved.value = 0;
                    moveBy = value(ap_lickport_wdraw_moveby);
                    if (moveBy ~= 0)                        
                        if (rand <= wdrawProb)
                            sPos = value(ap_lickport_motor_inreach_pos);
                            move_absolute(lickport_motors,sPos + moveBy,lickport_motors_properties.ap_motornum);                    
                        end                  
                    end
                 end
             end
        end
        
        return;
end


