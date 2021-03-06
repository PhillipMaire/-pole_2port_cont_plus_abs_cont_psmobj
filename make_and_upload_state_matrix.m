function [] = make_and_upload_state_matrix(obj, action)

GetSoloFunctionArgs;


switch action
    case 'init'
        clear global autosave_session_id;
        SoloParamHandle(obj, 'state_matrix');

        SoloParamHandle(obj, 'RealTimeStates', 'value', struct(...
            'bitcode', 0, ...
            'pretrial_pause', 0, ...
            'sampling_period', 0, ...
            'preanswer_pause',  0, ...
            'answer_epoch', 0,...
            'miss',0, ...
            'posttrial_pause', 0, ...
            'punish', 0, ...
            'reward_left', 0, ...
            'reward_right',0, ...
            'second_chance',0, ...
            'reward_cue',0, ...
            'reward_catch_trial', 0, ...
            'reward_collection',0, ...
            'withdraw_lickport_sample',0, ...
            'withdraw_lickport_preans',0, ...
            'present_lickport', 0, ...
            'restart_delay', 0));

        SoloFunctionAddVars('RewardsSection', 'ro_args', 'RealTimeStates');
        SoloFunctionAddVars('MotorsSection', 'ro_args', 'RealTimeStates');

        make_and_upload_state_matrix(obj, 'next_matrix');

        return;

    case 'next_matrix',
        % SAVE EVERYTHING before anything else is done!
        %autosave(value(MouseName));
        SavingSection(obj, 'autosave');

        % ----------------------------------------------------------------------
        % - Set parameters used in matrix:
        % ----------------------------------------------------------------------
        % Setting Andrew
        %    jnk = 2^0
        %    jnk =2^1
        %    Festo =2^2
        %    jnk =2^3
        %    jnk =2^4
        %    jnk =2^5
        %    watervalve=2^6
        %    watervalve=2^7
        % R = mouse right
        % L = mouse left
        %
        % Settings S. Peron
        wvLid = 2^7; % water valve (left) % CHANGE PER BOX!!!!
%         ledLid = 2^1; % lickport LED (left)
        wvRid = 2^6; % water valve (right) % CHANGE PER BOX!!!!
%         ledRid = 2^1; % lickport LED (right)
        pvid = 2^2; % Pneumatic (Festo) valve ID.
        etid = 2^4; % EPHUS (electrophysiology) trigger ID.
        slid = 2^5; % Signal line for signaling trial numbers and fiducial marks.
%         rcid = 2^3; % Reward cue ID (sound)
        absentNoRewid= 2^3;%MAKE SURE THIS ISN'T BEING USED SO IT DOES NOTHING
%         puffid = 0; % Airpuff valve ID. DISABLED
        sSound1 = 2^0; %sound1 for absent trial licks (stop licking indicator)
        sSound2 = 2^1; %sound2 for answer period indicator

        rwvtm = RWaterValveTime; % Defined in ValvesSection.m.
        lwvtm = LWaterValveTime; % Defined in ValvesSection.m.

        % Compute answer period time as 2 sec minus SamplingPeriodTime (from TimesSection.m) ,
        % unless SamplingPeriodTime is > 1 s (for training purposes), in which case it is 1 sec.
        % MOVED THIS PART TO TIMESECTION. - NX 4/9/09



        % program starts in state 40
        stm = [0 0 0 0 40 0.01  0 0];%stm is first defined here
        stm = [stm ; zeros(40-rows(stm), 8)];%this just makes stm 40 rows and 8 columns
        stm(36,:) = [35 35 35 35 35 1   0 0];%obviously just edits column 36 (all rows)
        b = rows(stm);

        % Alexis 9-3-14 MUST BE RENAMED BASED ON 2-PORT DISCRIM TO MAKE SURE HITS
        % AND INCORRECTS ARE MONITORED (refer to assignment of state times in
        % parse_trial in RewardsSection in state35 in state RPbox via the send
        % matrix line in this m-file--> this is done in the pstruct variable
        % struct).
        %
        % Assigned unused states names with unused state numbers (DO NOT ASSIGN <--
        % fills in areas of the pstruct that we do not currently want filled
        RealTimeStates.bitcode = b;
        RealTimeStates.pretrial_pause = b+1;
        RealTimeStates.sampling_period = b+2;
        RealTimeStates.preanswer_pause = b+3;
        RealTimeStates.answer_epoch = b+12;%4;
        RealTimeStates.miss = b+4;%5;
        RealTimeStates.posttrial_pause = b+5;%6;
        RealTimeStates.punish = b+6;%7;
        RealTimeStates.reward_left = b+7;%8;
        RealTimeStates.reward_right = b+8;%9;
        RealTimeStates.second_chance = b+13;%10; % if RewardOnWrong
        RealTimeStates.reward_cue = b+14;%11;
        RealTimeStates.reward_catch_trial = b+9;%12;
        RealTimeStates.reward_collection = b+10;%13;
        RealTimeStates.withdraw_lickport_sample = b+15;%14;
        RealTimeStates.withdraw_lickport_preans = b+16;%15;
        RealTimeStates.present_lickport = b+17;%16;
        RealTimeStates.restart_delay = b+11;%17;
        %psm below
        %testingpsm   RealTimeStates.correct_rejection = b+18;
        %psm above
        next_side = SidesSection(obj, 'get_next_side');

        % ----------------------------------------------------------------------
        % - Build matrix:
        % ----------------------------------------------------------------------
        switch SessionType % determined by SessionTypeSection

            case 'Water-Valve-Calibration'

            case 'Manual-Training'
               
                sBC = b; % bitcode
                sPrTP = b+1; % pretrial pause
                sPMS = b+2; % pole move & sample period
                sPrAP = b+3; % preanswer pause
                sLoMi = b+4; % log miss/ignore
                sPoTP = b+5; % posttrial pause
                sPun = b+6; % punish (be it airpuff, timeout, or whatev)
                sRwL = b+7; % reward left
                sRwR = b+8; % reward right
                sRCaT = b+9;%12; % to log unrewarded correct trals (catch)
                sRCol = b+10;%13; % give animal time to collect reward
                sRDel = b+11;%17 ; % restart delay


                % ---- assign gui variables
                ap_t = value(AnswerPeriodTime);
                sp_t = SamplingPeriodTime;
                eto_t = max(.01,ExtraITIOnError);
                pr_t = PoleRetractTime;
                pa_t = PreAnswerTime;
                prep_t = PreTrialPauseTime;
                postp_t = PostTrialPauseTime;
                rc_t = RewardCueTime;
                rcoll_t = RewardCollectTime;
                lpt_t = LickportTravelTime;
                msp_t = ManualSampPeriod;

                puff_t = AirpuffTime;

                wdraw_t = 0.3; % how long to stay in withdraw state and allow its detection

                % Check for (min,max) PreAnswerTime
                if (~isnumeric(pa_t)) % assume format is correct!
                    comIdx = find(pa_t == ',');
                    if (length(comIdx) > 0)
                        minValuePAT = str2num(pa_t(2:comIdx-1));
                        maxValuePAT = str2num(pa_t(comIdx+1:end-1));
                        pa_t = minValuePAT + ((maxValuePAT-minValuePAT)*rand(1));
                    else
                        disp('Bad format ; setting PreAnswerTime to 0');
                        pa_t= 0;
                    end
                end

                % puff disabled? can't have 0 time or RT freaks; set valve id
                % off instead
                if (puff_t == 0)
                    puff_t = 0.01;
                    puffid = 0;
                else
                    disp('*** At this time, airpuff is disabled. ***');
                end

                % Reward cue disabled? turn off valve, but must be at least
                % minimal time long
                if (rc_t == 0)
                    rc_t = 0.01;
                    rcid = 0;
                end

                % Adjust prepause based on bitcode, initial trigger
                prep_t = prep_t - 0.01 - 0.07; % 2 ms bit, 5 ms interbit, 10 bits = 70 ms = .07 s

                % Alexis 9-2-14 DECLARE A WRONG
                pps = sPoTP; % post-punish state default is post trial pause

                onLickS = sPMS;%helps define if there is lick during the sample period.

                % Adjust extra time out basd on airpuf tmie
%                 eto_t = eto_t - puff_t; %%%%tooke this out because dont
%                 use puff anymore, used 1 ms 0.001 time of the trigger for
%                 sound cue/restart time instead
                eto_t = max(eto_t,.01);

                if next_side=='r'; % 'r' means right trial.
                    onlickL = sPun; % incorrect
                    onlickR = sRwR; % correct
                    water_t = RWaterValveTime; % Defined in ValvesSection.m.
                elseif next_side=='l'; %left
                    onlickR = sPun; % punish
                    onlickL = sRwL; % water to left port
                    water_t = LWaterValveTime; % Defined in ValvesSection.m.
                else next_side=='a'; %for absent condition -psm
                    onlickL = sPun;
                    onlickR = sPun;
                    water_t = 0.01; %set so that forced state when absent wont give infinity water (if 0) or reward on absent
                end

                % Disable reward?  If so, do it by setting wv
                pas = sPoTP;
                if (FracNoReward > 0)
                    rvReward = rand;
                    if (rvReward < FracNoReward)
                        disp('Random cancellation of reward.')
                        wvLid = 0;
                        wvRid = 0;
                        pas = sRCaT;
                    end
                end
                rewVid=0;
                if onlickR==48
                    rewVid=wvRid;
                elseif onlickL==47
                    rewVid=wvLid;
                elseif onlickR==46 && onlickL==46
                    rewVid=absentNoRewid;
                else
                    error('L or R port id not 47 or 48 check make and upload state matrix')
                end

                 % Restart PreAnswer Period due to lick?
                if (strcmp(RestartPreAnsOnLick,'on'))
                    onLickP = sRDel;
                elseif (strcmp(RestartPreAnsOnLick,'off'))
                    onLickP = sPun;
                end

                stm = [stm ;
                    %LinSt   LoutSt   RinSt    RoutSt   TimeupSt Time      Dou      Aou  (Dou is bitmask format)
                    % line b (sBC = b)
                    sBC      sBC      sBC      sBC      101      .01        etid              0; ... %40 send bitcode
                    sPrTP    sPrTP    sPrTP    sPrTP    sPMS     prep_t     0                 0; ... %41 pretrial pause %Possibly sPMS -> sAns
                    onLickS  onLickS  onLickS  onLickS  sPrAP    pr_t+msp_t pvid              0; ... %42 Preanswer Pause
                    onlickL  onlickL  onlickR  onlickR  sLoMi    ap_t       pvid              0; ... %43 Check if correct lick
                    sLoMi    sLoMi    sLoMi    sLoMi    sPoTP    0.001      0                 0; ... %44 log miss/ignore
                    sPoTP    sPoTP    sPoTP    sPoTP    35       postp_t    0                 0; ... %45 posttrial pause
                    onLickP  onLickP  onLickP  onLickP  pps      eto_t      pvid              0; ... %46 punish
                    sRwL     sRwL     sRwL     sRwL     sRCol    water_t    pvid+wvLid        0; ... %47 reward left
                    sRwR     sRwR     sRwR     sRwR     sRCol    water_t    pvid+wvRid        0; ... %48 reward right
                    sRCaT    sRCaT    sRCaT    sRCaT    sPoTP    0.001      pvid              0; ... %49 to log unrewarded correct trials
                    sRCol    sRCol    sRCol    sRCol    sPoTP    rcoll_t    pvid              0; ... %50 give animal time to collect reward
                    sRDel    sRDel    sRDel    sRDel    sPun     0.001      pvid              0; ... %51 restart delay
                    52       52       52       52       sRCol    water_t    pvid+rewVid       0; ... %52 reward correct port
                    ];

                trialnum = n_done_trials + 1;

                bittm = 0.002; % bit time
                gaptm = 0.005; % gap (inter-bit) time
                numbits = 10; %2^10=1024 possible trial nums


                x = double(dec2binvec(trialnum)');%Convert decimal value to binary vector
                if length(x) < numbits
                    x = [x; repmat(0, [numbits-length(x) 1])];
                    %builds 0's on x so that it is
                    %the proper size for the bit code
                end
                % x is now 10-bit vector giving trial num, LSB first (at top).
                x(x==1) = slid;

                % Insert a gap state between bits, to make reading bit pattern clearer:
                x=[x zeros(size(x))]';
                x=reshape(x,numel(x),1);

                y = (101:(100+2*numbits))';
                t = repmat([bittm; gaptm],[numbits 1]);
                m = [y y y y y+1 t x zeros(size(y))];
                m(end,5) = sPrTP; % jump back to PREPAUSE.

                stm = [stm; zeros(101-rows(stm),8)];
                stm = [stm; m];



            case 'Licking'
                %Lin     Lout  Rin Rout   Tup  Tim  Dou Aou  (Dou is bitmask format)
                stm = [stm ;
                    b+1   b    b+2  b     35  999    0     0  ; ... % wait for lick  (This is state 40)
                    b+1   b+1  b+1  b+1   35  lwvtm  wvLid     0  ; ... % licked left -- reward left
                    b+2   b+2  b+2  b+2   35  rwvtm  wvRid   0  ; ... % licked right -- reward right
                    ];

            case 'Beam-Break-Indicator'
                stm = [stm ;
                    b+1   b  b+2 b    35  999  0      0  ; ...
                    b+1   b  b+1 b    35  999  ledLid  0  ; ...
                    b+2   b  b+2 b    35  999  ledRid  0  ; ...
                    ];


                % 2 port task
            case '2port-Discrim'

                % ---- labeling of states
                sBC = b; % bitcode
                sPrTP = b+1; % pretrial pause
                sPMS = b+2; % pole move & sample period
                sPrAP = b+3; % preanswer pause
                sLoMi = b+4; % log miss/ignore
                sPoTP = b+5; % posttrial pause
                sPun = b+6; % punish (be it airpuff, timeout, or whatev)
                sRwL = b+7; % reward left
                sRwR = b+8; % reward right
                sRCaT = b+9;%12; % to log unrewarded correct trals (catch)
                sRCol = b+10;%13; % give animal time to collect reward
                sRDel = b+11;%17 ; % restart delay


                % ---- assign gui variables
                ap_t = value(AnswerPeriodTime);
                sp_t = SamplingPeriodTime;
                eto_t = max(.01,ExtraITIOnError);
                pr_t = PoleRetractTime;
                pa_t = PreAnswerTime;
                prep_t = PreTrialPauseTime;
                postp_t = PostTrialPauseTime;
                rc_t = RewardCueTime;
                rcoll_t = RewardCollectTime;
                lpt_t = LickportTravelTime;

                puff_t = AirpuffTime;

                wdraw_t = 0.3; % how long to stay in withdraw state and allow its detection

                % Check for (min,max) PreAnswerTime
                if (~isnumeric(pa_t)) % assume format is correct!
                    comIdx = find(pa_t == ',');
                    if (length(comIdx) > 0)
                        minValuePAT = str2num(pa_t(2:comIdx-1));
                        maxValuePAT = str2num(pa_t(comIdx+1:end-1));
                        pa_t = minValuePAT + ((maxValuePAT-minValuePAT)*rand(1));
                    else
                        disp('Bad format ; settig PreAnswerTime to 0');
                        pa_t= 0;
                    end
                end

                % puff disabled? can't have 0 time or RT freaks; set valve id
                % off instead
                if (puff_t == 0)
                    puff_t = 0.01;
                    puffid = 0;
                else
                    disp('*** At this time, airpuff is disabled. ***');
                end

                % Reward cue disabled? turn off valve, but must be at least
                % minimal time long
                if (rc_t == 0)
                    rc_t = 0.01;
                    rcid = 0;
                end

                % Adjust prepause based on bitcode, initial trigger
                prep_t = prep_t - 0.01 - 0.07; % 2 ms bit, 5 ms interbit, 10 bits = 70 ms = .07 s

                % Alexis 9-2-14 DECLARE A WRONG
                pps = sPoTP; % post-punish state default is post trial pause

                onLickS = sPMS;%helps define if there is lick during the sample period.

                % Adjust extra time out basd on airpuf tmie
                eto_t = eto_t - puff_t;
                eto_t = max(eto_t,.01);
             
                 restartOnAllLicks = 1;
                if strcmp(BeepOn, 'all_wrong')
                    LRset  = 53;% 53 is the sound cue
                    ABSset = 53;
                elseif strcmp(BeepOn, 'abs_wrong')
                    LRset  = 46;% 46 is punish -psm
                    ABSset = 53;
                elseif strcmp(BeepOn, 'l-r_wrong')
                    LRset  = 53;
                    ABSset = 46;
                elseif strcmp(BeepOn, 'off')
                    LRset  = 46;
                    ABSset = 46;
                elseif strcmp(BeepOn, 'wrong_only_l-r')
                    LRset  = 53;
                    ABSset = 46;
                    restartOnAllLicks = 0;
                elseif strcmp(BeepOn, 'wrong_only_all')
                    LRset  = 53;% 53 is the sound cue
                    ABSset = 53;
                    restartOnAllLicks = 0;
                end

                if next_side=='r'; % 'r' means right trial.
                    onlickL = LRset; % incorrect
                    onLickLP = 53; 
                    onLickRP = 46;%stay on punish dont restart
                    onlickR = sRwR; % correct
                    water_t = RWaterValveTime; % Defined in ValvesSection.m.
                    restartState = LRset;
                elseif next_side=='l'; %left
                    onlickR  = LRset; % punish
                    onLickLP = 46; %stay on punish dont restart
                    onLickRP = 53; %make beep do restart. these depend on the 
                    %correct BeepOn settings being set
                    onlickL  = sRwL; % water to left port
                    water_t  = LWaterValveTime; % Defined in ValvesSection.m.
                    restartState = LRset; 
                else next_side=='a'; %for absent condition -psm
                    onlickL = ABSset; % 53 for sound cue dependent on absent trial
                    onLickLP = 53; 
                    onLickRP = 53;
                    onlickR = ABSset;
                    water_t = 0.01;
                    restartState = ABSset; 
                end

                % Disable reward?  If so, do it by setting wv
                pas = sPoTP;
                if (FracNoReward > 0)
                    rvReward = rand;
                    if (rvReward < FracNoReward)
                        disp('Random cancellation of reward.')
                        wvLid = 0;
                        wvRid = 0;
                        pas = sRCaT;
                    end
                end
                
                rewVid=0;
                if onlickR==48
                    rewVid=wvRid;
                elseif onlickL==47
                    rewVid=wvLid;
                elseif onlickR==46 && onlickL==46
                    rewVid=absentNoRewid;
                elseif onlickR==53 && onlickL==53
                    rewVid=absentNoRewid;
                else
                    error('L or R port id not 46 or 47 or 48 or 53 check make and upload state matrix')
                end

                % Restart PreAnswer Period due to lick?
                    %if onLickLP/onLickRP set earlier dont set now 
                    %restartOnAllLicks determines this. -psm
                if (strcmp(RestartPreAnsOnLick,'on')) && restartOnAllLicks
                    onLickLP = sRDel;
                    onLickRP = sRDel;
                elseif (strcmp(RestartPreAnsOnLick,'off'))
                    onLickLP = sPun;
                    onLickRP = sPun;
                end

               stm = [stm ;
                    %LinSt   LoutSt   RinSt    RoutSt   TimeupSt       Time      Dou      Aou  (Dou is bitmask format)
                    % line b (sBC = b)
                    sBC       sBC       sBC       sBC       101            .01         etid              0; ... %40 send bitcode
                    sPrTP     sPrTP     sPrTP     sPrTP     sPMS           prep_t      0                 0; ... %41 pretrial pause %Possibly sPMS -> sAns
                    onLickS   onLickS   onLickS   onLickS   54             pr_t+sp_t   pvid              0; ... %42 Preanswer Pause
                    onlickL   onlickL   onlickR   onlickR   sLoMi          ap_t        pvid              0; ... %43 Check if correct lick
                    sLoMi     sLoMi     sLoMi     sLoMi     sPoTP          0.001       0                 0; ... %44 log miss/ignore
                    sPoTP     sPoTP     sPoTP     sPoTP     35             postp_t     0                 0; ... %45 posttrial pause
                    onLickLP  onLickLP  onLickRP  onLickRP  pps            eto_t       pvid              0; ... %46 punish
                    sRwL      sRwL      sRwL      sRwL      sRCol          water_t     pvid+wvLid        0; ... %47 reward left
                    sRwR      sRwR      sRwR      sRwR      sRCol          water_t     pvid+wvRid        0; ... %48 reward right
                    sRCaT     sRCaT     sRCaT     sRCaT     sPoTP          0.001       pvid              0; ... %49 to log unrewarded correct trials
                    sRCol     sRCol     sRCol     sRCol     sPoTP          rcoll_t     pvid              0; ... %50 give animal time to collect reward
                    sRDel     sRDel     sRDel     sRDel     restartState   0.001       pvid              0; ... %51 restart delay
                    52        52        52        52        sRCol          water_t     pvid+rewVid       0; ... %52 reward correct port
                    53        53        53        53        sPun           0.005       pvid+sSound1      0; ... %53 lick on incorrect sound
                    54        54        54        54        sPrAP          0.005       pvid+sSound2      0; ... %54 answer period sound cue
                    ];


                trialnum = n_done_trials + 1;

                bittm = 0.002; % bit time
                gaptm = 0.005; % gap (inter-bit) time
                numbits = 10; %2^10=1024 possible trial nums


                x = double(dec2binvec(trialnum)');%Convert decimal value to binary vector
                if length(x) < numbits
                    x = [x; repmat(0, [numbits-length(x) 1])];
                    %builds 0's on x so that it is
                    %the proper size for the bit code
                end
                % x is now 10-bit vector giving trial num, LSB first (at top).
                x(x==1) = slid;

                % Insert a gap state between bits, to make reading bit pattern clearer:
                x=[x zeros(size(x))]';
                x=reshape(x,numel(x),1);

                y = (101:(100+2*numbits))';
                t = repmat([bittm; gaptm],[numbits 1]);
                m = [y y y y y+1 t x zeros(size(y))];
                m(end,5) = sPrTP; % jump back to PREPAUSE.

                stm = [stm; zeros(101-rows(stm),8)];
                stm = [stm; m];



            otherwise
                error('Invalid training session type')
        end

        stm = [stm; zeros(512-rows(stm),8)];


        rpbox('send_matrix', stm);
        state_matrix.value = stm;
        return;


    case 'reinit',
        % Delete all SoloParamHandles who belong to this object and whose
        % fullname starts with the name of this mfile:
        delete_sphandle('owner', ['^@' class(obj) '$'], ...
            'fullname', ['^' mfilename]);

        % Reinitialise
        feval(mfilename, obj, 'init');


    otherwise
        error('Invalid action!!');

end;

