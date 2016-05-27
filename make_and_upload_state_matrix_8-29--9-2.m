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
   
   % Settings S. Peron
   wvLid = 2^6; % water valve (left)
   ledLid = 2^1; % lickport LED (left)
   wvRid = 2^7; % water valve (right)
   ledRid = 2^0; % lickport LED (right)
   pvid = 2^2; % Pneumatic (Festo) valve ID.
   etid = 2^4; % EPHUS (electrophysiology) trigger ID.
   slid = 2^5; % Signal line for signaling trial numbers and fiducial marks.
   rcid = 2^3; % Reward cue ID (sound)
   
   puffid = 0; % Airpuff valve ID. DISABLED
   
   rwvtm = RWaterValveTime; % Defined in ValvesSection.m.  
   lwvtm = LWaterValveTime; % Defined in ValvesSection.m.  
   
   % Compute answer period time as 2 sec minus SamplingPeriodTime (from TimesSection.m) , 
   % unless SamplingPeriodTime is > 1 s (for training purposes), in which case it is 1 sec.
   % MOVED THIS PART TO TIMESECTION. - NX 4/9/09
      
      
      
   % program starts in state 40
   stm = [0 0 0 0 40 0.01  0 0];
   stm = [stm ; zeros(40-rows(stm), 8)];
   stm(36,:) = [35 35 35 35 35 1   0 0];
   b = rows(stm); 

   RealTimeStates.bitcode = b; 
   RealTimeStates.pretrial_pause = b+1;
   RealTimeStates.sampling_period = b+2;
   RealTimeStates.preanswer_pause = b+3;
   RealTimeStates.answer_epoch = b+4;
   RealTimeStates.miss = b+5;
   RealTimeStates.posttrial_pause = b+6;
   RealTimeStates.punish = b+7;
   RealTimeStates.reward_left = b+8;
   RealTimeStates.reward_right = b+9;     
   RealTimeStates.second_chance = b+10; % if RewardOnWrong
   RealTimeStates.reward_cue = b+11;
   RealTimeStates.reward_catch_trial = b+12;
   RealTimeStates.reward_collection = b+13;
   RealTimeStates.withdraw_lickport_sample = b+14;
   RealTimeStates.withdraw_lickport_preans = b+15;
   RealTimeStates.present_lickport = b+16;
   RealTimeStates.restart_delay = b+17;   
   
   next_side = SidesSection(obj, 'get_next_side');
   
   % ----------------------------------------------------------------------
   % - Build matrix:
   % ----------------------------------------------------------------------
   switch SessionType % determined by SessionTypeSection
        
      case 'Water-Valve-Calibration'
  
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
           sAns = b+4; % answer epoch
           sLoMi =b+5 ; % log miss/ignore
           sPoTP = b+6; % posttrial pause
           sPun = b+7; % punish (be it airpuff, timeout, or whatev)  
           sRwL = b+8; % reward left
           sRwR = b+9; % reward right
           sSC = b+10 ; % second chance
           sRC = b+11; % reward cue
           sRCaT = b+12; % to log unrewarded correct trals (catch)
           sRCol = b+13; % give animal time to collect reward
           sWDSa = b+14 ; % pull the lickport
           sWDPr = b+15 ; % pull the lickport           
           sLPP = b+16 ; % insert the lickport           
           sRDel = b+17 ; % restart delay         
           
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
           
           %% LIckport travel time - subtract from preans
%            pa_t = pa_t-lpt_t;
%            if (pa_t < 0)
%                disp('WARNING: your PreAnswerTime is too short -- must include the lickport withdraw time.');
%                disp('Right now, LickportWithdrawTime > PreAnswerTime.  Cannot BE!');
%                pa_t = .01;
%            end
                 
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
       
           % Allow wrong? i.e., if animal licks wrong port FIRST, and this
           % is yes, then post timeout state actually allows rewards
           re_wr = value(RewardOnWrong);
           pps = sPoTP; % post-punish state default is post trial pause
           if (strcmp(re_wr, 'yes'))
               % timeout must be ZERO in this case.
               if (eto_t > .01)
                   disp('RewardOnWrong yes requires timeout of 0; setting to 0.');
                   eto_t = 0.01;
               end
               pps = sSC; % second chance now
           end
           
           % pull the lickport motor mode
%            lpwdm = value(LPWithdrawMode);
%            onLickS = sPMS;
%            onLickP = sPrAP; % on lick in the pre-answer state
%            if (strcmp(lpwdm, 'postAndPremature'))
%                onLickS = sWDSa;
%                onLickP = sWDPr;
%            end
           
           % Adjust extra time out basd on airpuf tmie
           eto_t = eto_t - puff_t;
           eto_t = max(eto_t,.01);
           
           %% trial-type dependent variables
           onlickSL = sSC; % by default, if you lick from 2nd chance state, you stay
           onlickSR = sSC; 
           if next_side=='r' % 'r' means right trial.
               onlickL = sPun; % incorrect
               onlickR = sRwR; % correct
               water_t = RWaterValveTime; % Defined in ValvesSection.m.  
               
               % if RewardOnWrong, from punish staet lick right gives
               % reward
               if (strcmp(re_wr, 'yes')) ; onlickSR = sRwR; end
           else %left 
               onlickR = sPun; % punish
               onlickL = sRwL; % water to left port
               water_t = LWaterValveTime; % Defined in ValvesSection.m.  
               
               % if RewardOnWrong, from punish staet lick left gives
               % reward
               if (strcmp(re_wr, 'yes')) ; onlickSL = sRwL; end
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
           
           % Restart PreAnswer Period due to lick?
%            if (strcmp(RestartPreAnsOnLick,'on'))
%                onLickP = sRDel;
%            end
            
           
           stm = [stm ;
               %LinSt   LoutSt   RinSt    RoutSt   TimeupSt Time      Dou      Aou  (Dou is bitmask format)
               % line b (sBC = b)
               sBC      sBC      sBC      sBC      101      .01       etid       0; ... % send bitcode
               sPrTP    sPrTP    sPrTP    sPrTP    sPMS     prep_t    0          0; ... % pretrial pause
              % onLickS  onLickS  onLickS  onLickS  sPrAP    pr_t+sp_t pvid       0; ... % Preanswer Pause
             %  onLickP  onLickP  onLickP  onLickP  sLPP     pr_t+pa_t 0          0; ... % Insert the lickport
               onlickL  onlickL  onlickR  onlickR  sLoMi    ap_t      0          0; ... % Check if correct lick
               sLoMi    sLoMi    sLoMi    sLoMi    sPoTP    0.001     0          0; ... % log miss/ignore
               sPoTP    sPoTP    sPoTP    sPoTP    35       postp_t   0          0; ... % posttrial pause
               sPun     sPun     sPun     sPun     pps      eto_t     0          0; ... % punish
               sRwL     sRwL     sRwL     sRwL     sRCol    water_t   wvLid      0; ... % reward left
               sRwR     sRwR     sRwR     sRwR     sRCol    water_t   wvRid      0; ... % reward right              
               onlickSL onlickSL onlickSR onlickSR pas      ap_t      0          0; ... % Lick from 2nd chance state, you stay              
               sRC      sRC      sRC      sRC      sAns     rc_t      rcid       0; ... % reward cue              
               sRCaT    sRCaT    sRCaT    sRCaT    sPoTP    0.001     0          0; ... % to log unrewarded correct trials              
               sRCol    sRCol    sRCol    sRCol    sPoTP    rcoll_t   0          0; ... % give animal time to collect reward              
            %   sWDSa    sWDSa    sWDSa    sWDSa    sPMS     wdraw_t   pvid       0; ... % pull the lickport              
             %  sWDPr    sWDPr    sWDPr    sWDPr    sPrAP    wdraw_t   0          0; ... % pull lickport              
           %    sLPP     sLPP     sLPP     sLPP     sRC      lpt_t     0          0; ... % insert the lickport 
               sRDel    sRDel    sRDel    sRDel    sPrAP    0.001     0          0; ... % restart delay              
               ];
           
           %------ Signal trial number on digital output given by 'slid':
           % Requires that states 101 through 101+2*numbits be reserved
           % for giving bit signal.
           
           trialnum = n_done_trials + 1;
           
           %             trialnum = 511; %63;
           % Should maybe make following 3 SPHs in State Machine Control
           % GUI:
           bittm = 0.002; % bit time
           gaptm = 0.005; % gap (inter-bit) time
           numbits = 10; %2^10=1024 possible trial nums
           
           
           x = double(dec2binvec(trialnum)');
           if length(x) < numbits
               x = [x; repmat(0, [numbits-length(x) 1])];
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

   