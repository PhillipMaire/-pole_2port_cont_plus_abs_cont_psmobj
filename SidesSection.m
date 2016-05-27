% [x, y] = SidesSection(obj, action, x, y)
% Section that takes care of choosing the next correct side and keeping
% track of a plot of sides and hit/miss history.
%
% PARAMETERS:
% -----------
%
% obj      Default object argument.
%
% action   One of:
%            'init'      To initialise the section and set up the GUI
%                        for it; also calls 'choose_next_side' and
%                        'update_plot' (see below)
%
%            'reinit'    Delete all of this section's GUIs and data,
%                        and reinit, at the same position on the same
%                        figure as the original section GUI was placed.
%
%            'choose_next_side'  Picks what will be the next correct
%                        side. 
%
%            'get_next_side'  Returns either 'l' for left or 'r' for right.
%
%            'update_plot'    Update plot that reports on sides and hit
%                        history
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

function [x, y] = SidesSection(obj, action, x, y)
   
   GetSoloFunctionArgs;
   
   switch action
    
    case 'init',   % ------------ CASE INIT ----------------
      % Save the figure and the position in the figure where we are
      % going to start adding GUI elements:
      SoloParamHandle(obj, 'my_gui_info', 'value', [x y gcf]);

      % List of intended correct sides
      SoloParamHandle(obj, 'previous_sides', 'value', []);
      
      % Give read-only access to AnalysisSection.m:
      SoloFunctionAddVars('AnalysisSection', 'ro_args', 'previous_sides');
      SoloFunctionAddVars('MotorsSection', 'ro_args', 'previous_sides');      
      
      % 'Auto trainer' max # of FAs
      NumeditParam(obj, 'AutoTrainMinCorrect', 3, x, y); 
      next_row(y);
      
      NumeditParam(obj, 'AutoTrainMaxErrors', 3, x, y); 
      next_row(y);
      
      NumeditParam(obj, 'NumTrialsBiasCalc', 10, x, y,'TooltipString',...
          'How many trials to use -- per side -- in probabalistic biasing for motor or trial type seelction'); 
      SoloFunctionAddVars('MotorsSection', 'ro_args', 'NumTrialsBiasCalc');
      
      next_row(y);
      
      % If you are in brutal mode, this side is being used ...
      SoloParamHandle(obj, 'brutal_side', 'value', []);
      
      % Give read-only access to AnalysisSection.m:
      SoloFunctionAddVars('AnalysisSection', 'ro_args', 'brutal_side');
      
      % Autotrainer mode
      MenuParam(obj, 'AutoTrainMode', {'Off', 'Probabalistic','Alternate','Brutal'}, ...
                'Off', x, y);
      next_row(y);
                
      % Allow animal to correct reward from CORRECT port even if incorrect
      % port was licked?
      MenuParam(obj, 'RewardOnWrong', {'yes','no'}, ...
                'no', x, y, 'TooltipString', ...
                'If yes, animal can lick wrong port but still get reward if it licks correct port before time is up.');
      next_row(y);  
      % Give read-only access to make_and_upload_state_matrix.m:
      SoloFunctionAddVars('make_and_upload_state_matrix', 'ro_args', 'RewardOnWrong');
      
      % Max number of times same side can appear
      MenuParam(obj, 'MaxSame', {'1' '2' '3' '4' '5' '6' '7' '8' 'Inf'}, ...
                '3', x, y);
      next_row(y);
      % Prob of choosing left as correct side
      NumeditParam(obj, 'LeftPortProb', 0.3, x, y); 
      next_row(y, 1);
      NumeditParam(obj, 'RightPortProb', 0.3, x, y); 
      next_row(y, 1);
%have to set so that 
%if RightPortProb(1)+LeftPortProb(1)+OutReachProb(1)~=1
    %based on most recent updated gui of these three -psm
      
      
      SubheaderParam(obj, 'sidestitle', 'Trial Type Control', x, y);
      next_row(y);

      pos = get(gcf, 'Position');
      SoloParamHandle(obj, 'myaxes', 'saveable', 0, 'value', axes);
      set(value(myaxes), 'Units', 'pixels');
      set(value(myaxes), 'Position', [90 pos(4)-140 pos(3)-130 100]);
      set(value(myaxes), 'YTick', [1 2], 'YLim', [0.5 2.5], 'YTickLabel', ...
                        {'Right', 'Left'});
      NumeditParam(obj, 'ntrials', 100, x, y, ...
                   'position', [5 pos(4)-100 40 40], 'labelpos', 'top', ...
                   'TooltipString', 'How many trials to show in plot');
      set_callback(ntrials, {mfilename, 'update_plot'});      
      xlabel('trial number');
      SoloParamHandle(obj, 'previous_plot', 'saveable', 0);
      
      SidesSection(obj, 'choose_next_side');
      SidesSection(obj, 'update_plot');
      
      
      
    case 'choose_next_side', % --------- CASE CHOOSE_NEXT_SIDE -----
%called on by the primary program (with the same name as the @folder in the
%protocols -psm 
      % 108/l : nogo ; 114/r: go         111/o 
      %these are the just l and r's respective ASCII valies
      %ie char(114)==double('r'); char(108)==double('l'); char(111)==double('o');-psm 
      % --- autotrain mode is key to deciding what to do ...
      pickAtRandom = 0; % if 1, will simply use leftPortProb
      lpp = value(LeftPortProb); % this is changed by probabalistic autotrainer
      rpp = value(RightPortProb);
      oorp=1-rpp-lpp;
      ntbc = value (NumTrialsBiasCalc);
      
      switch lower(value(AutoTrainMode))
          % -- NO AUTOTRAINER -- just use maxSame and leftPortProb
          case 'off'
              pickAtRandom = 1;
              
          % -- Alternate: simple autotrainer where, after
          % AutoTrainMinCorrect licks are made, the autotrainer switches to
          %  the other side ; default is right
          case 'alternate'
              atmc = value(AutoTrainMinCorrect);
              next_side = 'r'; 

              % how many CORRECT licks in a row are we at?
              licks = find(hit_history == 1);
              if (length(licks) >= atmc)
                  last_side = previous_sides(licks(end));
                  cmp_mat = repmat(last_side, 1,atmc);
                  
                  % have we repeated last-side enough times?
                  if (sum(previous_sides(licks((end-atmc+1):end)) == cmp_mat) == atmc)
                      if (last_side == 'r') ; next_side = 'l'; end
                  else
                      next_side = last_side;
                  end
              end
              
          % -- Brutal: if animal makes MaxErrors in a row on either side,
          % present only that side until animal makes minCorect correct
          % responses
          case 'brutal'
              atmc = value(AutoTrainMinCorrect);
              atme = value(AutoTrainMaxErrors);
              next_side = '';
              
              % do we have enough CONSECUTIVE errors on eitherside
              % (consider only licks, not misses)
              cmp_mat = repmat(2,1,atme);
              rT = find(previous_sides == 'r');
              lT = find(previous_sides == 'l');
              
              % Are we already forcing one side? in this case, see if the
              %  minCorrect condition is met
              if (length(brutal_side) > 0)
                % if minCorrect condition is *not* met, next_side =
                % brutal_side
                next_side = value(brutal_side);

                if (brutal_side == 'r')
                    hh = hit_history(rT);
                else
                    hh = hit_history(lT);
                end
                
                % consider only licks
                val = find(hh >=0);
                
                % enuff trials?
                if (length(val) >= atmc)
                    % are ALL the last minCorrect responses correct?
                    if (sum(hh(val((end-atmc+1):end))) == atmc) 
                        next_side = ''; % this will disable brutality and resume pick @ random
                        disp('Disabling brutal mode');
                    else
                        disp(['Continuing brutal mode: ' num2str(next_side)]);
                    end
                end    
              
              % - NO brutal side
              else
                  % first see if we force right
                  if (length(rT) >= atme)
                      hhR = hit_history(rT);
                      val = find(hhR >= 0); % only consider lick trials
                      if (length(val) >= atme) % enough licks?
                          if (sum(hhR(val((end-atme+1):end))) == 0) % most recent all err?
                              next_side = 'r'; % then we need to repeat right.
                              brutal_side.value= 'r';
                              disp('Brutal mode R');                              
                          end
                      end
                  end

                  % no decision? then see if force left
                  if (length(next_side) == 0)
                       if (length(lT) >= atme)
                          if (length(hit_history) < max(lT))
                              hit_history(max(lT)) = 0;
                              disp('SidesSection::had to add a FALSE hit_history ; something was wrong.');
                          end
                          hhL = hit_history(lT);
                          val = find(hhL >= 0); % only consider lick trials
                          if (length(val) >= atme) % enough licks?
                              if (sum(hhL(val((end-atme+1):end))) == 0) % most recent all err?
                                  next_side = 'l'; % then we need to repeat left.
                                  brutal_side.value = 'l';
                                  disp('Brutal mode L');
                              end
                          end
                      end
                  end
              end
              
              % no forcing? then random
              if (length(next_side) == 0) ; pickAtRandom = 1; end
              
          % -- Probabalistic: pL = pL - (fL-fR)/2, where pL is probabiltiy
          % of left port, fL is fraction of left port correct and fR is
          % fraction of right port trials correct.  Force 0 < pL < 1 and
          % also note that must have 10 trials PER SIDE for this to work.
          case 'probabalistic'
              pickAtRandom = 1; % use it, but tweak it
              nTrials = ntbc; % how many per side to use
              
              rT = find(previous_sides == 'r');
              lT = find(previous_sides == 'l');            
              
              % determince frac correct for last nTrials each side
              if (length(rT) >= nTrials & length(lT) >= nTrials)
                  valL = find(hit_history(lT) >= 0);
                  valR = find(hit_history(rT) >= 0);
                  
                  % enuff VALID trials per side? i.e. with lick
                  if (length(valL) >= nTrials & length(valR) >= nTrials)
                      % restrict to nTrials
                      valL = valL(end-nTrials+1:end);
                      valR = valR(end-nTrials+1:end);
                      
                      % compute frac correct
                      fL = sum(hit_history(lT(valL)))/length(valL);
                      fR = sum(hit_history(rT(valR)))/length(valR);
                      
                      % bias
                      lpp = lpp - ((fL - fR)/2);
                      disp(['Using left probabiliy: ' num2str(lpp)]);
                      if (lpp < 0) ; lpp = 0; elseif (lpp > 1) ; lpp = 1 ; end
                  end
              end
      end
      
      %%%%%%%%%%%%%%%end of autotrainer switch%%%%%%%%%%%%%%%%

       % --- if you are to simply pick at random ...
       %pickAtRandom is set to 1 when  autotrainer is off (or using probabalistic autotrainer
       % which updates the lpp based on probabilitty), thus the 
       %below if statment just tests for autotrainer dropdown box selection
       if(pickAtRandom)

          brutal_side.value = ''; % in case there is something there
               
          % If MaxSame doesn't apply yet, choose at random
          if strcmp(value(MaxSame), 'Inf') | MaxSame > n_started_trials,
             
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-psm 
              randVar=rand(1); %needs to be a variable so that it can be used in two
              %consecutive if statemnts-psm
              if randVar<=lpp, next_side = 'l'; 
               elseif randVar>lpp && randVar<=lpp+rpp, next_side = 'r'; %psm
               else next_side = 'o' %psm%#change 'r' back to 'o'
             %psm else next_side = 'r'; 
             end;
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          else 
             % MaxSame applies, check for its rules:
             % If there's been a string of MaxSame guys all the same, force change:
             if all(previous_sides(n_started_trials-MaxSame+1:n_started_trials) == ...
                    previous_sides(n_started_trials))
                if previous_sides(n_started_trials)=='l' 
                    %%%%%-psm takes care of max same by choosing one of the other 
                    %two variables based on their relitive probabilities 
                     if rand(1)<=rpp/(rpp+oorp), next_side = 'r'; 
                     else next_side = 'o';%#change to 'o'
                     end
                elseif previous_sides(n_started_trials)=='r' 
                     if rand(1)<=lpp/(lpp+oorp), next_side = 'l'; 
                     else next_side = 'o';%#change to 'o'
                     end
                elseif previous_sides(n_started_trials)=='o' %#change to 'o'
                     if rand(1)<=lpp/(lpp+rpp), next_side = 'l'; 
                     else next_side = 'r';
                     end
                end;
             else
                % Haven't reached MaxSame limits yet, choose at random:
                randVar=rand(1);
                if randVar<=lpp, next_side = 'l'; 
                  elseif randVar>lpp && randVar<=lpp+rpp, next_side = 'r'; %psm
                  else next_side = 'o' %#change 'r' back to 'o'
                end;
             end;
          end;
       end

      
    %  session_type = SessionTypeSection(obj,'get_session_type'); 
      switch SessionType
          
          case {'Licking','Pole-conditioning'}
            next_side = 'r'; % Make it always the go-trial position, so mouse doesn't have to unlearn anything.
      end

      % logging globals (also used sometmies)
      previous_sides(n_started_trials+1) = next_side;
     
      lpph = left_port_prob_history(:);
      lpph(n_started_trials+1) = lpp;
      if (size(lpph,1) == 1) ; lpph = lpph'; end
      left_port_prob_history.value = [lpph];
      
    case 'get_next_side',   % --------- CASE GET_NEXT_SIDE ------
      if isempty(previous_sides),
         error('Don''t have next side chosen! Did you run choose_next_side?');
      end;
      x = previous_sides(length(previous_sides));
      return;
      
      
    case 'update_plot',     % --------- UPDATE_PLOT ------
      if ~isempty(value(previous_plot)), delete(previous_plot(:)); end;
      if isempty(previous_sides), return; end;

      % BLUE for upcoming/current
      ps = value(previous_sides);
      if ps(end)=='l', 
         hb = line(length(previous_sides), 2, 'Parent', value(myaxes));
      else                         
         hb = line(length(previous_sides), 1, 'Parent', value(myaxes));
      end;
      set(hb, 'Color', 'b', 'Marker', '.', 'LineStyle', 'none');
      
      % GEREN markers for correct
      xgreen = find(hit_history == 1);
      lefts  = find(previous_sides(xgreen) == 'l');
      rghts  = find(previous_sides(xgreen) == 'r');
      ygreen = zeros(size(xgreen)); ygreen(lefts) = 2; ygreen(rghts) = 1;
      hg = line(xgreen, ygreen, 'Parent', value(myaxes));
      set(hg, 'Color', 'g', 'Marker', '.', 'LineStyle', 'none'); 
 
      % RED markers for incorrect
      xred  = find(hit_history == 0);
      lefts = find(previous_sides(xred) == 'l');
      rghts = find(previous_sides(xred) == 'r');
      yred = zeros(size(xred)); yred(lefts) = 2; yred(rghts) = 1;
      hr = line(xred, yred, 'Parent', value(myaxes));
      set(hr, 'Color', 'r', 'Marker', '.', 'LineStyle', 'none'); 
      
      % BLACK x for miss/ignore
      xblack  = find(hit_history == -1);
      lefts = find(previous_sides(xblack) == 'l');
      rghts = find(previous_sides(xblack) == 'r');
      yblack = zeros(size(xblack)); yblack(lefts) = 2; yblack(rghts) = 1;
      hk = line(xblack, yblack, 'Parent', value(myaxes));
      set(hk, 'Color', 'k', 'Marker', 'x', 'LineStyle', 'none'); 

      previous_plot.value = [hb ; hr; hg ; hk];

      minx = n_done_trials - ntrials; if minx < 0, minx = 0; end;
      maxx = n_done_trials + 2; if maxx <= ntrials, maxx = ntrials+2; end;
      set(value(myaxes), 'Xlim', [minx, maxx]);
      drawnow;
      
      % Auto-trainer ----
      
      % increment counters
      psides = value(previous_sides);
    
%       if (length(psides) >1 & length(hit_history) > 0)
%           if(psides(end-1) == 'l') % CR or FA
%             if (hit_history(end) == 1) % CR
% 
%                 if (value(num_consecutive_fas) > 0)
%                     num_consecutive_fas.value = 0;
%                 end
%                 num_consecutive_crs.value = value(num_consecutive_crs)+1;
%             else % FA
%                 if (value(num_consecutive_crs) > 0)
%                     num_consecutive_crs.value = 0;
%                 end
%                 num_consecutive_fas.value = value(num_consecutive_fas)+1;
%             end
%           end
%       end
%       
%       disp(['consecutive FA: ' num2str(value(num_consecutive_fas)) ' CR: '  num2str(value(num_consecutive_crs))]);

      % implement punisher?
%       if(value(autotraining) == 0 & value(autotrain_max_fas) > 0 && ...
%          value(autotrain_max_fas) <= value(num_consecutive_fas))    
%         autotraining.value = 1;
%         num_consecutive_crs.value = 0;
%         disp('STARTING AUTOTRAINER');
%       end
      
      % disable autotrainer?
%       if (value(autotraining) == 1)
%           if (value(num_consecutive_crs) >= value(autotrain_min_crs))
%               autotraining.value = 0;
%               num_consecutive_fas.value = 0;
%               disp('STOPPING AUTOTRAINER');
%           end
%       end
      
      % ---- end autotrainer
      
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
   end;

   
      