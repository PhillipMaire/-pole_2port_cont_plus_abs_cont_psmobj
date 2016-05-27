% [x, y] = AnalysisSection(obj, action, x, y)
%
% For doing online analysis of behavior.
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

function [x, y] = AnalysisSection(obj, action, x, y)

GetSoloFunctionArgs;

switch action

    case 'init',   % ------------ CASE INIT ----------------
        % Save the figure and the position in the figure where we are
        % going to start adding GUI elements:
        SoloParamHandle(obj, 'my_gui_info', 'value', [x y gcf]); next_row(y,1.5);


        MenuParam(obj, 'analysis_show', {'view', 'hide'}, 'hide', x, y, 'label', 'Analysis', 'TooltipString', 'Online behavior analysis');
        set_callback(analysis_show, {mfilename,'hide_show'});

        next_row(y);
        SubheaderParam(obj, 'sectiontitle', 'Analysis', x, y);

        parentfig_x = x; parentfig_y = y;

        % ---  Make new window for online analysis
        SoloParamHandle(obj, 'analysisfig', 'saveable', 0);
        analysisfig.value = figure('Position', [875 50 400 330], 'Menubar', 'none',...
            'Toolbar', 'none','Name','Analysis','NumberTitle','off');

        x = 1; y = 1;
        
        % text stuff
        DispParam(obj, 'DprimeEQ',  '0.000  0.000  n/a', x, y, 'TooltipString', ...
            'D-prime considering only extremal quartiles of position range'); next_row(y);
        DispParam(obj, 'Dprime60',  '0.000  0.000  n/a', x, y, 'TooltipString', ...
            'D-prime considering last 60 trials only');  next_row(y);
        DispParam(obj, 'Dprime',    '0.000  0.000  n/a', x, y); next_row(y);
        DispParam(obj, 'PctCorrect','000.0  000.0  000.0', x, y); next_row(y);
        DispParam(obj, 'NumIgnores', '0000   0000   0000', x, y); next_row(y);
        DispParam(obj, 'NumRewards','0000   0000   0000', x, y); next_row(y);
        DispParam(obj, 'NumTrials', '0000   0000   0000', x, y); next_row(y);
        DispParam(obj, 'Label',     '  L      R    Tot  ' ,x, y); next_row(y);
        
        % Figure with position-dependent %correct
        pos = get(gcf, 'Position');
        SoloParamHandle(obj, 'poscoraxes', 'saveable', 0, 'value', axes);
        set(value(poscoraxes), 'Units', 'pixels');
        set(value(poscoraxes), 'Position', [50 pos(4)-110 pos(3)-60 100]);
        set(value(poscoraxes), 'YTick', [0 0.5 1], 'YLim', [0 1], 'YTickLabel', ...
                        {'0', '0.5', '1'});
        set(value(poscoraxes), 'XTick', 1:2:17, 'XLim', [0 18], 'XTickLabel', ...
                        {'1','3','5','7','9','11','13','15','17'});  
        xlabel('position ctr (10Ks of Zaber units; 20K bins)');
        ylabel('frac correct');
        SoloParamHandle(obj, 'previous_poscor_plot', 'saveable', 0);

        AnalysisSection(obj,'hide_show');
                
                
        x = parentfig_x; y = parentfig_y;
        set(0,'CurrentFigure',value(myfig));
        return;

   
    case 'update' 
        % --- gather relevant data
        correct = hit_history;
        nonIgnores = find(correct >= 0);

        % ALL trials
        sL = find(previous_sides(1:end-1) == 108); % 108 is char for l
        sR = find(previous_sides(1:end-1) == 114); % 114 is char for r
        sLNI = find(previous_sides(1:end-1) == 108 & correct' >= 0 ); % 108 is char for l
        sRNI = find(previous_sides(1:end-1) == 114 & correct' >= 0 ); % 114 r
        
        % only consider trials w/ respnoses
        if (length(nonIgnores) >= 61)
            nI60 = nonIgnores((end-60):(end-1));
            sL60 = find(previous_sides(nI60) == 108 ...
                & correct(nI60)' >= 0 ); % 108 is char for l
            sR60 = find(previous_sides(nI60) == 114 & ...
                correct(nI60)' >= 0); % 108 is char for 4
            sL60 = nI60(sL60);
            sR60 = nI60(sR60);
        else
            sL60 = [] ;sR60 = [];
        end
        
        % quartile estimator -- assumption is that each extremal quartile
        %  is for one lickport.  Only consider trials w/ responses
        posrange = 0;
        sQ1 = []; sQ4= [];
        if (length(pole_position_history) > 0 && range(pole_position_history) > 0)
          % ORIGINAL
          %posrange = range(pole_position_history);
          %posrangeq1 = min(pole_position_history)*[1 1] + [0 posrange/4];
          %posrangeq4 = posrangeq1(1)*[1 1] + [.75 1]*posrange;
          
          % new -- just use 25th and 75th quantile QED
          q25 = quantile(pole_position_history,.25);
          q75 = quantile(pole_position_history,.75);
          posrangeq1 =  min(pole_position_history)*[1 1] + [0 q25];
          posrangeq4 = [q75 max(pole_position_history)];
          
          eq1i = find(pole_position_history >= posrangeq1(1) & ...
                      pole_position_history <= posrangeq1(2) & ...
                      correct >= 0);
          eq4i = find(pole_position_history >= posrangeq4(1) & ...
                      pole_position_history <= posrangeq4(2) & ...
                      correct >=0);
          disp(['EQ posrange: ' num2str(posrangeq1) ' and ' num2str(posrangeq4)]);
          if (length(eq1i) > 1 && length(eq4i) > 1)
              sQ1 = eq1i;
              sQ4 = eq4i;
          end
        end
        
        % --- compute parameters
        
        % # trials
        nt =     [length(sL) ...
                  length(sR) ...
                  length(previous_sides(1:end-1))];
        ntNI =     [length(sLNI) ...
                  length(sRNI) ...
                  length(nonIgnores)];
        nt60 = [length(sL60) length(sR60)];
        ntEQ = [length(sQ4) length(sQ1)];
        
        % # rewards
        nr = [length(find(correct(sL) ==1)) length(find(correct(sR) ==1)) ...
              length(find(correct == 1))];
        nr60 = [length(find(correct(sL60) ==1)) length(find(correct(sR60) ==1))];
        nrEQ = [length(find(correct(sQ4) ==1)) length(find(correct(sQ1) ==1))];
         
        % # incorrects (DISTINCT from cases where he did not respond)
        ni = [length(find(correct(sL) ==0)) length(find(correct(sR) ==0)) ...
              length(find(correct == 0))];
        nig = [length(find(correct(sL) ==-1)) length(find(correct(sR) ==-1)) ...
              length(find(correct == -1))];
        ni60 = [length(find(correct(sL60) ==0)) length(find(correct(sR60) ==0))];
        niEQ = [length(find(correct(sQ4) ==0)) length(find(correct(sQ1) ==0))];
        
        % %correct
        pc= 100*nr./ntNI;
        
        % D-prime REGULAR
        
        % dprime (hit rate, false-alarm rate, # stim in hit/miss pos, #
        %          stim cr/fa pos) 
        % NOTE: we EXCLUDE ignore trials from d' calculation
        dpL = dprime(nr(1)/ntNI(1), ni(2)/ntNI(2), ntNI(1), ntNI(2));
        dpR = dprime(nr(2)/ntNI(2), ni(1)/ntNI(1), ntNI(2), ntNI(1));
        
        % D-prime last 60
        dpL60 = dprime(nr60(1)/nt60(1), ni60(2)/nt60(2), nt60(1), nt60(2));
        dpR60 = dprime(nr60(2)/nt60(2), ni60(1)/nt60(1), nt60(2), nt60(1));       
        
        % D-prime extreme quartiles
        dpLEQ = dprime(nrEQ(1)/ntEQ(1), niEQ(2)/ntEQ(2), ntEQ(1), ntEQ(2));
        dpREQ = dprime(nrEQ(2)/ntEQ(2), niEQ(1)/ntEQ(1), ntEQ(2), ntEQ(1));      
        
        % %correct in position bins
        posbin_count = zeros(9,1);
        posbin_frac_correct = nan*zeros(9,1);
        posbin = 0:20000:160000;        
        posbin_xvals = 1:2:17;
        for pb=1:9
            if (pb == length(posbin))
                idx = find(pole_position_history >= posbin(pb) & ...
                           pole_position_history <= posbin(pb)+20000);
            else
                idx = find(pole_position_history >= posbin(pb) & ...
                           pole_position_history < posbin(pb)+20000);
            end
            posbin_count(pb) = length(idx);
            if (length(idx)>0)
              n_correct = length(find(hit_history(idx) == 1));
              posbin_frac_correct(pb)= n_correct/length(idx);
            end
        end
        
        % --- update plot
        if ~isempty(value(previous_poscor_plot)), delete(previous_poscor_plot(:)); end;
        if isempty(previous_sides), return; end;

        % BLUE dots
        xvals = posbin_xvals;
        yvals = posbin_frac_correct;
        
        hb = line(xvals,yvals, 'Parent', value(poscoraxes));
        set(hb, 'Color', 'b', 'Marker', '.', 'LineStyle', 'none');
      %  arange = axis;
     %   
    %    plot([arange(1) arange(2)], [0.5 0.5], 'k:');
      
        previous_poscor_plot.value = [hb];

        drawnow;
      
        % --- update strings
        NumTrials.value = sprintf('%04d  %04d  %04d', nt(1), nt(2), nt(3));
        NumIgnores.value = sprintf('%04d  %04d  %04d', nig(1), nig(2), nig(3));   
        NumRewards.value = sprintf('%04d  %04d  %04d', nr(1), nr(2), nr(3));   
        PctCorrect.value = sprintf('%05.1f  %05.1f  %05.1f', pc(1), pc(2), pc(3));   
        Dprime.value = sprintf('%05.3f  %05.3f  n/a', dpL, dpR);
        Dprime60.value = sprintf('%05.3f  %05.3f  n/a', dpL60, dpR60);
        DprimeEQ.value = sprintf('%05.3f  %05.3f  n/a', dpLEQ, dpREQ);
        %---------------------------------------------------
 
        
    case 'hide_show'
        if strcmpi(value(analysis_show), 'hide')
            set(value(analysisfig), 'Visible', 'off');
        elseif strcmpi(value(analysis_show),'view')
            set(value(analysisfig),'Visible','on');
            AnalysisSection(obj, 'update_plot');
        end;
        return;


    case 'reinit'
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


