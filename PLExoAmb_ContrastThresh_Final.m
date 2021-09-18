function [stair1,stair2,ct1_ToPlot,ct2_ToPlot] = PLExoAmb_ContrastThresh_Final(Observer,Eyetrack,Diag,Block,Eye,stair1,stair2)

% last revised by MR on 8.27.18 
    % -incorporate eyetracking improvements: manual drift correction,
    % correct % correct, keeping track number of recalibrations, break with
    % eyes open after shorter breaks
    % - two interleaving staircases

% contrast thresholding procedure (range .01:.01:1) for 2AFC orientation discrimination from horizontal w/ 1 target and 1
% distractor using only NEUTRAL PRECUES; 4 deg tilted 6 cpd gabor at 4 deg
% eccentricity

% two simultaneous Best PEST staircases to estimate contrast at which observer performing at 75% and 88%  

% Observer: Observer initials, remember to input initials in single quotes,'xx'
% Eyetrack: eyetracking on = 1, eyetracking off = 0,
% Diag: which diagonal this observer will train on (counterbalanced across
    % observers); 1 = upper-right & lower-left quadrants, 2 = upper-left &
    % lower-right quadrants
% Block: number block (1-4)
% Eye: The fellow eye is the left(L) or right(R) eye
% stair1: sets initial staircase values for estimate 1 (88% accuracy) (Block 1 based on uniform prior) that will be transformed
    % into stair after one trial
% stair2: sets initial staircase values for estimate 2 (75% accuracy) (Block 1 based on uniform prior) that will be transformed
    % into stair after one trial



global stimulus;
global staircase1;
global staircase2;
global ct_estimates1; 
global ct_estimates2; 
global task;
global fb_count;
global recalib_count;

thisdir = pwd; % function should be run inside qCSF folder

% make an Observer directory if necessary - 
    datadirname = fullfile(thisdir,'MR_PL_Amb_ExoAtt_data',Observer);
    if ~isdir(datadirname);
        disp(sprintf('Making Observer directory %s',datadirname));
        mkdir(datadirname);
    end
    
% to ensure observer folder is findable    
addpath(genpath('/Users/purplab/Desktop/Mariel/')); 

    disp(sprintf('[ DATA ] saving data in: %s',datadirname));
    

% NEW for interleaved staircases

ct_estimates1 = []; % blank matrix to store trial by trial estimates of ct for plotting; row 1 = 88% stair, row 2 = 75% stair
ct_estimates2 = []; % blank matrix to store trial by trial estimates of ct for plotting; row 1 = 88% stair, row 2 = 75% stair

staircase1 = stair1;
staircase2 = stair2;

% NEW for interleaved staircases

stimulus = [];
stimulus.x_fixcorrect = 0;
stimulus.y_fixcorrect = 0;

fb_count = 0;
recalib_count = 0;

% clearing old variables:
clear task myscreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initalize the screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.EyeTrack=Eyetrack;
myscreen.datadir = datadirname;
myscreen.allowpause = 0;
myscreen.saveData = -2;
myscreen.background=.5;
myscreen.keyboard.nums = [84 85 50]; % 1 & 2 on keyboard; 1 = CCW, 2 = CW

% initalize the screen
myscreen = initScreen(myscreen);
if stimulus.EyeTrack
    myscreen = eyeCalibDisp(myscreen); % function to calibrate eyetracker
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

task{1}.waitForBacktick = 1;

% 1: ITI 2: Wait until fixation 3: drift correction 4: Precue 5: ISI 6: Target + Response Cue
% 7: ISI + Response Cue 8: Response Window + Response Cue (fixation cross
% white and Observer can respond) 9: eyetracking 10: present text for break screen
% 11: present tone at end of break 12: text display to open eyes and wait for next trial
% 13: end of block feedback on screen 14: end of block feedback in command line
task{1}.segmin =     [1 Inf .01 .06 .06 .12 .7 5 1 15 .2 5 5 .1]; % 15 sec break
task{1}.segmax =     [1 Inf .01 .06 .06 .12 .7 5 1 15 .2 5 5 .1];
task{1}.getResponse = [0 0 0 0 0 0 0 1 0 0 0 0]; % Observer can only input response in segment 7; keyboard presses during other segments not recognized

n_repeats = 20; % n_repeats * 2 possible target locations (2) * 2 staircases = 40 trials/staircase

if Diag == 1 % target can only be in UR or LL quadrant during training for that observer
    [location, staircase, repeat] = ndgrid([1,3],[1,2],1:n_repeats); % two training locations, 2 staircases, *n repetitions
else Diag == 2; % target can only be in UL or LR quadrant during training for that observer
     [location, staircase, repeat] = ndgrid([1,3],[1,2],1:n_repeats); % two training locations, 2 staircases, *n repetitions,
end

task{1}.numTrials = length(location(:)); 
random_order = randperm(task{1}.numTrials);
task{1}.randVars.targetLocation = location(random_order);
random_order = randperm(task{1}.numTrials);
task{1}.randVars.staircase = staircase(random_order);
task{1}.randVars.uniform.targetOrientation = 1:2; % can be CW or CCW
task{1}.randVars.uniform.distractorOrientation1 = 1:2; % only one Gabor distractor on each trial
task{1}.randVars.uniform.distractorOrientation2 = 1:2;
task{1}.randVars.uniform.distractorOrientation3 = 1:2;
task{1}.randVars.len_ = task{1}.numTrials;
stimulus.trialend = 0;
stimulus.trialnum = 1;
stimulus.FixationBreak = zeros(1,length(location(:)));
stimulus.Recalibration = zeros(1,length(location(:)));
stimulus.LocationIndices = unique(location);

task{1}.random = 1;
[task{1}, myscreen] = initTask(task{1},myscreen,@StartSegmentCallback,@DrawStimulusCallback,@responseCallback);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myscreen = initStimulus('stimulus',myscreen);
stimulus = myInitStimulus(stimulus,myscreen,task);

% myscreen = eyeCalibDisp(myscreen);
myscreen.eyetracker.savedata = true;%%%%% TO ADD FOR ONLINE EYETRACKING
myscreen.eyetracker.data = [1 1 1 0];%%%%% TO ADD FOR ONLINE EYETRACKING

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
    
    % update the task
    % automatically runs the task, you only need to change: StartSegmentCallback,DrawStimulusCallback,responseCallback
    [task, myscreen, phaseNum] = updateTask(task,myscreen,phaseNum);

    % flip screen
    myscreen = tickScreen(myscreen,task);

end
stair1 = staircase1; % this is the final staircase information output
stair2 = staircase2; % this is the final staircase information output

clear stimulus.tmp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = endTask(myscreen,task);

ct1_ToPlot = ct_estimates1;
ct2_ToPlot = ct_estimates2;
     
    % saves with correct file name according to participant, eye tested,
    % staircase block number
    filename_mat = sprintf('ct_%s_%s_%d.mat',Observer,Eye,Block);

if task{1, 1}.trialnum >= 80 % ensures incomplete blocks do not override full blocks bc same name
    save(fullfile(datadirname,filename_mat)); 
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 1: function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = StartSegmentCallback(task, myscreen)
    global stimulus;
    global fb_count;
    global recalib_count;
    
    if (task.thistrial.thisseg == 1) % ITI
        stimulus.trialend = stimulus.trialend + 1;

    elseif (task.thistrial.thisseg == 2) % fixation
        
 % NEW for automatic drift correction %
       stimulus.eyeParams.fixPos = []; 
 % NEW for automatic drift correction 
 
        task.thistrial.fixationbreak = 0; %%%%% TO ADD FOR ONLINE EYETRACKING %%%
        stimulus.tmp.targetLocation  = stimulus.eccentricity*[stimulus.locations{task.thistrial.targetLocation}];
        stimulus.tmp.distractorIndices=stimulus.LocationIndices(stimulus.LocationIndices~=task.thistrial.targetLocation);
        for Locs=1:length(stimulus.tmp.distractorIndices)
            stimulus.tmp.distractorLocations{Locs}= stimulus.eccentricity*[stimulus.locations{stimulus.tmp.distractorIndices(Locs)}];
        end
        stimulus.FixationStarted=0;

        %response cue
        stimulus.tmp.respcueLocation=stimulus.respcueLocation{task.thistrial.targetLocation};
        stimulus.CueSounded=0;
        
   elseif (task.thistrial.thisseg == 3) % drift correction
        
 % NEW for automatic drift correction %
 
            ep = myscreen.eyetracker.eyepos;
 
       % compute average fixation position from during fixation period
       % column-wise
            avgPos = nanmean(stimulus.eyeParams.fixPos,1); % used by MJ
                        
        % perform automatic drift correction if applicable
            if (sqrt(avgPos(end,1)^2+avgPos(end,2)^2))<=stimulus.TrialStartFixDist && stimulus.FixationStarted
               stimulus.eyeParams.fixRef = avgPos;
               stimulus.x_fixcorrect = -stimulus.eyeParams.fixRef(1);
               stimulus.y_fixcorrect = -stimulus.eyeParams.fixRef(2);

            else % if avg fix position not within the fixation window, then reset the fixRef and fixPos to center
               stimulus.eyeParams.fixRef = [0 0];
               stimulus.eyeParams.fixPos = [];
               stimulus.x_fixcorrect = 0;
               stimulus.y_fixcorrect = 0;
            end  
            
    elseif (task.thistrial.thisseg == 8); % response %%%%
        if ~task.thistrial.gotResponse 
            mglPlaySound(stimulus.CueSound); % play response window start tone
   
     % NEW: for manual calibration
      
        else task.thistrial.gotResponse & task.thistrial.buttonState(3) == 1
            
            task.randVars.fixBreak(stimulus.trialnum) = task.thistrial.fixationbreak;
            stimulus.Recalibration(stimulus.trialnum) = task.thistrial.fixationbreak;

            recalib_count = recalib_count + 1; % add recalibration to counter
            
            stimulus.eyeParams.fixRef = [0 0];
            stimulus.x_fixcorrect = 0;
            stimulus.x_fixcorrect = 0;
               
      % recreate params for fixation break trials and create and add new trial to end of block
        
        %% targetLocation
        if  task.randVars.fixBreak(stimulus.trialnum) == 1;
            temp = zeros(1,size(task.randVars.targetLocation,2)+1);
            temp(1:size(task.randVars.targetLocation,2)) = task.randVars.targetLocation;
            temp(end) = task.randVars.targetLocation(stimulus.trialnum);
            task.randVars.targetLocation = temp;
        end
        
        %% which staircase to run
        if task.randVars.fixBreak(stimulus.trialnum) == 1;
            temp = zeros(1,size(task.randVars.staircase,2)+1);
            temp(1:size(task.randVars.staircase,2)) = task.randVars.staircase;
            temp(end) = task.randVars.staircase(stimulus.trialnum);
            task.randVars.staircase = temp;
        end
        
        %% increase number of trials
        if  task.randVars.fixBreak(stimulus.trialnum) == 1;
            % task.randVars.len_ is length of randomization vector in number of trials for each RandVar, the
            % default is 250 but you can overide in program
            task.randVars.len_ = task.randVars.len_ + 1;
            
            % since in test with fix breaks .numTrials can > default of
            % 250, need to manually tell vector length of ALL randVars (including
            % randVars and randVars.uniform) to increase by one
            task.randVars.varlen_(1) = task.randVars.varlen_(1) + 1;
            task.randVars.varlen_(2) = task.randVars.varlen_(2) + 1;
            task.randVars.varlen_(3) = task.randVars.varlen_(3) + 1;
            task.randVars.varlen_(4) = task.randVars.varlen_(4) + 1;
            task.numTrials = task.numTrials + 1;
            stimulus.FixationBreak = [stimulus.FixationBreak 0]; %%%%%%%%%!!!! Added to increase size of stimulus.FixationBreak matrix each added trial
            stimulus.Recalibration = [stimulus.Recalibration 0]; %%%%%%%%%!!!! Added to increase size of stimulus.Recalibration matrix each added trial

        end
        
     % NEW: for manual calibration

        end
        
    elseif (task.thistrial.thisseg == 9); %%%%% TO ADD FOR ONLINE EYETRACKING
            
       % New for manual calibration
    if task.thistrial.buttonState(3) == 1
        stimulus.starColorFeedback(stimulus.trialnum)=NaN;
    end
       % New for manual calibration
        
        task.randVars.fixBreak(stimulus.trialnum) = task.thistrial.fixationbreak;
        
        fb_count = fb_count + 1; % add fb to counter
        
        if stimulus.trialnum >= 6
            
            recent_fb = sum(task.randVars.fixBreak(stimulus.trialnum-5:stimulus.trialnum-1)); % count how many fix breaks happened in the previous 5 trials 
        
            if fb_count >= 3 && recent_fb >= 3
                
               myscreen = eyeCalibDisp(myscreen); % return back to initial calibration display  
               stimulus.Recalibration(stimulus.trialnum) = task.thistrial.fixationbreak;

               recalib_count = recalib_count + 1; % add recalibration to counter

               fb_count = 0; % reset counter back to 0
               
            end 
            
        end
            
        % recreate params for fixation break trials and create and add new trial to end of block
        
        % targetLocation
        if task.randVars.fixBreak(stimulus.trialnum) == 1;
            temp = zeros(1,size(task.randVars.targetLocation,2)+1);
            temp(1:size(task.randVars.targetLocation,2)) = task.randVars.targetLocation;
            temp(end) = task.randVars.targetLocation(stimulus.trialnum);
            task.randVars.targetLocation = temp;
        end
        
        % which staircase to run
        if task.randVars.fixBreak(stimulus.trialnum) == 1;
            temp = zeros(1,size(task.randVars.staircase,2)+1);
            temp(1:size(task.randVars.staircase,2)) = task.randVars.staircase;
            temp(end) = task.randVars.staircase(stimulus.trialnum);
            task.randVars.staircase = temp;
        end

        %% increase number of trials
        if task.randVars.fixBreak(stimulus.trialnum) == 1;
            task.randVars.len_ = task.randVars.len_ + 1;
            task.randVars.varlen_(1) = task.randVars.varlen_(1) + 1;
            task.randVars.varlen_(2) = task.randVars.varlen_(2) + 1;
            task.randVars.varlen_(3) = task.randVars.varlen_(3) + 1;
            task.numTrials = task.numTrials + 1;
            stimulus.FixationBreak = [stimulus.FixationBreak 0]; %%%%%%%%%!!!! Added to increase size of stimulus.FixationBreak matrix each added trial
        end
        stimulus.trialnum = stimulus.trialnum + 1;

    end

mglClearScreen(stimulus.grayColor);
setGammaTable(1);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % TASK 1: function that gets called to draw the stimulus each frame
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = DrawStimulusCallback(task, myscreen)
global stimulus;
global staircase1;
global staircase2;
global ct_estimates1;
global ct_estimates2;
       
mglClearScreen(stimulus.grayColor); %###
    
if (task.thistrial.thisseg == 1) % ITI
        mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.black); % note that I changed fix cross color from default white to black 
             
elseif (task.thistrial.thisseg == 2) % FIXATION
        mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.black);
        
        if stimulus.EyeTrack
            
                % defines what to call an eyetracker position sample
                ep = myscreen.eyetracker.eyepos;
             
        % NEW: drift correction

            % creates a matrix of eye position samples across the entire
            % fixation period (obtains a new sample every time screen is
            % refreshed at 1000 Hz) 
            stimulus.eyeParams.fixPos = [stimulus.eyeParams.fixPos; ep];       
                
             % if trial fixation distance (calculated as radius or hypotenuse relative to [0,0] is less than size of allowed window
            if (sqrt((ep(end,1)+stimulus.x_fixcorrect)^2+(ep(end,2)+stimulus.y_fixcorrect)^2))<=stimulus.TrialStartFixDist && ~stimulus.FixationStarted
                stimulus.FixationStart=mglGetSecs;
                stimulus.FixationStarted=1;                
                % don't start trial if fixation distance is greater than window
            elseif (sqrt((ep(end,1)+stimulus.x_fixcorrect)^2+(ep(end,2)+stimulus.y_fixcorrect)^2))>=stimulus.TrialStartFixDist && stimulus.FixationStarted
                stimulus.FixationStarted=0;
                % if the elapsed time since the trial started is >250 ms
                % and fixation distance is still less than window size jump
                % to next segment
            elseif (sqrt((ep(end,1)+stimulus.x_fixcorrect)^2+(ep(end,2)+stimulus.y_fixcorrect)^2)) <= stimulus.TrialStartFixDist
                stimulus.FixationDur=mglGetSecs(stimulus.FixationStart);
                if stimulus.FixationDur >=stimulus.TrialStartFixDur
                    task = jumpSegment(task);
                end
            end
        else  % if eyetracking is off
                mglWaitSecs(.24) % waits 240 ms before moving into drift correction period for last 10 ms
                task = jumpSegment(task);
        end
        
elseif (task.thistrial.thisseg == 3) % DRIFT CORRECTION 
    
    mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.black);
       
elseif (task.thistrial.thisseg == 4) % Neutral cue
    mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.black);
            
    % precue for neutral training group 

           mglGluDisk(0,1,stimulus.PrecueSize,stimulus.white); % supposed to draw white filled dot of .16 size .5 deg from edge of cross (.5 tall total, so .25 from center of screen)
           mglGluDisk(0,-1,stimulus.PrecueSize,stimulus.white); % draw identical white dot below cross
           mglGluDisk(1,0,stimulus.PrecueSize,stimulus.white); % draw identical white dot to right of cross
           mglGluDisk(-1,0,stimulus.PrecueSize,stimulus.white); % draw identical white dot to left of cross
           
    if stimulus.EyeTrack;%%%%% TO ADD FOR ONLINE EYETRACKING
        ep=myscreen.eyetracker.eyepos;
        if (sqrt((ep(end,1)+stimulus.x_fixcorrect)^2+(ep(end,2)+stimulus.y_fixcorrect)^2))>stimulus.TrialStartFixDist
           stimulus.FixationBreak(stimulus.trialnum)=1;
           task.thistrial.fixationbreak = 1;
           task = jumpSegment(task);
        end
    end
                
elseif (task.thistrial.thisseg == 5) % ISI
       mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.black);
                
       if stimulus.EyeTrack%%%%% TO ADD FOR ONLINE EYETRACKING
          ep=myscreen.eyetracker.eyepos;
          if task.thistrial.fixationbreak == 1
             task = jumpSegment(task);
          end;
          if (sqrt((ep(end,1)+stimulus.x_fixcorrect)^2+(ep(end,2)+stimulus.y_fixcorrect)^2))>stimulus.TrialStartFixDist
             stimulus.FixationBreak(stimulus.trialnum)=1;
             task.thistrial.fixationbreak = 1;
             task = jumpSegment(task);
          end
       end
                
elseif (task.thistrial.thisseg == 6) % GABOR STIMULI & RESPONSE CUE
       mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.black);
       
       % DRAWS RESPONSE CUE
       mglLines2(stimulus.tmp.respcueLocation(1), stimulus.tmp.respcueLocation(2), stimulus.tmp.respcueLocation(3), stimulus.tmp.respcueLocation(4),2.5,stimulus.white);
   
       if stimulus.EyeTrack%%%%% TO ADD FOR ONLINE EYETRACKING
          ep=myscreen.eyetracker.eyepos;
          if task.thistrial.fixationbreak == 1
             task = jumpSegment(task);
          end;
          if (sqrt((ep(end,1)+stimulus.x_fixcorrect)^2+(ep(end,2)+stimulus.y_fixcorrect)^2))>stimulus.TrialStartFixDist
           stimulus.FixationBreak(stimulus.trialnum)=1;
           task.thistrial.fixationbreak = 1;
           task = jumpSegment(task);
          end
       end
    
   % NEW for interleaved staircases  
       
     %%% draws target gabor - updates stim contrast on
     %%% every trial according to which staircase is being run 
     
     if task.thistrial.staircase == 1
        drawGabor(staircase1.xCurrent,stimulus.tmp.targetLocation, stimulus.rotation(task.thistrial.targetOrientation), 1); 
        
         % draws distractor gabor
            for Loc=1:length(stimulus.tmp.distractorIndices)
                 eval(sprintf('drawGabor(staircase1.xCurrent,stimulus.tmp.distractorLocations{Loc}, stimulus.rotation(task.thistrial.distractorOrientation%g), 1);',Loc));
            end
             
     else task.thistrial.staircase == 2;
        drawGabor(staircase2.xCurrent,stimulus.tmp.targetLocation, stimulus.rotation(task.thistrial.targetOrientation), 1);
        
            for Loc=1:length(stimulus.tmp.distractorIndices)
                 eval(sprintf('drawGabor(staircase2.xCurrent,stimulus.tmp.distractorLocations{Loc}, stimulus.rotation(task.thistrial.distractorOrientation%g), 1);',Loc));
            end
     end
     
   % NEW for interleaved staircases  
                
elseif (task.thistrial.thisseg == 7) % ISI-2
       mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.black);
       
       % DRAWS RESPONSE CUE
       mglLines2(stimulus.tmp.respcueLocation(1), stimulus.tmp.respcueLocation(2), stimulus.tmp.respcueLocation(3), stimulus.tmp.respcueLocation(4),2.5,stimulus.white);
       
       if task.thistrial.fixationbreak == 1
          task = jumpSegment(task);
       end;
       
elseif (task.thistrial.thisseg == 8) % RESPONSE WINDOW
    
       mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);

       % DRAWS RESPONSE CUE
       mglLines2(stimulus.tmp.respcueLocation(1), stimulus.tmp.respcueLocation(2), stimulus.tmp.respcueLocation(3), stimulus.tmp.respcueLocation(4),2.5,stimulus.white);
       
       % RESPONSE WINDOW TONE
       if ~stimulus.CueSounded
          mglClearScreen(stimulus.grayColor);
          mglPlaySound(stimulus.CueSound);
          stimulus.CueSounded=1;
       end
       
       % IF OBSERVER INPUT RESPONSE, SKIP TO NEXT SEGMENT
       if task.thistrial.gotResponse || task.thistrial.buttonState(3) == 1

           mglClearScreen(stimulus.grayColor);
           
           if task.thistrial.staircase == 1   
               % update vector of ct estimates as function of trial for this trial
               ct_estimates1 = [ct_estimates1 staircase1.xCurrent];
               
           else task.thistrial.staircase == 2;
               ct_estimates2 = [ct_estimates2 staircase2.xCurrent];

           end
           
          task = jumpSegment(task);
          
       end;

       % IF OBSERVER HAD FIXATION BREAK, SKIP TO NEXT SEGMENT
       if task.thistrial.fixationbreak
          mglClearScreen(stimulus.grayColor);
          task = jumpSegment(task);
       end
                
elseif (task.thistrial.thisseg == 9) %%%%% TO ADD FOR ONLINE EYETRACKING
    
     % NEW: for manual calibration
       if ~task.thistrial.fixationbreak | task.thistrial.buttonState(3) == 1
          mglClearScreen(stimulus.grayColor);
          task = jumpSegment(task);
       end
     % NEW: for manual calibration  
     
       if task.thistrial.fixationbreak & task.thistrial.buttonState(3) == 0
          mglTextSet('Courier',50,stimulus.black);
          mglTextDraw('Please fixate',[0 0]);
       end
       
elseif (task.thistrial.thisseg == 10) % different points through block, text will appear onscreen to ask observer to take break
    
    if mod(task.numTrials,2) == 0
    
        halfblock = round((task.numTrials/2))+1;
        
        if stimulus.trialnum == halfblock && ~task.thistrial.fixationbreak
           mglTextSet('Courier',50,stimulus.black);
           mglTextDraw('Please take a break and close your eyes',[0,0]);
        else
            task = jumpSegment(task); % if it's not halfway through block, skip this segment
        end
        
    else mod(task.numTrials,2) == 1;
        
        halfblock = round((task.numTrials/2))+1;
        
        if stimulus.trialnum == halfblock && ~task.thistrial.fixationbreak 
           mglTextSet('Courier',50,stimulus.black);
           mglTextDraw('Please take a break and close your eyes',[0,0]);
        else
            task = jumpSegment(task); % if it's not halfway through block, skip this segment
        end 
    end
    
elseif (task.thistrial.thisseg == 11)  % after break ends, tone is presented so observer opens eyes
    
    if mod(task.numTrials,2) == 0 % with fixation breaks the total num trials is even

        halfblock = round((task.numTrials/2))+1;
        
        if stimulus.trialnum == halfblock && ~task.thistrial.fixationbreak
                mglPlaySound(stimulus.CueSound);
        else
            task = jumpSegment(task); % if it's not halfway through block, skip this segment
        end
        
   else mod(task.numTrials,2) == 1; % with fixation breaks the total num trials is odd

        halfblock = round((task.numTrials/2))+1;

        if stimulus.trialnum == halfblock && ~task.thistrial.fixationbreak
                mglPlaySound(stimulus.CueSound);
        else
            task = jumpSegment(task); % if it's not halfway through block, skip this segment
        end
    end
    
elseif (task.thistrial.thisseg == 12) % ask observer to keep eyes open to adjust to illumination of screen
    
   % NEW for manual calibration

        halfblock = round((task.numTrials/2))+1;
    
        if stimulus.trialnum == halfblock && ~task.thistrial.fixationbreak
           mglTextSet('Courier',50,stimulus.black);
           mglTextDraw('Please keep your eyes open',[0,0]);
        else
            task = jumpSegment(task); % if it's not halfway through block, skip this segment
        end
        
    % NEW for manual calibration
        
       
elseif (task.thistrial.thisseg == 13) % end of block feedback on screen for observer 

    if stimulus.trialnum<=task.numTrials % End of block Feedback
        task = jumpSegment(task);
    else
        
% NEW CORRECT CALCULATION OF % CORRECT
        recalib_index = find(stimulus.Recalibration == 1);
        stimulus.FixationBreak(recalib_index) = NaN; % replace 1s for NaN where participant recalibrated but did not break fixation
        
        CorrectTrials = find(stimulus.starColorFeedback==2);
        CorrectCount = length(stimulus.starColorFeedback(CorrectTrials));
        TotalTrials_ind = find(stimulus.starColorFeedback==2 |stimulus.starColorFeedback==3);
        TotalTrials = length(stimulus.starColorFeedback(TotalTrials_ind));
     
        PercentCorrect=([num2str(round(100*(CorrectCount/(TotalTrials)))) num2str('%')]);
        NumFixBreak=num2str(nansum(stimulus.FixationBreak));
      
        Recalibs=num2str(nansum(stimulus.Recalibration));
      
        mglTextSet('Courier',50,stimulus.black);
        eval(sprintf('mglTextDraw(''%s correct'',[0 1]);',PercentCorrect));
        eval(sprintf('mglTextDraw(''%s fixation breaks'',[0 0]);',NumFixBreak));
        eval(sprintf('mglTextDraw(''%s recalibration(s)'',[0 -1]);',Recalibs));
        
% NEW CORRECT CALCULATION OF % CORRECT
    end   

else (task.thistrial.thisseg == 14) % End of block Feedback presented to experimenter in command window
 
    if stimulus.trialnum<=task.numTrials % 
        task = jumpSegment(task);
    else   
        
% NEW for manual calibration
   
        CorrectTrials = find(stimulus.starColorFeedback==2);
        CorrectCount = length(stimulus.starColorFeedback(CorrectTrials));
        TotalTrials_ind = find(stimulus.starColorFeedback==2 |stimulus.starColorFeedback==3);
        TotalTrials = length(stimulus.starColorFeedback(TotalTrials_ind));
     
        PercentCorrect=([num2str(round(100*(CorrectCount/(TotalTrials)))) num2str('%')]);
        NumFixBreak=num2str(nansum(stimulus.FixationBreak));
      
        Recalibs=num2str(nansum(stimulus.Recalibration));
      
        % fix to make accurate to account for fixation breaks
        disp(sprintf('Observer performance was %s correct',PercentCorrect));
        disp(sprintf('Observer broke fixation %s times',NumFixBreak));
        disp(sprintf('Observer recalibrated %s times',Recalibs)); % note that this counts the number of times recalibration screen appeared, but participant may have opted out of recalibration

% NEW for manual calibration
        
    end
   
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % function to get the Observer's response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = responseCallback(task, myscreen)
global stimulus;
global staircase1;
global staircase2; 

mglClearScreen(stimulus.grayColor);
if ~task.thistrial.gotResponse
    
  % NEW: add manual calibration option at response window
  
  if task.thistrial.buttonState(1) == 1 || task.thistrial.buttonState(2) == 1

    % check response correct or not
    stimulus.tmp.response = task.thistrial.whichButton == task.thistrial.targetOrientation; %1 for CCW and 2 for CW

    % if correct
    if stimulus.tmp.response
       mglPlaySound(stimulus.CorrectSound);
       stimulus.starColorFeedback(stimulus.trialnum)=2;
          % NEW for simultaneous staircases
          if task.thistrial.staircase == 1
            %%%%%%%%%%% update staircase1 according to whether trial correct or incorrect
            staircase1 = usePalamedesStaircase(staircase1,stimulus.tmp.response);
            %%%%%%%%%%%  
          else task.thistrial.staircase == 2;
            %%%%%%%%%%% update staircase2 according to whether trial correct or incorrect
            staircase2 = usePalamedesStaircase(staircase2,stimulus.tmp.response);
            %%%%%%%%%%%   
          end
       task = jumpSegment(task);
    
    else
       mglPlaySound(stimulus.IncorrectSound);
       stimulus.starColorFeedback(stimulus.trialnum)=3;
           % NEW for simultaneous staircases
          if task.thistrial.staircase == 1
            %%%%%%%%%%% update staircase1 according to whether trial correct or incorrect
            staircase1 = usePalamedesStaircase(staircase1,stimulus.tmp.response);
            %%%%%%%%%%%  
          else task.thistrial.staircase == 2;
            %%%%%%%%%%% update staircase2 according to whether trial correct or incorrect
            staircase2 = usePalamedesStaircase(staircase2,stimulus.tmp.response);
            %%%%%%%%%%%   
          end
       task = jumpSegment(task);
    end
  
  else task.thistrial.buttonState(3) == 1
       mglClearScreen(stimulus.grayColor);
     
       myscreen = eyeCalibDisp(myscreen);
       
       % New for manual calibration
        stimulus.starColorFeedback(stimulus.trialnum)=NaN;
      
        stimulus.FixationBreak(stimulus.trialnum)=1;
        stimulus.Recalibration(stimulus.trialnum)=1;

        task.thistrial.fixationbreak = 1;

        % New for manual calibration
                           
        task = jumpSegment(task);
        mglClearScreen(stimulus.grayColor);

   end
        
   % NEW: add manual calibration option at response window %  
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to draw the gabor stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawGabor(desiredContrast,position,orientation,sf)
% drawGaborPedCdeltaC
%
%        $Id: drawGabor.m, v 1 2007/01/18 19:40:56 ?? ?? $
%      usage: drawGabor(desiredContrast,position,orientation,sf)
%    purpose: draw a gabor stimulus on the screen with a specified contrast
%             within a given clut range (it finds the closest contrast to
%             the requested one within the available clut range)

global stimulus;

% now find closest matching contrast we can display with this gamma table
displayContrastNum = min(round(stimulus.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast),stimulus.nDisplayContrasts);
% disp(sprintf('Desired contrast: %0.4f Actual contrast: %0.4f Difference: %0.4f',desiredContrast,stimulus.currentMaxContrast*(displayContrastNum/stimulus.nDisplayContrasts),desiredContrast-stimulus.currentMaxContrast*(displayContrastNum/stimulus.nDisplayContrasts)));
    if round(stimulus.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast)>stimulus.nDisplayContrasts
        disp(sprintf('[drawGabor] Desired contrast out of range %0.2f > %0.2f',desiredContrast,stimulus.currentMaxContrast));
        keyboard
    end

mglBltTexture(stimulus.tex{sf}(displayContrastNum+1),position,0,0,orientation); %mglBltTexture(texture,position,hAlignment,vAlignment,rotation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to create a gamma table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setGammaTable(maxContrast)
global stimulus;

% set the reserved colors
gammaTable(1:size(stimulus.reservedColors,1),1:size(stimulus.reservedColors,2))=stimulus.reservedColors;
% create the gamma table
cmax = 0.5+maxContrast/2;cmin = 0.5-maxContrast/2;
luminanceVals = cmin:((cmax-cmin)/(stimulus.nGratingColors-1)):cmax;
% now get the linearized range
redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');

% check the table
% plot(stimulus.linearizedGammaTable.redTable,'k');
% hold on
% plot(256*(0.25:0.5/250:0.75),redLinearized,'r');

% set the gamma table
gammaTable((stimulus.minGratingColors+1):256,:)=[redLinearized;greenLinearized;blueLinearized]';
% set the gamma table
mglSetGammaTable(gammaTable);
% remember what the current maximum contrast is that we can display
stimulus.currentMaxContrast = maxContrast;
                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task)
global MGL;


% let's get the linearized gamma table
stimulus.linearizedGammaTable = mglGetGammaTable;
stimulus.linearizedGammaTable.redTable(1:3) = 0; % this is just to provisionally deal with what appears to be some bug: the first value in each of these gamma tables is a NaN
stimulus.linearizedGammaTable.greenTable(1:3) = 0;
stimulus.linearizedGammaTable.blueTable(1:3) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stimulus parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stimulus.init = 1;

% gabors
stimulus.width = 4*.8;%stimulus.gaussSdx*7 = 3.2;             % in deg
stimulus.height = 4*.8;%stimulus.gaussSdy*7;            % in deg
stimulus.gaussSdx = stimulus.width/7;                % in deg
stimulus.gaussSdy = stimulus.height/7;               % in deg
stimulus.Tilt = 4;
stimulus.rotation = [stimulus.Tilt -stimulus.Tilt];
stimulus.DistractorRotation = 0; % this is the tilt orientation of the gabor stimulus from vertical in Degrees
stimulus.sf = 6;   % in cpd
stimulus.orientation = 0; % 90 = horizontal ref angle?; stimulus.orientation = 0 is vertical; % in deg
stimulus.phase = 0; % in deg
stimulus.eccentricity = 4;  % in deg 


% fixation
stimulus.FCwidth = 1.0; %stimulus.FCwidth = 0.5;
stimulus.FClinewidth = 1.5;
stimulus.TrialStartFixDist= 1; % 1 degree radius in which to fixate before trial starts
stimulus.TrialStartFixDur=.25;
stimulus.cornerDist=(stimulus.width/2); % used to calculate placeholder distance 
stimulus.placeHolderSize=.05;
stimulus.PrecueSize=.16;
stimulus.ExoCueDist=(stimulus.width/2)+.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frames and locations:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.locations = {[cosd(1),cosd(1)],[cosd(1), -cosd(1)],[-cosd(1), -cosd(1)],[-cosd(1), cosd(1)]};  %4 locations: x,y coordinates are specified here. starts in N and moves clockwise

% stimulus.locations ={[1,1],[1,-1],[-1,-1],[-1,1]}; places stimuli and
% response cue on diagonals

% uncomment if add placeholders back in
%     for i=1:length(stimulus.locations)
%         stimulus.placeholders{i}= [stimulus.eccentricity*stimulus.locations{i}];
%         stimulus.placeholders{i}(2,:)=[stimulus.placeholders{i}(1,1:2)]+[0 stimulus.cornerDist];
%         stimulus.placeholders{i}(3,:)=[stimulus.placeholders{i}(1,1:2)]+[stimulus.cornerDist 0];
%         stimulus.placeholders{i}(4,:)=[stimulus.placeholders{i}(1,1:2)]+[0 -stimulus.cornerDist];
%         stimulus.placeholders{i}(5,:)=[stimulus.placeholders{i}(1,1:2)]+[-stimulus.cornerDist 0];
%     end
    
% Note: locations should specify a vector of length 1, i.e.,
%                       sum(locations.^2)==1
stimulus.frameThick = .08;
stimulus.reservedColors = [0 0 0; 1 1 1; 0 .6 0];

% Response cue
stimulus.rcue.XLocation{1} = [-.6; 0; -0.4; 0];
stimulus.rcue.YLocation{1} = [ 0; .6; 0; -.6];
stimulus.rcue.XLocation{2} = [.7; 0; 0.4; 0];
stimulus.rcue.YLocation{2} = [ 0; .6; 0; -.6];
stimulus.rcue.color        = [0; .6; 0]; % green
row=[4 5 2 3];
signs{1}=[-1 -1];signs{2}=[-1 1];signs{3}=[1 1];signs{4}=[1 -1];
    for i=1:length(stimulus.LocationIndices)
%         stimulus.EndocueLocation{i}(1:2)=.8*[stimulus.locations{i}];
%         stimulus.EndocueLocation{i}(3:4)=1.6*[stimulus.locations{i}];
%         stimulus.NeutralcueLocation{i}(1:2)=.8*[stimulus.locations{i}];
%         stimulus.NeutralcueLocation{i}(3:4)=1.2*[stimulus.locations{i}];
        stimulus.ExocueLocation{i}=(stimulus.eccentricity+stimulus.ExoCueDist)*stimulus.locations{i};
    end
    
% for i=1:length(stimulus.LocationIndices)
for i = 1:length(stimulus.locations)
    stimulus.respcueLocation{i}(1:2)=.8*[stimulus.locations{i}];
    stimulus.respcueLocation{i}(3:4)=1.6*[stimulus.locations{i}];
end

%% !!! not sure about this section
% stimulus.ExocueLocation{length(stimulus.LocationIndices)+1}=[0 0];
% stimulus.respcueLocation{i}= [1 .5;1 .7]; %[0.8785 0;1.1714 0];
% stimulus.respcueLocation{2}= [-.5 -1;-.7 -1];%[-0.8785 0;-1.1714 0];
% stimulus.respcueLocation{3}= [1 -.5;1 -.7];%[0.8785 0;1.1714 0];
% stimulus.respcueLocation{4}= [-.5 1;-.7 1];%[-0.8785 0;-1.1714 0];

% tentative_threshold = 0.1;
% log_steps = 1; % i.e., step=1/2 means every 2 steps is equal to doubling the contrast. step=1 means every step doubles the contrast, step=1/3 means every 3 steps double the contrast!
% % step and number of contrast levels determine the range. We want a range
% % that covers very low (i.e., performance is at chance) and very high
% % (i.e., performance at ceiling).
%
% % % calculate a logarithmic scale, centered around the average threshold
% stimulus.contrasts = 2.^((-3:3)*log_steps)*tentative_threshold;
% MG changed: stimulus.contrasts =[0.10 0.20 0.30 0.40 0.60 0.80 0.90];
%stimulus.contrasts =.3; %changed to use 1 contrast level

% some computations:
% for contrast
% idea: for low contrast stimuli, take advantage of the color-resolution of the
% colormap (10 bits) instead of the limited resolution of the grayspace (8 bit)
% that is available. This minimized the quatization error (causes
% artifacts, in particular in the periphery of the Gabor) by a factor of 4
% (2 bits). The same idea is achievable by dithering (adding random noise)
stimulus.nReservedColors=size(stimulus.reservedColors,1);
stimulus.nGratingColors = 256-(2*floor(stimulus.nReservedColors/2)+1);
stimulus.minGratingColors = 2*floor(stimulus.nReservedColors/2)+1;
stimulus.midGratingColors = stimulus.minGratingColors+floor(stimulus.nGratingColors/2);
stimulus.maxGratingColors = 255;
stimulus.deltaGratingColors = floor(stimulus.nGratingColors/2);

% to set up color values
stimulus.black = [0 0 0];
stimulus.white = [1/255 1/255 1/255];
stimulus.green = [0 6/255 0];
stimulus.grey = [.025 .025 .025];
stimulus.background = [stimulus.midGratingColors/255 stimulus.midGratingColors/255 stimulus.midGratingColors/255];
stimulus.fixation.color = [0; .6; 0]'; % green

% calculate a grating, a gaussian envelope (gaussian is in the alpha
% channel), and a mask (for now, just high-contrast random noise)
    for thisSF = 1:length(stimulus.sf)      %only one spatial frequency
        gratingMatrix{thisSF} = mglMakeGrating(stimulus.width,stimulus.height,stimulus.sf(thisSF),stimulus.orientation,stimulus.phase);
    end

grating(:,:,4) = 255*mglMakeGaussian(stimulus.width,stimulus.height,stimulus.gaussSdx,stimulus.gaussSdy);

% making the texture for all the Gabor stimuli:
disppercent(-inf,'Calculating gabors');
    for thisSF = 1:length(stimulus.sf)
        for thisContrast = 0:stimulus.deltaGratingColors
        %% stimulus.texture
            grating(:,:,1) = stimulus.midGratingColors+gratingMatrix{thisSF}*thisContrast;
            grating(:,:,2) = grating(:,:,1);
            grating(:,:,3) = grating(:,:,1);
            stimulus.tex{thisSF}(thisContrast+1) = mglCreateTexture(grating);
            disppercent(thisContrast/stimulus.deltaGratingColors);
        end
    end
disppercent(inf);
stimulus.nDisplayContrasts = stimulus.deltaGratingColors;
disppercent(inf);


% calculate gray color
stimulus.grayColor = stimulus.background; %stimulus.midGratingColors/255;

% sounds
stimulus.CueSound = find(strcmp(MGL.soundNames,'Ping'));
stimulus.CorrectSound = find(strcmp(MGL.soundNames,'Submarine'));
stimulus.IncorrectSound = find(strcmp(MGL.soundNames,'Basso'));
