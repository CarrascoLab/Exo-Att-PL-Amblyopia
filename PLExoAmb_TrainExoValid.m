function myscreen = PLExoAmb_TrainExoValid(Observer,Eyetrack,Contrasts,Diag,Day,Block)

% 2AFC training task of orientation discrimination 4 deg tilt from horizontal w/ 1 target and 1
% distractor using VALID EXOGENOUS PERIPHERAL PRECUES

% Observer: Observer initials, remember to input initials in single quotes,'xx'
% Eyetrack: eyetracking on = 1, eyetracking off = 0,
% Contrasts: individual's 2 contrast levels [c50 rmax] in 1x2 vector
% Diag: which diagonal this observer will train on (counterbalanced across
% observers); 1 = upper-right & lower-left quadrants, 2 = upper-left &
% lower-right quadrants
% Day: number training day (1-10)
% Block: number training block (1-8)

global stimulus;
global MGL;
global fb_count;

fb_count = 0;

thisdir = pwd;

% make an Observer directory if necessary
    datadirname = fullfile(thisdir,'MR_PL_Amb_ExoAtt_data',Observer);
    if ~isdir(datadirname);
        disp(sprintf('Making Observer directory %s',datadirname));
        mkdir(datadirname);
    end

disp(sprintf('[ DATA ] saving data in: %s',datadirname));
addpath(genpath('/Users/purplab/Desktop/Mariel/'));

stimulus = [];
stimulus.contrasts=Contrasts;
stimulus.Tilt =4; % degrees from horizontal - input value from Day tilt thresholding procedure

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
myscreen.keyboard.nums = [84 85]; % 1 & 2 on keyboard; 1 = CCW, 2 = CW

% initalize the screen
myscreen = initScreen(myscreen);
if stimulus.EyeTrack
    myscreen = eyeCalibDisp(myscreen);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

task{1}.waitForBacktick = 1;
task{1}.segmin =     [1 Inf .06 .06 .12 .7 5 1 30 5 .1]; % 1: ITI 2: Wait until fixation 3: Precue 4: ISI 5: Target + Response Cue 6: ISI + Response Cue 7: Response Window + Response Cue (fixation cross white and Observer can respond) 8: eyetracking 
task{1}.segmax =     [1 Inf .06 .06 .12 .7 5 1 30 5 .1];
task{1}.getResponse = [0 0 0 0 0 0 1 0 0 0 0]; % Observer can only input response in last segment 7 

n_repeats = 20; % n_repeats * stimulus conditions (4) = 80 trials with break halfway

if Diag == 1 % target can only be in UR or LL quadrant during training for that observer
    [Contrasts, location, repeat] = ndgrid(1:2,[1,3],1:n_repeats); % 2 contrasts: rmax & c50, two training locations, *n repetitions)
else Diag == 2 % target can only be in UL or LR quadrant during training for that observer
    [Contrasts, location, repeat] = ndgrid(1:2,[2,4],1:n_repeats); % 2 contrasts: rmax & c50, two training locations, *n repetitions)
end

task{1}.numTrials = length(location(:)); 
random_order = randperm(task{1}.numTrials);
task{1}.randVars.contrast = Contrasts(random_order);
task{1}.randVars.targetLocation = location(random_order);
task{1}.randVars.uniform.targetOrientation = 1:2;
task{1}.randVars.uniform.distractorOrientation1 = 1:2; % only one Gabor distractor on each trial
task{1}.randVars.uniform.distractorOrientation2 = 1:2;
task{1}.randVars.uniform.distractorOrientation3 = 1:2;
task{1}.randVars.len_ = task{1}.numTrials;
stimulus.trialend = 0;
stimulus.trialnum = 1;
stimulus.FixationBreak = zeros(1,length(location(:)));
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
    % runs automatically the task, you only need to change: StartSegmentCallback,DrawStimulusCallback,responseCallback
    [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end
clear stimulus.tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = endTask(myscreen,task);

filename_mat = sprintf('%s_D%d_B%d.mat',Observer,Day,Block);

if task{1, 1}.trialnum >= 80 % ensures incomplete blocks do not override full blocks bc same name
    save(fullfile(datadirname,filename_mat)); 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 1: function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = StartSegmentCallback(task, myscreen);
    % segments: 1:ITI,   2:fixation,    3:stimulus, 4:response
    global stimulus;
    global fb_count;

    
    if (task.thistrial.thisseg == 1) % ITI
        stimulus.trialend = stimulus.trialend + 1;
        
    elseif (task.thistrial.thisseg == 2) % fixation
        task.thistrial.fixationbreak = 0;%%%%% TO ADD FOR ONLINE EYETRACKING %%%
        stimulus.tmp.targetLocation  = stimulus.eccentricity*[stimulus.locations{task.thistrial.targetLocation}];
        stimulus.tmp.distractorIndices=stimulus.LocationIndices(stimulus.LocationIndices~=task.thistrial.targetLocation);
        for Locs=1:length(stimulus.tmp.distractorIndices)
            stimulus.tmp.distractorLocations{Locs}= stimulus.eccentricity*[stimulus.locations{stimulus.tmp.distractorIndices(Locs)}];
        end
        stimulus.FixationStarted=0;

        %response cue
        stimulus.tmp.respcueLocation=stimulus.respcueLocation{task.thistrial.targetLocation};
        stimulus.CueSounded=0;
        
    elseif (task.thistrial.thisseg == 8);%%%%% TO ADD FOR ONLINE EYETRACKING
        
        task.randVars.fixBreak(stimulus.trialnum) = task.thistrial.fixationbreak;
        
        fb_count = fb_count + 1; % add fb to counter
        
        if stimulus.trialnum >= 6
            
            recent_fb = sum(task.randVars.fixBreak(stimulus.trialnum-5:stimulus.trialnum-1)); % count how many fix breaks happened in the previous 5 trials 
        
            if fb_count >= 3 && recent_fb >= 3
                
               myscreen = eyeCalibDisp(myscreen); % return back to initial calibration display     

               fb_count = 0; % reset counter back to 0
               
            end 
            
        end
        
        % recreate params for fixation break trials and create and add new trial to end of block
        %% Contrast
        if task.randVars.fixBreak(stimulus.trialnum) == 1;
            temp = zeros(1,size(task.randVars.contrast,2)+1);
            temp(1:size(task.randVars.contrast,2)) = task.randVars.contrast;
            temp(end) = task.randVars.contrast(stimulus.trialnum);
            task.randVars.contrast = temp;
        end
        
        %% targetLocation
        if task.randVars.fixBreak(stimulus.trialnum) == 1;
            temp = zeros(1,size(task.randVars.targetLocation,2)+1);
            temp(1:size(task.randVars.targetLocation,2)) = task.randVars.targetLocation;
            temp(end) = task.randVars.targetLocation(stimulus.trialnum);
            task.randVars.targetLocation = temp;
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
        
    elseif (task.thistrial.thisseg == 7); % response %%%%
        if ~task.thistrial.gotResponse
            mglPlaySound(stimulus.CueSound);
        end;
        
    end


mglClearScreen(stimulus.grayColor);
setGammaTable(1);


%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % TASK 1: function that gets called to draw the stimulus each frame
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = DrawStimulusCallback(task, myscreen)
global stimulus;
    
mglClearScreen(stimulus.grayColor); %###
    
if (task.thistrial.thisseg == 1) % ITI
        mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.black); % note that I changed fix cross color from default white to black 
        
        % !!! decide whether to add in placeholders
%         for Gabor=stimulus.LocationIndices
%             for corner=2:5
%                 mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeHolderSize,stimulus.black)
%             end
%         end
     
elseif (task.thistrial.thisseg == 2) % FIXATION
        mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.black);
        
        % !!! decide whether to add in placeholders
%         for Gabor=stimulus.LocationIndices
%             for corner=2:5
%                 mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeHolderSize,stimulus.black)
%             end
%         end
        
        if stimulus.EyeTrack
                ep=myscreen.eyetracker.eyepos;
            if (sqrt(ep(end,1)^2+ep(end,2)^2))<=stimulus.TrialStartFixDist && ~stimulus.FixationStarted
                stimulus.FixationStart=mglGetSecs;
                stimulus.FixationStarted=1;
            elseif (sqrt(ep(end,1)^2+ep(end,2)^2))>=stimulus.TrialStartFixDist && stimulus.FixationStarted
                stimulus.FixationStarted=0;
            elseif (sqrt(ep(end,1)^2+ep(end,2)^2)) <= stimulus.TrialStartFixDist
                stimulus.FixationDur=mglGetSecs(stimulus.FixationStart);
                if stimulus.FixationDur >=stimulus.TrialStartFixDur
                    task = jumpSegment(task);
                end
            end
        else 
                mglWaitSecs(.25) % waits 250 ms before exo cue comes on
                task = jumpSegment(task);
        end
       
elseif (task.thistrial.thisseg == 3) % Exo cue
    mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.black);
    
    % !!! decide whether to add in placeholders
%     for Gabor=stimulus.LocationIndices
%         for corner=2:5
%             mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeHolderSize,stimulus.black)
%         end
%     end 
            
%     % precue for 100% valid exogenous training group !!!!! 

% from DPF exo att - exo cue of placeholders
% if task.thistrial.ExoCueCondition==2
%        for Gabor=stimulus.LocationIndices
%            for corner=2:5
%                mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.ExoCueSize,stimulus.white)
%            end
%        end
%     elseif task.thistrial.ExoCueCondition==1
       for Gabor=task.thistrial.targetLocation
           for corner=2:5
               mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.PrecueSize,stimulus.white)
           end
       end
%     end

%         mglGluDisk(stimulus.tmp.targetLocation*.58,stimulus.tmp.targetLocation*.58,stimulus.PrecueSize,stimulus.white);
%         mglGluDisk(stimulus.tmp.targetLocation*1.42,stimulus.tmp.targetLocation*1.42,stimulus.PrecueSize,stimulus.white);
%         mglGluDisk(stimulus.ExocueLocation(stimulus.LocationIndices(1)),stimulus.ExocueLocation(stimulus.LocationIndices(1)),stimulus.PrecueSize,stimulus.white);
%         mglGluDisk(stimulus.ExocueLocation(stimulus.LocationIndices(2)),stimulus.ExocueLocation(stimulus.LocationIndices(2)),stimulus.PrecueSize,stimulus.white); 
    
    if stimulus.EyeTrack;%%%%% TO ADD FOR ONLINE EYETRACKING
        ep=myscreen.eyetracker.eyepos;
        if (sqrt(ep(end,1)^2+ep(end,2)^2))>stimulus.TrialStartFixDist
           stimulus.FixationBreak(stimulus.trialnum)=1;
           task.thistrial.fixationbreak = 1;
           task = jumpSegment(task);
        end
    end
                
elseif (task.thistrial.thisseg == 4) % ISI
       mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.black);
       
       % !!! decide whether to add in placeholders
%        for Gabor=stimulus.LocationIndices
%            for corner=2:5
%                mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeHolderSize,stimulus.black)
%            end
%        end
                
       if stimulus.EyeTrack%%%%% TO ADD FOR ONLINE EYETRACKING
          ep=myscreen.eyetracker.eyepos;
          if task.thistrial.fixationbreak == 1
             task = jumpSegment(task);
          end;
          if (sqrt(ep(end,1)^2+ep(end,2)^2))>stimulus.TrialStartFixDist
             stimulus.FixationBreak(stimulus.trialnum)=1;
             task.thistrial.fixationbreak = 1;
             task = jumpSegment(task);
          end
       end
                
elseif (task.thistrial.thisseg == 5) % GABOR STIMULI & RESPONSE CUE
       mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.black);
       
       % DRAWS RESPONSE CUE
       mglLines2(stimulus.tmp.respcueLocation(1), stimulus.tmp.respcueLocation(2), stimulus.tmp.respcueLocation(3), stimulus.tmp.respcueLocation(4),2.5,stimulus.white);
       
       % !!! decide whether to add in placeholders
%        for Gabor=stimulus.LocationIndices
%            for corner=2:5
%                mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeHolderSize,stimulus.black)
%            end
%        end
                
       if stimulus.EyeTrack%%%%% TO ADD FOR ONLINE EYETRACKING
          ep=myscreen.eyetracker.eyepos;
          if task.thistrial.fixationbreak == 1
             task = jumpSegment(task);
          end;
          if (sqrt(ep(end,1)^2+ep(end,2)^2))>stimulus.TrialStartFixDist
           stimulus.FixationBreak(stimulus.trialnum)=1;
           task.thistrial.fixationbreak = 1;
           task = jumpSegment(task);
          end
       end
    
     % draws target gabor  
     drawGabor(stimulus.contrasts(task.thistrial.contrast),stimulus.tmp.targetLocation, stimulus.rotation(task.thistrial.targetOrientation), 1);
     % draws distractor gabor
     for Loc=1:length(stimulus.tmp.distractorIndices)
         eval(sprintf('drawGabor(stimulus.contrasts(task.thistrial.contrast),stimulus.tmp.distractorLocations{Loc}, stimulus.rotation(task.thistrial.distractorOrientation%g), 1);',Loc));
     end
                
elseif (task.thistrial.thisseg == 6) % ISI-2
       mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.black);
       
       % DRAWS RESPONSE CUE
       mglLines2(stimulus.tmp.respcueLocation(1), stimulus.tmp.respcueLocation(2), stimulus.tmp.respcueLocation(3), stimulus.tmp.respcueLocation(4),2.5,stimulus.white);
       
       % !!! decide whether to add in placeholders
%        for Gabor=stimulus.LocationIndices
%            for corner=2:5
%                mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeHolderSize,stimulus.black)
%            end
%        end

       if task.thistrial.fixationbreak == 1
          task = jumpSegment(task);
       end;
       
elseif (task.thistrial.thisseg == 7) % RESPONSE WINDOW
    
    %!!! check if color of fixation cross changes from black to white !!!!
       mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
       % !!! decide whether to add in placeholders
%        for Gabor=stimulus.LocationIndices
%            for corner=2:5
%                mglGluDisk(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeHolderSize,stimulus.black)
%            end
%        end
%        
       % DRAWS RESPONSE CUE
       mglLines2(stimulus.tmp.respcueLocation(1), stimulus.tmp.respcueLocation(2), stimulus.tmp.respcueLocation(3), stimulus.tmp.respcueLocation(4),2.5,stimulus.white);
       
       % RESPONSE TONE FEEDBACK
       if ~stimulus.CueSounded
          mglPlaySound(stimulus.CueSound);
          stimulus.CueSounded=1;
       end
                
       if task.thistrial.gotResponse
          task = jumpSegment(task);
       end;

       if task.thistrial.fixationbreak
          task = jumpSegment(task);
       end
                
elseif (task.thistrial.thisseg == 8) %%%%% TO ADD FOR ONLINE EYETRACKING
       if ~task.thistrial.fixationbreak
          task = jumpSegment(task);
       end
       if task.thistrial.fixationbreak
          mglTextSet('Courier',50,stimulus.black);
          mglTextDraw('Please fixate',[0 0]);
       end

elseif (task.thistrial.thisseg == 9) % halfway through block, text will appear onscreen to ask observer to take 1 min break

  if mod(task.numTrials,2) == 0
    
        halfblock = ((task.numTrials/2)+1);
        
        if stimulus.trialnum == halfblock
           mglTextSet('Courier',50,stimulus.black);
           mglTextDraw('Please take a short break and close your eyes!',[0,0]);
        else
            task = jumpSegment(task); % if it's not halfway through block, skip this segment
        end
        
    else mod(task.numTrials,2) == 1
        
        halfblock = ((task.numTrials+1)/2)+1;
        
        if stimulus.trialnum == halfblock %% NOTE: will only work if num trials is even number, so if with added fix break trials it becomes odd, the text will not show
           mglTextSet('Courier',50,stimulus.black);
           mglTextDraw('Please take a short break!',[0,0]);
        else
            task = jumpSegment(task); % if it's not halfway through block, skip this segment
        end        
  end

elseif (task.thistrial.thisseg == 10) % end of block feedback on screen for observer 

    if stimulus.trialnum<=task.numTrials % End of block Feedback
        task = jumpSegment(task);
    else
        CorrectCount=0;
        IncorrectCount=0;
        for Stars=1:length(stimulus.starColorFeedback)
            if stimulus.starColorFeedback(Stars)==2
                CorrectCount=CorrectCount+1;
            elseif stimulus.starColorFeedback(Stars)==3
                IncorrectCount=IncorrectCount+1;
            end
        end
        PercentCorrect=([num2str(round(100*(CorrectCount/(task.numTrials)))) num2str('%')]);
        NumFixBreak=num2str(sum(stimulus.FixationBreak));
      
        mglTextSet('Courier',50,stimulus.black);
        eval(sprintf('mglTextDraw(''%s correct'',[0 0]);',PercentCorrect));
        eval(sprintf('mglTextDraw(''%s fixation breaks'',[0 -1]);',NumFixBreak));
        
    end  

elseif (task.thistrial.thisseg == 11) % End of block Feedback presented to experimenter in command window
 
    if stimulus.trialnum<=task.numTrials % 
        task = jumpSegment(task);
    else   
        CorrectCount=0;
        IncorrectCount=0;
        for Stars=1:length(stimulus.starColorFeedback)
            if stimulus.starColorFeedback(Stars)==2
                CorrectCount=CorrectCount+1;
            elseif stimulus.starColorFeedback(Stars)==3
                IncorrectCount=IncorrectCount+1;
            end
        end
        PercentCorrect=([num2str(round(100*(CorrectCount/(task.numTrials)))) num2str('%')]);
        NumFixBreak=num2str(sum(stimulus.FixationBreak));
      
        disp(sprintf('Observer performance was %s correct',PercentCorrect));
        disp(sprintf('Observer broke fixation %s times',NumFixBreak));
    end 

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % function to get the Observer's response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = responseCallback(task, myscreen)
global stimulus;
mglClearScreen(stimulus.grayColor); %###
if ~task.thistrial.gotResponse

    % check response correct or not
    stimulus.tmp.response = [task.thistrial.whichButton == (task.thistrial.targetOrientation)]; %1 for CCW and 2 for CW

    % give feedback:
    if stimulus.tmp.response
       mglPlaySound(stimulus.CorrectSound);
       stimulus.starColorFeedback(stimulus.trialnum)=2;
       task = jumpSegment(task);
    else
       mglPlaySound(stimulus.IncorrectSound);
       stimulus.starColorFeedback(stimulus.trialnum)=3;
       task = jumpSegment(task);
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to draw the gabor stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawGabor(desiredContrast,position,orientation,sf);
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
% gabors
stimulus.width = 4*.8;%stimulus.gaussSdx*7 = 3.2;             % in deg
stimulus.height = 4*.8;%stimulus.gaussSdy*7;            % in deg
% stimulus.gaussSdx = 1; %0.8;  %0.3; %0.5; %1; % stimulus.width/7;                % in deg
% stimulus.gaussSdy = 1; %0.8;  %0.3; %0.5; %1; % stimulus.height/7;               % in deg
stimulus.gaussSdx = stimulus.width/7;                % in deg
stimulus.gaussSdy = stimulus.height/7;               % in deg

stimulus.rotation = [stimulus.Tilt -stimulus.Tilt]; % this is the tilt orientation of the gabor stimulus from vertical in Degrees
stimulus.DistractorRotation = 0; % this is the tilt orientation of the gabor stimulus from vertical in Degrees
% stimulus.DistractorRotation = [stimulus.Tilt -stimulus.Tilt]; % this is the tilt orientation of the gabor stimulus from vertical in Degrees
stimulus.init = 1;

stimulus.sf = 6;   % in cpd
stimulus.orientation = 0;      % in deg
% original: stimulus.orientation = 90;      % in deg
stimulus.phase = 0;             % in deg

stimulus.eccentricity = 4;  % in deg !!! may change
% original stimulus.eccentricity = 8*.8;  % 6.4 deg


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
    for i=1:length(stimulus.locations)
        stimulus.placeholders{i}= [stimulus.eccentricity*stimulus.locations{i}];
        stimulus.placeholders{i}(2,:)=[stimulus.placeholders{i}(1,1:2)]+[0 stimulus.cornerDist];
        stimulus.placeholders{i}(3,:)=[stimulus.placeholders{i}(1,1:2)]+[stimulus.cornerDist 0];
        stimulus.placeholders{i}(4,:)=[stimulus.placeholders{i}(1,1:2)]+[0 -stimulus.cornerDist];
        stimulus.placeholders{i}(5,:)=[stimulus.placeholders{i}(1,1:2)]+[-stimulus.cornerDist 0];
    end
    
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
        stimulus.ExocueLocation{i}(1:2)=(stimulus.eccentricity+stimulus.ExoCueDist)*stimulus.locations{stimulus.LocationIndices(i)};
        stimulus.ExocueLocation{i}(3:4)=(stimulus.eccentricity+stimulus.ExoCueDist)*stimulus.locations{stimulus.LocationIndices(i)};
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
