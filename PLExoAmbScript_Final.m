%% PL Exo Amb Script %%

% Main experimental script or perceptual learning experiment with exogenous
% attention cues in amblyopes

% Last updated 9.18.21

% First created by Mariel Roberts on 8.16.17

% Important: run this cell FIRST each day to set screen params, add paths, and set directory %%

    % add main folder to path
    addpath(genpath('/Users/purplab/Desktop/Mariel/'));

    % redirect to be inside PLExoAmb experiment folder
    cd /users/purplab/Desktop/Mariel/PLExoAmb

    % check using correct screen params
    mglEditScreenParams;
    % check using correct eye calibration parameters ( HV5 and calib area
    % of .48
    mglEyelinkParams;       
    
    
    %% Days 2 OR 15: qCSF pre- OR post-test measurement %%

    % Obtain pre-test OR post-test quick contrast response function 
    
    % 4AFC discrimination task (V,+45,H,-45) at fovea using 1-4 keys on top
    
    % Note that Gabor contrast and SF will change on every trial; by
    % design, on some trials you may not see the target at all! 

    %% Practice %%         
    
    % First, we want them to learn the procedural demands of the task with 
    % one block (70 total trials) of easier (i.e. higher contrast) practice trials

% IMPORTANT: Observers should perform the practice BINOCULARLY

% If post-test, run a few practice trials to get them reacquainted with
% task, but do not need to run full block
        
        clear;
        close all;
        clc;
        
        Observer = input('Observer initials (eg, XX):  ','s'); % subject initials, in single quotes
        Eyetrack = input('Eyetracking on? 1 = Y, 0 = N:  '); % 0 = no eyetracking, 1 = online eyetracking
        
        PLExoAmb_PRACTICEqCSF_4AFC_FOVEA_PLACEHOLDERS(Observer,Eyetrack); 
            
% IMPORTANT: If observer performance is at chance (<65%) because they are finding it hard to
% map stimulus to motor response, run the practice block again
    
     %% Main qCSF pre- or post-test measurement %%
    
        % Obtain the real qCSF estimates (70 trials/block)
        % 3x each for each eye, interleaved in separate blocks
    
% IMPORTANT: MAKE SURE THE OBSERVER COMPLETES THE 6 BLOCKS UNDER THE
% CORRECT MONOCULAR VIEWING CONDITIONS IN THIS SPECIFIC ORDER
          
        Observer = input('Observer initials (eg, XX):  ','s'); 
        Eyetrack = input('Eyetracking on? 1 = Y, 0 = N:  '); % 0 = no eyetracking, 1 = online eyetracking
        Session = input('Pre-test (0) or Post-test (1)?:  '); 

        FE_qCSF_AULCSF_vals = NaN(1,3);
        AE_qCSF_AULCSF_vals = NaN(1,3);
        
        % estimates of model params:
        % [(1) peak gain (2) peak SF (3) bandwidth (width at half height) (4) truncation
        
        FE_qCSF_modelparams_vals = NaN(3,4);
        AE_qCSF_modelparams_vals = NaN(3,4);   
        
        % estimate of full CSF function
        FE_qCSF_vals = NaN(3,12);
        AE_qCSF_vals = NaN(3,12);   
        

        %% FELLOW (non-amblyopic eye) with amblyopic eye patched %%

        Eye = input('Eye tested? F = Fellow, A = Amblyopic:  ','s');   
        Block = input('Block number? 1-3:  '); 
        PLExoAmb_qCSF_4AFC_FOVEA_PLACEHOLDERS(Observer,Eyetrack,Session,Eye,Block);  
        
        FE_qCSF_AULCSF_vals(1,1) = qcsf.data.estAULCSF(end);
        FE_qCSF_modelparams_vals(1,:) = qcsf.data.estCSF(end,:);
        FE_qCSF_vals(1,:) = qcsf.data.estSensitivity(end,:);

        %% AMBLYOPIC with fellow eye patched %%
        
        Eye = input('Eye tested? F = Fellow, A = Amblyopic:  ','s');   
        Block = input('Block number? 1-3:  '); 
        PLExoAmb_qCSF_4AFC_FOVEA_PLACEHOLDERS(Observer,Eyetrack,Session,Eye,Block);
        
        AE_qCSF_AULCSF_vals(1,1) = qcsf.data.estAULCSF(end);
        AE_qCSF_modelparams_vals(1,:) = qcsf.data.estCSF(end,:);
        AE_qCSF_vals(1,:) =  qcsf.data.estSensitivity(end,:);

        %% FELLOW (non-amblyopic eye) with amblyopic eye patched %%
        
        Eye = input('Eye tested? F = Fellow, A = Amblyopic:  ','s');   
        Block = input('Block number? 1-3:  '); 
        PLExoAmb_qCSF_4AFC_FOVEA_PLACEHOLDERS(Observer,Eyetrack,Session,Eye,Block);  
        
        FE_qCSF_AULCSF_vals(1,2) = qcsf.data.estAULCSF(end);
        FE_qCSF_modelparams_vals(2,:) = qcsf.data.estCSF(end,:);
        FE_qCSF_vals(2,:) = qcsf.data.estSensitivity(end,:);
        
%% IMPORTANT: Ask observer to take a ~5 min break ! Ask if they want water or to use the bathroom. %%
        
        %% AMBLYOPIC with fellow eye patched %%
        
        Eye = input('Eye tested? F = Fellow, A = Amblyopic:  ','s');   
        Block = input('Block number? 1-3:  '); 
        PLExoAmb_qCSF_4AFC_FOVEA_PLACEHOLDERS(Observer,Eyetrack,Session,Eye,Block);  
        
        AE_qCSF_AULCSF_vals(1,2) = qcsf.data.estAULCSF(end);
        AE_qCSF_modelparams_vals(2,:) = qcsf.data.estCSF(end,:);
        AE_qCSF_vals(2,:) = qcsf.data.estSensitivity(end,:);

        %% FELLOW (non-amblyopic eye) with amblyopic eye patched %%

        Eye = input('Eye tested? F = Fellow, A = Amblyopic:  ','s');   
        Block = input('Block number? 1-3:  '); 
        PLExoAmb_qCSF_4AFC_FOVEA_PLACEHOLDERS(Observer,Eyetrack,Session,Eye,Block); 
        
        FE_qCSF_AULCSF_vals(1,3) = qcsf.data.estAULCSF(end);
        FE_qCSF_modelparams_vals(3,:) = qcsf.data.estCSF(end,:);
        FE_qCSF_vals(3,:) = qcsf.data.estSensitivity(end,:);
        
        %% AMBLYOPIC with fellow eye patched %%
        
        Eye = input('Eye tested? F = Fellow, A = Amblyopic:  ','s');   
        Block = input('Block number? 1-3:  '); 
        PLExoAmb_qCSF_4AFC_FOVEA_PLACEHOLDERS(Observer,Eyetrack,Session,Eye,Block);  
        
        AE_qCSF_AULCSF_vals(1,3) = qcsf.data.estAULCSF(end);
        AE_qCSF_modelparams_vals(3,:) = qcsf.data.estCSF(end,:);
        AE_qCSF_vals(3,:) = qcsf.data.estSensitivity(end,:);

    %% calculate mean ± SD est AULCSF for both eyes
        M_FE_qCSF_AULCSF = geomean(FE_qCSF_AULCSF_vals);
        SD_FE_qCSF_AULCSF = nanstd(FE_qCSF_AULCSF_vals);
        
        M_AE_qCSF_AULCSF = geomean(AE_qCSF_AULCSF_vals(2:3));
        SD_AE_qCSF_AULCSF = nanstd(AE_qCSF_AULCSF_vals(2:3));
 

        M_FE_qCSF_modelparams = nanmean(FE_qCSF_modelparams_vals);
        SD_FE_qCSF_modelparams = nanstd(FE_qCSF_modelparams_vals);
        
        M_AE_qCSF_modelparams = nanmean(AE_qCSF_modelparams_vals);
        SD_AE_qCSF_modelparams = nanstd(AE_qCSF_modelparams_vals);
        
        
        M_FE_qCSF_values = nanmean(FE_qCSF_vals);
        SD_FE_qCSF_values = nanstd(FE_qCSF_vals);
        
        M_AE_qCSF_values = nanmean(AE_qCSF_vals);
        SD_AE_qCSF_values = nanstd(AE_qCSF_vals);
        
  
        % save all calculated values in mat file in observer folder
        
        if Session == 0
           Session = 'pre'
        else Session == 1
           Session = 'post'
        end
        
        qcsf_filename = sprintf('%s_qCSF_%s',Observer,Session);
        datadirname = fullfile(pwd,'MR_PL_Amb_ExoAtt_data',Observer);
        save(fullfile(datadirname, qcsf_filename));
        
        % display all values in command line for ease of uploading to wiki
        
%% Day 3: Contrast thresholding procedure (Best PEST) %%%

% Run two simultaneous staircases to measure threshold estimates at 75% and 88% accuracy 

    % Set observer and clear previous staircases %

        clear;
        clc;
        
    % subject initials, in single quotes
        Observer = input('Observer initials (eg, XX):  ','s'); 
    % 0 = no eyetracking, 1 = online eyetracking
        Eyetrack = input('Eyetracking on? 1 = Y, 0 = N:  '); 

    %% Practice block: 40 trials 2AFC discrim of foveal target, 64% contrast, 2 cpd Gabor 
    
        PLExoAmb_PracticeFovea_Final(Observer,Eyetrack);
       
     %% Practice block: 80 trials 2AFC discrim, 64% contrast, 6 cpd Gabor, 4 deg tilt @ 4 ecc along trained and untrained diagonal once
    
% IMPORTANT: The observer should perform the practice with BOTH EYES
% along each diagonal once
    
        PLExoAmb_Practice_Final(Observer,Eyetrack,1);
        PLExoAmb_Practice_Final(Observer,Eyetrack,2);
         
% IMPORTANT: If observer performance is at chance (<65%), run the
% practice block repeatedly until they are at >85% correct
   
%% Main contrast threshold procedure     
          
    % run Best PEST contrast thresholding procedure to estimate observer's contrast thresholds at 75% and 88% accuracy simultaneously 
    % Run four full blocks per eye (first block for each eye excluded from estimate since uniform prior) of 40 trials each (~30-40 mins total) 
    
% IMPORTANT: The observer should FIRST perform this task with their FELLOW EYE
% along the diagonal at which they'll be trained (randomly pre-assigned in wiki)
    
    % observers train on 1 = UR/LL quadrants or 2 = UL/LR quadrants
    Diag = input('Trained diagonal? 1 = UR/LL, 2 = UL/LR quadrants:  '); 
    Eye = input('The fellow eye is the left(L) or right(R) eye?:  ','s');

 % initialize the 1st staircase parameters - estimate threshold at 88% %

        stairParams1.whichStair = 1; % 1 = best PEST; 2 = QUEST
        stairParams1.alphaRange = .01:.01:1; % set range of possible contrast thresholds
        stairParams1.fitBeta = 2; % slope parameter
        stairParams1.fitLambda = 0.01; % lapse rate parameter (distance from asymptote)
        stairParams1.fitGamma = 0.5; % floor parameter, chance performance is 50%
        stairParams1.PF = 'arbWeibull'; % PF will be the Weibull function but changed to have an arbritrary performance level
        
        %%% expected staircase convergence is 88% performance - typical for 2AFC %%%
        stairParams1.threshPerformance = 0.75; 
        
 % initialize the 2nd staircase parameters - estimate threshold at 75% %

        stairParams2.whichStair = 1; % 1 = best PEST; 2 = QUEST
        stairParams2.alphaRange = .01:.01:1; % set range of possible contrast thresholds
        stairParams2.fitBeta = 2; % slope parameter
        stairParams2.fitLambda = 0.01; % lapse rate parameter (distance from asymptote)
        stairParams2.fitGamma = 0.5; % floor parameter, chance performance is 50%
        stairParams2.PF = 'arbWeibull'; % PF will be the Weibull function but changed to have an arbritrary performance level
        
        %%% expected staircase convergence is 75% performance - typical for 2AFC %%%
        stairParams2.threshPerformance = 0.88;         

 % set uniform priors across tilt threshold range for initial staircases %

        uniform_prior1 = ones(size(stairParams1.alphaRange));
        uniform_prior2 = ones(size(stairParams2.alphaRange));
       
    % create empty matrix to store each contrast threshold estimate
    % row 1 is estimate for 88% and row 2 is estimate is for 75%
    % columns are different blocks (averaging just the final 3)
        contrastthresh_est = NaN(2,4);
        
    % run function 4 times, and store that block's ct estimate
    
% IMPORTANT: if this full loop is not completed, you must manually check which
% block is missing and change function input accordingly
    
    for Block = 1:4
        
        % if first block, then use uniform prior
        if Block == 1
            
         stairParams1.lastPosterior = uniform_prior1;
         stairParams2.lastPosterior = uniform_prior2;
         
        % if already run one block, use the last posterior as the new prior 
        else 
            
         stairParams1.lastPosterior = stair1.pdf;
         stairParams2.lastPosterior = stair2.pdf;
         
        end
        
     % initialize the staircases %
    [stair1] = usePalamedesStaircase(stairParams1);
    [stair2] = usePalamedesStaircase(stairParams2);
        
    % run contrast thresholding program for block
    [stair1,stair2,ct1_ToPlot,ct2_ToPlot] = PLExoAmb_ContrastThresh_Final(Observer,Eyetrack,Diag,Block,Eye,stair1,stair2);
        
    % place ct estimates for that block in matrix
    contrastthresh_est(1,Block) = stair1.xCurrent;
    contrastthresh_est(2,Block) = stair2.xCurrent;
        
    end
        
           
    if ~isnan(contrastthresh_est) % if none of the values are NaN

       if ~any(contrastthresh_est == 0) % if none of the values are 0

           if size(contrastthresh_est) == [2 4] % if the matrix is size 2x3

               % calculate and print observer contrast threshold estimate %
               % take geometric mean (since measures are not independent and less biased to outliers) of last three separate estimates
                Contrast_thresh_est = geomean(contrastthresh_est(:,2:4),2);
                std_Contrastthresh_est1 = std(contrastthresh_est(1,2:4));
                std_Contrastthresh_est2 = std(contrastthresh_est(2,2:4));

                % save the contrast threshold for this observer in their folder %
                msg = sprintf('The contrast thresholds for 75 and 85 accuracy in the fellow eye of observer %s are %g ± %g and %g ± %g.', Observer, Contrast_thresh_est(1), std_Contrastthresh_est1,Contrast_thresh_est(2), std_Contrastthresh_est2);
                disp(msg); 
                cte_file = sprintf('%s_%s_Contrast_thresh_est',Observer,Eye);        
                datadirname = fullfile(pwd,'MR_PL_Amb_ExoAtt_data',Observer);
                save(fullfile(datadirname,cte_file)); % contrast threshold estimates are saved in Observer file

            else  
                msg = sprintf('The observer has not completed all blocks');
                disp(msg);
            end

        else  
            msg = sprintf('The observer has not completed all blocks');
            disp(msg);
       end

    else  
    msg = sprintf('The observer has not completed all blocks');
    disp(msg);

    end

%% IMPORTANT: The observer should next perform this task with their AMBLYOPIC EYE
% along the diagonal at which they'll be trained (randomly pre-assigned in wiki)

        Eye = input('The amblyopic eye is the left(L) or right(R) eye?:  ','s');

     % initialize the 1st staircase parameters - estimate threshold at 88% %

            stairParams1.whichStair = 1; % 1 = best PEST; 2 = QUEST
            stairParams1.alphaRange = .01:.01:1; % set range of possible contrast thresholds
            stairParams1.fitBeta = 2; % slope parameter
            stairParams1.fitLambda = 0.01; % lapse rate parameter (distance from asymptote)
            stairParams1.fitGamma = 0.5; % floor parameter, chance performance is 50%
            stairParams1.PF = 'arbWeibull'; % PF will be the Weibull function but changed to have an arbritrary performance level

            %%% expected staircase convergence is 88% performance - typical for 2AFC %%%
            stairParams1.threshPerformance = 0.75; 

     % initialize the 2nd staircase parameters - estimate threshold at 75% %

            stairParams2.whichStair = 1; % 1 = best PEST; 2 = QUEST
            stairParams2.alphaRange = .01:.01:1; % set range of possible contrast thresholds
            stairParams2.fitBeta = 2; % slope parameter
            stairParams2.fitLambda = 0.01; % lapse rate parameter (distance from asymptote)
            stairParams2.fitGamma = 0.5; % floor parameter, chance performance is 50%
            stairParams2.PF = 'arbWeibull'; % PF will be the Weibull function but changed to have an arbritrary performance level

            %%% expected staircase convergence is 75% performance - typical for 2AFC %%%
            stairParams2.threshPerformance = 0.88;         

     % set uniform priors across tilt threshold range for initial staircases %

            uniform_prior1 = ones(size(stairParams1.alphaRange));
            uniform_prior2 = ones(size(stairParams2.alphaRange));

        % create empty matrix to store each contrast threshold estimate
        % row 1 is estimate for 88% and row 2 is estimate is for 75%
        % columns are different blocks (averaging just the final 3)
            contrastthresh_est = NaN(2,4);

        % run function 4 times, and store that block's ct estimate

    % IMPORTANT: if this full loop is not completed, you must manually check which
    % block is missing and change function input accordingly

        for Block = 1:4
            
            % if first block, then use uniform prior
            if Block == 1

             stairParams1.lastPosterior = uniform_prior1;
             stairParams2.lastPosterior = uniform_prior2;

            % if already run one block, use the last posterior as the new prior 
            else 

             stairParams1.lastPosterior = stair1.pdf;
             stairParams2.lastPosterior = stair2.pdf;

            end

         % initialize the staircases %
        [stair1] = usePalamedesStaircase(stairParams1);
        [stair2] = usePalamedesStaircase(stairParams2);

        % run contrast thresholding program for block
        [stair1,stair2,ct1_ToPlot,ct2_ToPlot] = PLExoAmb_ContrastThresh_Final(Observer,Eyetrack,Diag,Block,Eye,stair1,stair2);

        % place ct estimates for that block in matrix
        contrastthresh_est(1,Block) = stair1.xCurrent;
        contrastthresh_est(2,Block) = stair2.xCurrent;
          
        end

      if ~isnan(contrastthresh_est) % if none of the values are NaN

       if ~any(contrastthresh_est == 0) % if none of the values are 0 % NOTE: if you skip a block, the est will input as .5, so still need to check all values

           if size(contrastthresh_est) == [2 4] % if the matrix is size 2x3

               % calculate and print observer contrast threshold estimate %
               % take geometric mean (since measures are not independent and less biased to outliers) of last three separate estimates
                Contrast_thresh_est = geomean(contrastthresh_est(:,2:4),2);
                std_Contrastthresh_est1 = std(contrastthresh_est(1,2:4));
                std_Contrastthresh_est2 = std(contrastthresh_est(2,2:4));

                % save the contrast threshold for this observer in their folder %
                msg = sprintf('The contrast thresholds for 75 and 85 accuracy in the amblyopic eye of observer %s are %g ± %g and %g ± %g.', Observer, Contrast_thresh_est(1), std_Contrastthresh_est1,Contrast_thresh_est(2), std_Contrastthresh_est2);
                disp(msg); 
                cte_file = sprintf('%s_%s_Contrast_thresh_est',Observer,Eye);        
                datadirname = fullfile(pwd,'MR_PL_Amb_ExoAtt_data',Observer);
                save(fullfile(datadirname,cte_file)); % contrast threshold estimates are saved in Observer file

            else  
                msg = sprintf('The observer has not completed all blocks');
                disp(msg);
            end

        else  
            msg = sprintf('The observer has not completed all blocks');
            disp(msg);
       end

    else  
    msg = sprintf('The observer has not completed all blocks');
    disp(msg);

    end
        
%% Days 3 OR 14: Main task pre- or post-test measurement

% 2AFC task (+4 or -4 from H), 6 cpd Gabor at 4 dva ecc

% Measure d' for 2 contrast levels (individual's d' contrast level at ~75% and 88% correct)
% x 2 diagonals (trained & untrained) x 2 eyes (1 using fellow and amblyopic eye - different ct per eye) 

% 160 trials/block * 4 blocks = 640 trials in session (~45 mins-1 hr) = 80
% trials per condition w/ breaks every 40 trials, switching eyes every main
% block of 160 trials
   
        clear;
        close all;
        clc;
        
        Observer = input('Observer initials (eg, XX):  ','s'); 
        Eyetrack = input('Eyetracking on? 1 = Y, 0 = N:  ');
        Session = input('Pre-test (0) or Post-test (1)?:  ');   

%% first test using FELLOW eye (with AMBLYOPIC eye patched) at diagonal 1

        Test_Eye = input('Testing which eye? 1 = FE, 2 = AE:  '); 
        Eye = input('Fellow eye is left(L) or right(R) eye?:  ','s');
        
        % automatically load observer's contrast threshold estimates
          load(sprintf('%s_%s_Contrast_thresh_est.mat',Observer,Eye));
          
          % clear unecessarily saved variables from ct 
          clear cte_file;
          clear datadirname;
          clear msg; 
          
        % !! change to reflect two different contrast estimates for 75% and
        % 88%
           
        Contrasts = [Contrast_thresh_est(1) Contrast_thresh_est(2)]; 
        
       % !! to add: spit error if ct threshold estimate is saved as 0 or NaN
        % test at diag 1 with AE
        
        PLExoAmb_Test_Final(Observer,Session,Eyetrack,Contrasts,1,Test_Eye);
        
%% now ask them to move patch to perform task with AMBLYOPIC eye at diagonal 1
  
        Test_Eye = input('Testing which eye? 1 = FE, 2 = AE:  ');
        Eye = input('Amblyopic eye is left(L) or right(R) eye?:  ','s');

          % need to automatically load contrast threshold for amblyopic eye
          load(sprintf('%s_%s_Contrast_thresh_est.mat',Observer,Eye));
          
          % clear unecessarily saved variables from ct 
          clear cte_file;
          clear datadirname;
          clear msg; 
          
        % !! change to reflect two different contrast estimates for 75% and
        % 88%
           
        Contrasts = [Contrast_thresh_est(1) Contrast_thresh_est(2)]; 
        
       % !! to add: spit error if ct threshold estimate is saved as 0 or NaN

        % test at diag 1 with AE
        PLExoAmb_Test_Final(Observer,Session,Eyetrack,Contrasts,1,Test_Eye);
        
%% now move patch back to perform task with FELLOW eye at diagonal 2

        Test_Eye = input('Testing which eye? 1 = FE, 2 = AE:  ');
        Eye = input('Fellow eye is left(L) or right(R) eye?:  ','s');

          % need to automatically load contrast threshold for amblyopic eye
          load(sprintf('%s_%s_Contrast_thresh_est.mat',Observer,Eye));
          
          % clear unecessarily saved variables from ct 
          clear cte_file;
          clear datadirname;
          clear msg; 
           
        % !! change to reflect two different contrast estimates for 75% and
        % 88%
           
        Contrasts = [Contrast_thresh_est(1) Contrast_thresh_est(2)]; 
        
        % !! to add: spit error if ct threshold estimate is saved as 0 or NaN

        PLExoAmb_Test_Final(Observer,Session,Eyetrack,Contrasts,2,Test_Eye); 
        
%% now move patch back to perform task with AMBLYOPIC eye at diagonal 2

        Test_Eye = input('Testing which eye? 1 = FE, 2 = AE:  ');
        Eye = input('Amblyopic eye is left(L) or right(R) eye?:  ','s');

          % need to automatically load contrast threshold for amblyopic eye
          load(sprintf('%s_%s_Contrast_thresh_est.mat',Observer,Eye));
          
          % clear unecessarily saved variables from ct 
          clear cte_file;
          clear datadirname;
          clear msg; 
           
        % !! change to reflect two different contrast estimates for 75% and
        % 88%
           
        Contrasts = [Contrast_thresh_est(1) Contrast_thresh_est(2)]; 
        
        % !! to add: spit error if ct threshold estimate is saved as 0 or NaN

        PLExoAmb_Test_Final(Observer,Session,Eyetrack,Contrasts,2,Test_Eye); 
        
%% Days 4-13 %%

    % Training Task % 

% IMPORTANT: Before they leave, check that they have completed at least 6
% blocks! 80 trials/block * 6 blocks = 480 trials/session * 10 sessions =
% 4800 trials

    Observer = input('Observer initials (eg, XX):  ','s');
    Eyetrack = input('Eyetracking on? 1 = Y, 0 = N:  '); % 0 = no eyetracking, 1 = online eyetracking
    Eye = input('Amblyopic eye is left(L) or right(R) eye?:  ','s');

          % need to automatically load contrast threshold for amblyopic eye
          load(sprintf('%s_%s_Contrast_thresh_est.mat',Observer,Eye));
          
          % clear unecessarily saved variables from ct 
          clear cte_file;
          clear datadirname;
          clear msg; 
           
          ct = Contrast_thresh_est; 
        
        % based on prior literature, all observers should be performing at asymptote at 64% contrast 
        dmax = .64; 
        Contrasts = [ct dmax];
    
    Diag = input('Trained diagonal? 1 = UR/LL, 2 = UL/LR quadrants:  ');    
    Day = input('Training Day (1-10)?:  ');                            

    
% IMPORTANT: for each session, each observer should complete 6 training blocks to keep amount of training consistent                                   
% IF, for some unplanned emergency or you run out of time in the hour, then
% make sure to write down EXACTLY how many trials completed per day in wiki

    n_blocks = 2; % change this to number of blocks you want to set to run automatically per half session

    % made two different functions according to type of precue
    AttCondition = input('With valid exo cues? 1 = Y, 0 = N:  '); 
    Startingblock = input('Start with which number block?  ');
    
    for Block = Startingblock:(Startingblock+n_blocks-1)
        if AttCondition == 0
            PLExoAmb_TrainNeutral(Observer,Eyetrack,Contrasts,Diag,Day,Block);
        else AttCondition == 1
            PLExoAmb_TrainExoValid(Observer,Eyetrack,Contrasts,Diag,Day,Block);
        end
    end
    
%% !!! goes back to home screen; remind them to take a longer break and come get the experimenter !!! 
        
    Startingblock = input('Start with which number block?  ');

 for Block = Startingblock:(Startingblock+n_blocks-1)
        if AttCondition == 0
            PLExoAmb_TrainNeutral(Observer,Eyetrack,Contrasts,Diag,Day,Block);
        else AttCondition == 1
            PLExoAmb_TrainExoValid(Observer,Eyetrack,Contrasts,Diag,Day,Block);
        end
 end
