%% MEGruns: 'rst' = preplay, 'lci' = lexicallyCuedImagination, 'str' = Structure Learning, 'rst' = resting state
%%          'rwd' = reward learning,  'seq' = sequence code, 'pos' = position code

%% list of artifacts that appeared in PCA:
%  1. most runs from most subjects have sensor 252 (MRT42) alone as the first principle component. it has some kind of massive,
%     low-frequency drift, independent of all the other sensors. but superimposed on this, it seems to contain good signal.
%  2. subjects 6,15,26,30 have the artifact i first noticed, a clump of right temporal/parietal/occipital sensors
%     anticorrelated with each other. 26 has it worst.
%  3. subjects 16,17,18,19,22,24,25,29 have a milder artifact on sensor 51 (MLF61).
%  4. subject 21 has a MASSIVE artifact on the posterior occipital sensors.
%  5. subjects 5, 23, and 28 all have different kinds of artifacts.
%  6. subject 16 on one trial has a huge artifact.

is.nSubj = 1;

%% [1] 
%% Note: Head position - normal, a bit towards frontal,the position might be too up, blockade some eyesight
%% Note: Task order - 1
%% Note: Training - Take 2 Rounds to finish (10 sessions)
%% Note: Age: 34
%% Note: tired,eye moment might be noisy;  
%% Note: tired, especially after the first session, 
%% Note: loss battry in the third session, also a lot of movement caused by the noise in the last run

is.fnDate{is.nSubj} = '2018-02-12'; is.fnSID{is.nSubj} = '01';
is.fnBehav{is.nSubj} = {'0001'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05206_Yunzhe_20180212'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'};  
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %% [2] Paradkar, Swarali - NO TRIGGER?
% %% Note: Head position - normal, a bit towards frontal
% %% Note: Task order - 1
% %% Note: Training - Take 1 Rounds to finish (5 sessions)
% %% Note: Age: 22
% %% Note: INTERFERENCE FROM THE EEG!  
% %% Note: THE PROJECTOR is not in the CENTRAL, lost left side VISION.   
% %% Note: didn't drink coffee, only mild tea, 
% %% Note: very bad in self maintaience 
% 
% is.fnDate{is.nSubj} = '2018-02-17'; is.fnSID{is.nSubj} = '02';
% is.fnBehav{is.nSubj} = {'0002'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
% is.fnMEG{is.nSubj} = 'MG05210_Yunzhe_20180217'; 
% is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'};  % the lci didn't start MEG recording until a few trials in.
% is.taskversion{is.nSubj}=2; 
% is.nSubj = is.nSubj + 1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [3]  MEDAL IN HER BRA !!!: Signal may be noisy
%% Note: Head position - normal
%% Note: Task order - 1
%% Note: Training - Take 1 Rounds to finish (5 sessions)
%% Note: Age: 22
%% Note: INTERFERENCE FROM THE EEG!  
%% Note: THE PROJECTOR is not in the CENTRAL, lost left side VISION.   
%% Note: short eye sight,  cannot see very clear; EYE GAZE is noisy

is.fnDate{is.nSubj} = '2018-02-18'; is.fnSID{is.nSubj} = '03';
is.fnBehav{is.nSubj} = {'0003'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05211_Yunzhe_20180218'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'};  
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [4]  go to pee, in the resting state after third learning, 
%% Note: Head position - normal
%% Note: Task order - 1
%% Note: Training - Take 1 Rounds to finish (5 sessions)
%% Note: Age: 22
%% Note: lost eye data for the third Associative Learning and onward
%% Note: 2nd resting is noisy, >30mm head movement 

is.fnDate{is.nSubj} = '2018-02-19'; is.fnSID{is.nSubj} = '04';
is.fnBehav{is.nSubj} = {'0004'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05212_Yunzhe_20180219'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'};  
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [5]  VERY BAD!, go to pee after the 2nd functional localizer, data might be corrupted
%% Note: Head position - normal,a bit towards back
%% Note: Task order - 1
%% Note: Training - Take 2  Rounds to finish (10  sessions)
%% Note: Age: 38!!!
%% Note: tons of movement, 
%% Note: fitting error up to > 80% in the third functional localizer. 

is.fnDate{is.nSubj} = '2018-02-19'; is.fnSID{is.nSubj} = '05';
is.fnBehav{is.nSubj} = {'0005'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05213_Yunzhe_20180219'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'};  
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [6] Cannot figure out the sequence - at least after 2 learing block, claim figure out the seq after the third run
%% Note: Head position - normal
%% Note: Task order - 2
%% Note: Training - Take 2  Rounds to finish (10  sessions)
%% Note: Age: 24, left handed
%% Note: pretty big, so we changed the chair position; 
%% Note: some muscle tension on his back and neck; 
%% Note: huge fitting error in the end of 2nd functional localizer. 
%% Note: Eye tracking data is bit noisy

is.fnDate{is.nSubj} = '2018-02-20'; is.fnSID{is.nSubj} = '06';
is.fnBehav{is.nSubj} = {'0006'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05214_Yunzhe_20180220'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [7] didn't get the sequence after the first learning block, claim to get it after the second
%% Note: Head position - normal, a bit of frontal
%% Note: Task order - 2
%% Note: Training - Take 1 Rounds to finish (5  sessions)
%% Note: Age: 19
%% Note: drink multiple coffees
%% Note: mixed 1 with 2 - NEED TO DO REVERSE CODING ON HIS BEHAVIOR RESPONSE!!!

is.fnDate{is.nSubj} = '2018-02-20'; is.fnSID{is.nSubj} = '07';
is.fnBehav{is.nSubj} = {'0007'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05215_Yunzhe_20180220'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [8] lab member
%% Note: Head position - normal
%% Note: Task order - 2
%% Note: Training - Take 1 Rounds to finish (5  sessions)
%% Note: Age: 30, left handed
%% Note: tired doing functional localizer; didn't drink coffee

is.fnDate{is.nSubj} = '2018-02-21'; is.fnSID{is.nSubj} = '08';
is.fnBehav{is.nSubj} = {'0008'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05217_Yunzhe_20180221'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [9] 
%% Note: Head position - normal, a bit of frontal, towards left
%% Note: Task order - 2
%% Note: Training - Take 1 Rounds to finish (5  sessions)
%% Note: Age: 23
%% Note: didn't drink coffee. Wear contact lenses, Eye tracking data is bit noisy
%% Note: head position is off for the third functional localizer, >10mm

is.fnDate{is.nSubj} = '2018-02-22'; is.fnSID{is.nSubj} = '09';
is.fnBehav{is.nSubj} = {'0009'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05219_Yunzhe_20180222'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'};  
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [10]
%% Note: Head position - normal, a bit of frontal
%% Note: Task order - 2
%% Note: Training - Take 1 Rounds, but did not finish, don't have practice data!
%% Note: Age: 34
%% Note: didn't dirnk coffee, cough a lot
%% Note: did not look at the pictures after first run, but still try to find sequence

is.fnDate{is.nSubj} = '2018-02-23'; is.fnSID{is.nSubj} = '10';
is.fnBehav{is.nSubj} = {'0010'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05220_Yunzhe_20180223'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'};  
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [11] 
%% Note: Head position - normal. A bit towards back
%% Note: Task order - 2
%% Note: Training - Take 2  Rounds to finish (10  sessions)
%% Note: Age: 24
%% Note: didn't drink coffee, head position is A bit towards back in the beginning,
%% Note: and then gradually move to NORMAL?. After the 2nd StrLearn

is.fnDate{is.nSubj} = '2018-02-24'; is.fnSID{is.nSubj} = '11';
is.fnBehav{is.nSubj} = {'0011'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05221_Yunzhe_20180224'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'};  
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [12]
%% Note: Head position - normal
%% Note: Task order - 3
%% Note: Training - Take 1  Rounds to finish (5  sessions)
%% Note: Age: 31
%% Note: didn't drink coffee, signal is a little noisy…, lose eye tracker in the 2nd resting

is.fnDate{is.nSubj} = '2018-03-12'; is.fnSID{is.nSubj} = '12';
is.fnBehav{is.nSubj} = {'0012'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05230_Yunzhe_20180312'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [13]
%% Note: Head position - normal, a bit of frontal
%% Note: Task order - 3
%% Note: Training - Take 1  Rounds to finish (5  sessions)
%% Note: Age: 25
%% Note: lab member, 
%% Note: didn't drink coffee, signal is a little noisy, 
%% Note: eye tracker very unstable, wear contact lens, 
%% Note: a lot of musle atrifects - want to go to the loo, so a lot of musle tension

is.fnDate{is.nSubj} = '2018-03-13'; is.fnSID{is.nSubj} = '13';
is.fnBehav{is.nSubj} = {'0013'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05231_Yunzhe_20180313'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [14] 
%% Note: Head position - normal, a bit of frontal
%% Note: Task order - 3
%% Note: Training - Take 1  Rounds to finish (5  sessions)
%% Note: Age: 24, (indian)
%% Note: Professional participants,
%% Note: drink black tea before the experiment; 
%% Note: More head movemnt and eye blink than previous participants 

is.fnDate{is.nSubj} = '2018-03-16'; is.fnSID{is.nSubj} = '14';
is.fnBehav{is.nSubj} = {'0014'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05237_Yunzhe_20180316'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [15]  left hand
%% Note: Head position - normal, a bit of frontal
%% Note: Task order - 3
%% Note: Training - Take 1  Rounds to finish (5  sessions)
%% Note: Age: 26
%% Note: drink tea,  
%% Note: Less than 70% on average after the applied learning
%% Note: In the position code, have a huge 200mm once and more than 95% fitting error, but move back!

is.fnDate{is.nSubj} = '2018-03-19'; is.fnSID{is.nSubj} = '15';
is.fnBehav{is.nSubj} = {'0015'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05239_Yunzhe_20180319'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [16] 
%% Note: Head position - normal
%% Note: Task order - 3
%% Note: Training - Take 2 Rounds to finish (10  sessions)
%% Note: Age: 24
%% Note: drink tea,  
%% Note: Less than 70% after the applied learning, around 33% ACC after first run of applied learning
%% Note: lost eye tracker in the middle of the last run - position code

is.fnDate{is.nSubj} = '2018-03-19'; is.fnSID{is.nSubj} = '16';
is.fnBehav{is.nSubj} = {'0016'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05240_Yunzhe_20180319'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [17] 
%% Note: Head position - normal
%% Note: Task order - 3
%% Note: Training - Did not finish
%% Note: Age: 20
%% Note: Quite impatient, late for 30min;
%% Note: Noisy eye movement; wear contact lens.
%% Note: compressor error - maynot have good data in the first session (preplay), 
%% Note: maybe problem with saving the second session; 
%% Note: Lost eye tracker in the middle of the 9th  run - seq  code

is.fnDate{is.nSubj} = '2018-03-20'; is.fnSID{is.nSubj} = '17';
is.fnBehav{is.nSubj} = {'0017'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05242_Yunzhe_20180320'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [18] 
%% Note: Head position - normal? a bit of frontal - quite frontal actually
%% Note: Task order - 4
%% Note: Training - Take 2 Rounds to finish (10  sessions)
%% Note: Age: 23
%% Note: lab member, 
%% Note: eye tracker data is really noisy!
%% Note: quite a lot of musle artifects , 
%% Note: chair location maybe wrong - but this is as the best as we can possibly get, 

is.fnDate{is.nSubj} = '2018-03-21'; is.fnSID{is.nSubj} = '18';
is.fnBehav{is.nSubj} = {'0018'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05244_Yunzhe_20180321'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [19] 
%% Note: Head position - normal
%% Note: Task order - 4
%% Note: Training - Did not finish
%% Note: Age: 28? - cannot be sure
%% Note: very causal - have difficultity in following instructions, 
%% Note: wear makeups, seems tired - noisy eye movement
%% Note: Preplay: huge movement once - may not be a good data
%% Note: lose eye tracker data in the second run of applied learning

is.fnDate{is.nSubj} = '2018-03-22'; is.fnSID{is.nSubj} = '19';
is.fnBehav{is.nSubj} = {'0019'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05248_Yunzhe_20180322'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [20] 
%% Note: Head position - normal
%% Note: Task order - 4
%% Note: Training - Did not finish
%% Note: Age: 30
%% Note: Environmental noise, due to the coil; 
%% Note: short sighted; 
%% Note: a little nervious, want to pee, but hold up to the end, 
%% Note: close eyes a lot, and some huge movement occocassionly,  
%% Note: A LOT OF MUSULE MOVEMENT!, 
%% Note: EYE TRACKER DATA MESSY! HUGE ALPHA & THINK A LOT BEFORE ANY TESTING…

is.fnDate{is.nSubj} = '2018-03-23'; is.fnSID{is.nSubj} = '20';
is.fnBehav{is.nSubj} = {'0020'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05251_Yunzhe_20180323'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [21] 
%% Note: Head position - normal
%% Note: Task order - 4
%% Note: Training - Take 2  Rounds to finish (10  sessions)
%% Note: Age: 27
%% Note: Environmental noise, due to the coil; 
%% Note: A LOT OF MUSULE MOVEMENT; 
%% Note: noisy eyement 

is.fnDate{is.nSubj} = '2018-03-24'; is.fnSID{is.nSubj} = '21';
is.fnBehav{is.nSubj} = {'0021'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05252_Yunzhe_20180324'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [22]  the first to reach 90%
%% Note: Head position - normal
%% Note: Task order - 4
%% Note: Training - Take 2  Rounds to finish (10  sessions)
%% Note: Age: 20
%% Note: Environmental noise, due to the coil; 
%% Note: not very enthusiastic
%% Note: one huge head movement >30mm in the second resting (replay)

is.fnDate{is.nSubj} = '2018-03-24'; is.fnSID{is.nSubj} = '22';
is.fnBehav{is.nSubj} = {'0022'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05253_Yunzhe_20180324'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [23]  the ACC in applied learning <60%
%% Note: Head position - normal, a bit of frontal
%% Note: Task order - 4
%% Note: Training - Take 1 Rounds to finish (5  sessions)
%% Note: Age: 20
%% Note: Environmental noise, due to the coil; 
%% Note: noisy eyement

is.fnDate{is.nSubj} = '2018-03-25'; is.fnSID{is.nSubj} = '23';
is.fnBehav{is.nSubj} = {'0023'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05254_Yunzhe_20180325'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [24]  
%% Note: Head position - normal, a bit of frontal
%% Note: Task order - 4
%% Note: Training - Take 2 Rounds to finish (10  sessions)
%% Note: Age: 27
%% Note: Environmental noise, due to the coil; 
%% Note: a lot OF MUSULE MOVEMENT

is.fnDate{is.nSubj} = '2018-03-25'; is.fnSID{is.nSubj} = '24';
is.fnBehav{is.nSubj} = {'0024'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05255_Yunzhe_20180325'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [25] 
%% Note: Head position - normal
%% Note: Task order - 5
%% Note: Training - Take 1 Rounds to finish (5  sessions)
%% Note: Age: 30
%% Note: A SUDDEN fitting error in PREPLAY!!!
%% Note: noise eye tracker in the pos coding for 2min, then been asked to sit a bit up - but cause head movement to >10mm

is.fnDate{is.nSubj} = '2018-04-07'; is.fnSID{is.nSubj} = '25';
is.fnBehav{is.nSubj} = {'0025'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05269_Yunzhe_20180407'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [26]  (pretty bad!, ACC <50%)
%% Note: Head position - normal
%% Note: Task order - 5
%% Note: Training - Take 2  Rounds to finish (10  sessions) - BUT WRITE IT DOWN!
%% Note: Age: 33
%% Note: old + poor memory, 
%% Note: go the bathroom after the functional localizer, 
%% Note: always move, huge movement in the second resting!!!

is.fnDate{is.nSubj} = '2018-04-08'; is.fnSID{is.nSubj} = '26';
is.fnBehav{is.nSubj} = {'0026'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05270_Yunzhe_20180408'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [27]  
%% Note: Head position - normal
%% Note: Task order - 5
%% Note: Training - Take 2 Rounds to finish (10  sessions) 
%% Note: Age: 23
%% Note: eye movement is a little messy
%% Note: go the bathroom after the functional localizer, 
%% Note: always move, huge movement in the second resting!!!

is.fnDate{is.nSubj} = '2018-04-15'; is.fnSID{is.nSubj} = '27';
is.fnBehav{is.nSubj} = {'0027'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05283_Yunzhe_20180415'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [28]  
%% Note: Head position - normal
%% Note: Task order - 6
%% Note: Training - Did not finish
%% Note: Age: 23
%% Note: lab member, 
%% Note: eye movement is a little messy
%% Note: nervous, drink a lot of coffee
%% Note: DID NOT save the 2-min resting after the SECOND & THIRD learning session...

is.fnDate{is.nSubj} = '2018-04-19'; is.fnSID{is.nSubj} = '28';
is.fnBehav{is.nSubj} = {'0028'}; is.fnTrainClass{is.nSubj} = {'Train_0001', 'StrReplay_0001'};
is.fnMEG{is.nSubj} = 'MG05291_Yunzhe_20180419'; 
is.MEGruns{is.nSubj} = {'rst' 'lci' 'lci' 'lci' 'str' 'str' 'str' 'rst' 'seq' 'pos'}; 
is.taskversion{is.nSubj}=2; 
is.nSubj = is.nSubj + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



is.nSubj = is.nSubj - 1;