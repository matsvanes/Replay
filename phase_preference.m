% This script loads in probabilities of the decoding model applied to
% resting state data. Based on source parcellations, the hypothesis is
% tested whether the highest probabilities are phase clustered (e.g. in the
% theta/alpha band), w.r.t. the lowest probabilities (this is a good
% control because classifiers have a trivial phase preference).
% 2021 - Mats van Es

%% Define some variables
whichstudy = 2;
if whichstudy==1
  MFStartup;
else
  MFStartup_studyII;
  load('/Volumes/T5_OHBA/data/replay/StrLearn_MEGexp/BehavData/Exptrial.mat')
  is.dirname='StrLearn_MEGexp/';
  savebase = '/Users/matsvanes/Data/YunzheData/PhasePreference/StrLearn_MEGexp/';
end
datadir = [is.studydir,'Neuron2020Analysis/Study2/bfnew_1to45hz/'];
if ~isfolder(is.AnalysisPath),   mkdir(is.AnalysisPath), end
is.usePrecomputed = true;
is.iRun = 2;
is.topPercentile=1;
fsample=250;
t_window=fsample/2;
parc=parcellation('fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
FOI = 'alpha'; % can be 'alpha', 'theta'

%% Get alpha peak frequency
if strcmp(FOI,'alpha')
  for iSj=1:length(is.goodsub)
    sprintf('Subject %d', iSj)
    RST = spm_eeg_load([datadir, sprintf('sfold_giles_symmetric_f_session%d', iSj*2)]);
    doplot=0;
    [~, peakfreq(iSj)] = compute_peakfreq(RST, 2:0.1:20, [],[], doplot);
  end
end

%% Get phases
[~,~,triggerpoints,goodsamples]=getSubjectMasks(datadir);
triggerpoints = triggerpoints(2:2:end); % only need post-task data
goodsamples = goodsamples(2:2:end);

for iSj=1:length(is.goodsub)
  sprintf('Subject %d', iSj)
  RST = spm_eeg_load([datadir, sprintf('sfold_giles_symmetric_f_session%d', iSj*2)]);

  % get FOI phase
  if strcmp(FOI, 'alpha')
    [phase{iSj}, pow{iSj}, maxidx(iSj)] = getphasetimecourse(RST, peakfreq(iSj));
  elseif strcmp(FOI,'theta')
    [phase{iSj}, pow{iSj}, maxidx(iSj)] = getphasetimecourse(RST, 6,'dpss',2);
  end
  phase{iSj} = phase{iSj}(:,triggerpoints{iSj});
  pow{iSj} = pow{iSj}(:,triggerpoints{iSj});
end



%% phase preference for reactivation (top vs bottom percentile)
for iSj=1:length(is.goodsub)
  tmp=load([is.AnalysisPath, 'classifiers/TestResting4Cell/', sprintf('prob%d_%d',iSj,is.iRun)]);
  probAll{iSj} = tmp.prob;
end

for iSj=1:length(is.goodsub)
  % take the phase of these points
  [topphase{iSj}, botphase{iSj}] = get_percentilePhase(probAll{iSj}', phase{iSj}, t_window, is.topPercentile);

  % compare concentration values
  for iParc = 1:38
    % [~,fval(iSj, iParc)] = circ_ktest(topphase{iSj}(iParc,:), botphase{iSj}(iParc,:));
    % compute non-uniformity of the phase data (Rayleigh test)
    [~,z_top(iSj,iParc)]=circ_rtest(topphase{iSj}(iParc,:));
    [~,z_bot(iSj,iParc)]=circ_rtest(botphase{iSj}(iParc,:));
  end
end
GA_diff = mean(z_top-z_bot);
filename = [savebase 'phasepreference_reactivation_', sprintf('%s', FOI)];
% parc.savenii(GA_diff, filename);
% osl_render4D([filename, '.nii.gz'], 'visualise', true);

f1=figure('WindowState','maximized');
plot_phaseConcentration(topphase, botphase, z_top, z_bot, FOI,'reactivation')
saveas(f1,[filename, '.png'])
% stat = stat_phasePreference(z_top, z_bot);

%% phase preference for replay (top vs bottom percentile)
clear z_top z_bot botphase topphase GA_diff
if whichstudy==1
  load([wd,'GenericReplayData/STUDYI_ReplayOnset/RwdReplay2ndOnset'],'ToRall');
  replayScores = ToRall;
else
  load([is.AnalysisPath,'classifiers/Sequence_by_Training_4Cell/', 'STUDYII_ReplayOnset'],'ToRall');
  replayScores = ToRall;
end

for iSj=1:length(is.goodsub)
  % take the phase of these points
  [topphase{iSj}, botphase{iSj}] = get_percentilePhase(replayScores(iSj,:), phase{iSj}, t_window, is.topPercentile);

  for iParc = 1:38
    % [~,fval(iSj, iParc)] = circ_ktest(topphase{iSj}(iParc,:), botphase{iSj}(iParc,:));
    % compute non-uniformity of the phase data (Rayleigh test)
    [~,z_top(iSj,iParc)]=circ_rtest(topphase{iSj}(iParc,:));
    [~,z_bot(iSj,iParc)]=circ_rtest(botphase{iSj}(iParc,:));
  end
end
GA_diff = mean(z_top-z_bot);
filename = [savebase 'phasepreference_replay_', sprintf('%s', FOI)];
% parc.savenii(GA_diff, filename);
% osl_render4D([filename, '.nii.gz'], 'visualise', true);

f2=figure('WindowState','maximized');
plot_phaseConcentration(topphase, botphase, z_top, z_bot, FOI,'replay')
saveas(f2,[filename, '.png'])
%% Legacy
%{
    for iC=1:nStim
      % use regression models
      Rreds{iRun}(:,iC) = -(data * gnf{iC}.beta + gnf{iC}.Intercept);   % work on the X*beta space! MVE: the minus is because of a weird matlab definition. Now high positive values mean high probability
    end
    prob{iRun} = 1./(1+exp(Rreds{iRun})); % the Rreds are the negative of how they usually are defined. That's why we're using +Rreds instead of minus (which is normally in the logsigmoid function)
    
    % find highest percentile probabilities
    x=prob{iRun}; x(isnan(x))=-Inf;
    [val, idx] = sort(prob{iRun}, 'descend', 'MissingPlacement', 'last');
    n = round((size(val,1)-sum(isnan(data(:,1))))*perc/100);
    idx_top = idx(1:n,:);
    val_top{iRun} = val(1:n,:);

    
    
    % find lowest percentile probabilities
    [val, idx] = sort(prob{iRun}, 'ascend', 'MissingPlacement', 'last');
    n = round((size(val,1)-sum(isnan(data(:,1))))*perc/100);
    idx_bottom = idx(1:n,:);
    val_bottom{iRun} = val(1:n,:);
    
    % Get alpha phase from parcel time courses.
    S=[];
    S.outfile = fullfile(dir_temp, 'parcelphase');
    S.D = RST;
    RST_copy = spm_eeg_copy(S);
    [phase, pow, maxidx] = getphasetimecourse(RST_copy);
    phase = phase(maxidx, :); % take the parcel with maximum power
    % TODO: Continue here!
    time = 0:1/RST.fsample:(size(data,1)-1)/RST.fsample;
    f=10;
    clear i
    a = cos(2*pi*f*time) + i*sin(2*pi*f*time); % raw alpha signal
    phase = angle(a);
    phase = repmat(phase, [8 1])';
    for k=1:8
      % take the phases of the highest/lowest probabilities
      phase_top{iRun}(:,k) = phase(idx_top(:,k));
      phase_bottom{iRun}(:,k) = phase(idx_bottom(:,k));
      % probably we can use the parametric test, but this assumes max 1
      % mode, and a "gaussian" distribution around it.
      % alternatives are the nonparametric Hodges-Ajne test (circ_otest),
      % but this doesn't return a test statistic. Other options is to
      % compare the distribution with a uniform distribution using the
      % Kuiper test (circ_kuipertest) or Winson's U2 test.
      [pval(iRun,k), z(iRun,k)] = circ_rtest(phase_top{iRun}(:,k));
      %       [pval(iRun,k), m(iRun,k)] = circ_otest(phase_top{iRun}(:,k));
      subplot(2,8,(iRun-1)*8+k);
      rose(phase_top{iRun}(:,k)); hold on
      rose(phase_bottom{iRun}(:,k));
      if (iRun-1)*8+k == 4
        title('pre-task Resting State')
      elseif (iRun-1)*8+k == 12
        title('post-task Resting State')
      end
    end
    for k=1:8
      [~,fval(iRun,k)] = circ_ktest(phase_top{iRun}(:,k), phase_bottom{iRun}(:,k));
    end
  end
  % optionally - use Watson's U2 statistic to compare the two circular
  % distributions. See https://uk.mathworks.com/matlabcentral/fileexchange/43543-pierremegevand-watsons_u2
  % or possibly use circ_ktest for equal concentration parameter?
end
%}
