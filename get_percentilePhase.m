function [topphase, botphase] = get_percentilePhase(prob, phase, t_window, percentile)

% get top percentile probabilities
topperc = prob > prctile(prob(:),100-percentile);
topperc = [zeros(size(topperc,1),1),diff(topperc,[],2)]==1;% eliminate adjacent points
topperc(:,1:t_window)=0;topperc(:,end-t_window:end)=0;% eliminate border points:
topidx = upsample_replayTimes(topperc, phase, t_window);

if nargout>1
  % get bottom percentile probabilities
  botperc = prob < prctile(prob(:),percentile);
  botperc = [zeros(size(topperc,1),1),diff(botperc,[],2)]==1;% eliminate adjacent points
  botperc(:,1:t_window)=0;botperc(:,end-t_window:end)=0;
  botidx = upsample_replayTimes(botperc, phase, t_window);


  % get same amount of points for replay and control
  n = min([numel(botidx), numel(topidx)]);
  ix1=randperm(numel(topidx));
  ix2=randperm(numel(botidx));
  topidx = topidx(ix1(1:n));
  botidx = botidx(ix2(1:n));
  
  % take the phase of these points
  botphase = phase(:,botidx);
end
topphase = phase(:,topidx);
