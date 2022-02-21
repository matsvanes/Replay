function [topphase, ZtopPhase] = phase_nonuniformity(probability, phase, is, fname, mask)

if is.usePrecomputed
  try
    tmp = load(fname);
    topphase = tmp.topphase;
    ZtopPhase = tmp.ZtopPhase;
  catch
    is.usePrecomputed = false;
    if ~is.doPermutation
      fprintf('Loading in data was not possible - recomputing now \n')
    end
    [topphase, ZtopPhase] = phase_nonuniformity(probability, phase, is, fname, mask);
  end
else
  if ~is.doPermutation
    fprintf('Computing phase non-uniformity \n')
  end
  topphase = get_percentilePhase(probability', phase, is.t_window, is.topPercentile, mask);

  [~,ZtopPhase] = circ_rtest_multivar(topphase, [], [], 3);

  if exist('fname', 'var') && ~isempty(fname)
    save(fname, 'topphase', 'ZtopPhase')
  end
end
