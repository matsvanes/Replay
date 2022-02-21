function gnfperm = permute_classifierweights(gnf, is, randstate)
% This function loops over lambdas, and permutes the classifier weights
% (betas) over stimuli in a pairwise fassion, i.e. the weight for
% MEG channel 1 is permuted over stimuli. Permutation is the same for each
% lambda.

if exist('randstate','var') && ~isempty(randstate)
  rng(randstate);
end

nChan = length(gnf{1}.beta);
perms = zeros(nChan,is.nStim);
for iChan=1:nChan
  perms(iChan,:) = randperm(is.nStim);
end

gnfperm=gnf;
for iL=1:size(gnf,2)
  for iStim=1:is.nStim
    for iChan=1:nChan
      gnfperm{iStim,iL}.beta(iChan,1) = gnf{perms(iChan,iStim),iL}.beta(iChan,1);
    end
  end
end
