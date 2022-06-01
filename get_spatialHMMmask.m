% get spatial mask - using net_mean from NeuronFig2Analyses.m
function mask = get_spatialHMMmask
K=12;
template_string='5usingtemplate';
hmmdir='/Users/matsvanes/Data/YunzheData/Neuron2020Analysis/Study2/hmm_1to45hz/';
nnmf_outfileWB = fullfile( hmmdir, ['embedded_HMM_K',int2str(K),template_string,'_nnmfWB']);
load(nnmf_outfileWB)

net_mean = zeros(38,K);
for k = 1:K
    net_mean(:,k) = squeeze(nnmfWB_res.nnmf_psd_maps(k,1,:))';
end

mask=zeros(38,12);
for kk=1:12
  toplot = net_mean(:,kk);%-mean(net_mean,2);
  psdthresh = prctile(abs(toplot),50);
  %CL = psdthresh*[-1.25,1.25];
  %CL = [min(toplot(:)), max(toplot(:))];
  toplot(abs(toplot)<psdthresh) = NaN;
  toplot(~isnan(toplot))=1;
  mask(:,kk)=toplot;
end