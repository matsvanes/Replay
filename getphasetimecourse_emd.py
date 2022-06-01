import emd
import mne
import numpy as np
import osl
from scipy.signal import butter, lfilter
from scipy.io import savemat
from scipy import ndimage, stats
import pickle
import scipy
import matplotlib.pyplot as plt
import sails

sample_rate = 250
max_imfs = 4
nbins = 49
W = 125  # 0.5 seconds at 250 Hz
t = np.arange(-0.5, 0.5 + 1 / 250, 1 / 250)


# %%
def butter_lowpass(cutoff, sample_rate, order=5):
    nyq = 0.5 * sample_rate
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, sample_rate, order=5):
    b, a = butter_lowpass(cutoff, sample_rate, order=order)
    y = lfilter(b, a, data)
    return y


for iSes in range(1, 44, 2):
    D = osl.utils.spmio.SPMMEEG(
        f"/Users/matsvanes/Data/YunzheData/Neuron2020Analysis/Study2/bfnew_1to45hz/sfold_giles_symmetric_f_session{iSes + 1}.mat")
    time_vect = D.time
    IMF = np.zeros((38, len(time_vect), max_imfs))
    IP = np.zeros((38, len(time_vect), max_imfs))
    IF = np.zeros((38, len(time_vect), max_imfs))
    IA = np.zeros((38, len(time_vect), max_imfs))
    HHT = np.zeros((38, nbins, len(time_vect)))
    for iParc in range(0, 38):
        print(f"Subject {(iSes + 1) / 2}, Parcel {iParc + 1}")
        x = D.get_data()[iParc, :]
        x = butter_lowpass_filter(x, 30, sample_rate)

        # Iterative Masked EMD
        tmpIMF = emd.sift.iterated_mask_sift(x, sample_rate=sample_rate, max_imfs=max_imfs)
        tmpIP, tmpIF, tmpIA = emd.spectra.frequency_transform(tmpIMF, sample_rate, 'nht')
        f, tmpHHT = emd.spectra.hilberthuang(tmpIF, tmpIA, (1, 50, nbins), sum_time=False)

        # collect results for all parcels
        IMF[iParc, :, :] = tmpIMF
        IP[iParc, :, :] = tmpIP
        IF[iParc, :, :] = tmpIF
        IA[iParc, :, :] = tmpIA
        HHT[iParc, :, :] = tmpHHT
        spctrm = scipy.signal.hilbert(IMF, axis=1)

    savemat(f"/Users/matsvanes/Data/YunzheData/Mats/Study2/phase/emd_session{iSes + 1}.mat",
            {"imf": IMF, "IP": IP, "IF": IF, "IA": IA, "hht": HHT, "f": f, "spctrm": spctrm})

# %%
for iSes in range(1, 44, 2):
    Q = scipy.io.loadmat(f"/Users/matsvanes/Data/YunzheData/Mats/Study2/phase/emd_session{iSes + 1}.mat")
    imf = Q['imf']
    spctrm = scipy.signal.hilbert(imf, axis=1)
    savemat(f"/Users/matsvanes/Data/YunzheData/Mats/Study2/phase/emd_spctrm_session{iSes + 1}.mat", {"spctrm": spctrm})

# %%
tmp = scipy.io.loadmat('/Users/matsvanes/Data/YunzheData/Mats/Study2/Triggerpoints.mat')
triggerpoints = np.squeeze(tmp['triggerpoints'])
tmp = scipy.io.loadmat('/Users/matsvanes/Data/YunzheData/Mats/Study2/replay_idx_mask1.mat')
topidx = np.squeeze(tmp['topidx'])
W = 125
nsubs = 22
group = {"imf_tl": np.zeros((38, 4, 251, 22)), "IF_tl": np.zeros((38, 4, 251, 22)), "IA_tl": np.zeros((38, 4, 251, 22)),
         "IP_tl": np.zeros((38, 4, 251, 22)), "hht_tl": np.zeros((38, nbins, 251, 22)),
         'imf': np.zeros((38, 4, 75000, nsubs)), 'IF': np.zeros((38, 4, 75000, nsubs)),
         'IP': np.zeros((38, 4, 75000, nsubs)), 'IA': np.zeros((38, 4, 75000, nsubs)),
         'hht': np.zeros((38, nbins, 75000, nsubs))}
for iSes in range(1, 44, 2):
    iSj = int((iSes + 1) / 2)
    print(f'{iSj}')
    Q = scipy.io.loadmat(f"/Users/matsvanes/Data/YunzheData/Mats/Study2/phase/emd_session{iSes + 1}.mat")
    imf = Q['imf'][:, np.where(np.concatenate(triggerpoints[iSes]))[0], :]
    IP = Q['IP'][:, np.where(np.concatenate(triggerpoints[iSes]))[0], :]
    IF = Q['IF'][:, np.where(np.concatenate(triggerpoints[iSes]))[0], :]
    IA = Q['IA'][:, np.where(np.concatenate(triggerpoints[iSes]))[0], :]
    hht = Q['hht'][:, :, np.where(np.concatenate(triggerpoints[iSes]))[0]]
    f = Q['f']
    group['imf'][:, :, :, int(iSj - 1)] = np.transpose(imf, (0, 2, 1))
    group['IA'][:, :, :, int(iSj - 1)] = np.transpose(IA, (0, 2, 1))
    group['IF'][:, :, :, int(iSj - 1)] = np.transpose(IF, (0, 2, 1))
    group['IP'][:, :, :, int(iSj - 1)] = np.transpose(IP, (0, 2, 1))
    group['hht'][:, :, :, int(iSj - 1)] = hht

    # create epochs
    idx = np.squeeze(topidx[int(iSj - 1)])
    trls = np.zeros((len(idx), 2))
    for i in np.arange(len(idx)):
        trls[i, 0] = int(idx[i] - W)
        trls[i, 1] = int(idx[i] + W)

    imf_tl = np.zeros((38, 4, 2 * W + 1, len(idx)))
    IF_tl = np.zeros((38, 4, 2 * W + 1, len(idx)))
    IP_tl = np.zeros((38, 4, 2 * W + 1, len(idx)))
    IA_tl = np.zeros((38, 4, 2 * W + 1, len(idx)))
    hht_tl = np.zeros((38, nbins, 2 * W + 1, len(idx)))
    mask = np.zeros((75000))
    for i in np.arange(trls.shape[0]):
        mask[int(trls[i, 0]):int(trls[i, 1] + 1)] = 1
        imf_tl[:, :, :, i] = np.transpose(imf[:, int(trls[i, 0]):int(trls[i, 1] + 1), :], axes=[0, 2, 1])
        IF_tl[:, :, :, i] = np.transpose(IF[:, int(trls[i, 0]):int(trls[i, 1] + 1), :], axes=[0, 2, 1])
        IP_tl[:, :, :, i] = np.transpose(IP[:, int(trls[i, 0]):int(trls[i, 1] + 1), :], axes=[0, 2, 1])
        IA_tl[:, :, :, i] = np.transpose(IA[:, int(trls[i, 0]):int(trls[i, 1] + 1), :], axes=[0, 2, 1])
        hht_tl[:, :, :, i] = hht[:, :, int(trls[i, 0]):int(trls[i, 1] + 1)]

    group['imf_tl'][:, :, :, int(iSj - 1)] = np.mean(imf_tl, axis=3)
    group['IA_tl'][:, :, :, int(iSj - 1)] = np.mean(IA_tl, axis=3)
    group['IF_tl'][:, :, :, int(iSj - 1)] = np.mean(IF_tl, axis=3)
    group['IP_tl'][:, :, :, int(iSj - 1)] = emd.imftools.ip_circular_mean(IP_tl, axis=3)
    group['hht_tl'][:, :, :, int(iSj - 1)] = np.mean(hht_tl, axis=3)

    # mask_cycles = emd.cycles.get_cycle_vector(IP[0, :, :], return_good=False, mask=mask > 0)

filehandler = open("".join(('/Users/matsvanes/Data/YunzheData/Mats/Study2/emd/', "emd_group.p")), "wb")
pickle.dump({"group": group}, filehandler)
filehandler.close()

# Look at frequency difference from mean
tmp = scipy.io.loadmat("/Users/matsvanes/Documents/Werk/scripts/Replay/Giles38_Adjacency.mat")
adjacency = tmp['adjacency']
IF_tl_normalized = group["IF_tl"] - np.transpose(np.tile(np.expand_dims(np.mean(group["IF_tl"], axis=2), axis=3), 251),
                                                 [0, 1, 3, 2])
D = osl.utils.spmio.SPMMEEG(
    f"/Users/matsvanes/Data/YunzheData/Neuron2020Analysis/Study2/bfnew_1to45hz/sfold_giles_symmetric_f_session2.mat")

info = mne.create_info(ch_names=[D.channels[k].label for k in range(len(D.channels))], ch_types=['dipole'] * 38,
                       sfreq=D.fsample)
null = IF_tl_normalized * 0

X = np.transpose(IF_tl_normalized, [3, 2, 1, 0])
for i in range(3):
    T_obs, clusters, cluster_p_values, H0 = clu = \
        mne.stats.spatio_temporal_cluster_1samp_test(X[:, :, i, :], adjacency=scipy.sparse.csr_matrix(adjacency),
                                                     n_jobs=1, tail=0, n_permutations=1000, buffer_size=None,
                                                     verbose=True, out_type='mask')
    good_cluster_inds = np.where(cluster_p_values < 0.05)[0]
    if len(good_cluster_inds) > 0:
        for k in range(len(good_cluster_inds)):
            plt.figure()
            plt.pcolormesh(t, np.arange(38), np.transpose(clusters[good_cluster_inds[k]]))
            plt.title(f'IMF {i + 1} ({np.round(np.mean(group["IF_tl"][:, i, :, :]))} Hz) - cluster {k + 1}')

plt.figure()
for k in np.arange(4):
    plt.subplot(2, 2, k + 1)
    plt.plot(t, np.transpose(np.mean(
        group["IF_tl"] - np.transpose(np.tile(np.expand_dims(np.mean(group["IF_tl"], axis=2), axis=3), 251),
                                      [0, 1, 3, 2]), axis=3)[:, k, :]))
    plt.plot(t, np.mean(np.transpose(np.mean(
        group["IF_tl"] - np.transpose(np.tile(np.expand_dims(np.mean(group["IF_tl"], axis=2), axis=3), 251),
                                      [0, 1, 3, 2]), axis=3)[:, k, :]), axis=1), 'k', linewidth=2)
    plt.title(f'IMF {k + 1} ({np.round(np.mean(group["IF_tl"][:, k, :, :]))} Hz)')
    plt.ylabel('IF w.r.t. mean')
    plt.xlabel('Time (s)')

fig = plt.figure()
vmin = [-0.8, -0.3, -0.15, -0.1]
vmax = [0.8, 0.3, 0.15, 0.1]
for k in np.arange(4):
    ax = plt.subplot(2, 2, k + 1)
    c = plt.pcolormesh(t, np.arange(38), np.mean(
        group["IF_tl"] - np.transpose(np.tile(np.expand_dims(np.mean(group["IF_tl"], axis=2), axis=3), 251),
                                      [0, 1, 3, 2]), axis=3)[:, k, :], cmap='inferno', vmin=vmin[k], vmax=vmax[k])
    # plt.pcolormesh(t, np.arange(38), np.mean(
    #    group["IF_tl"], axis=3)[:, k, :])
    plt.title(f'IMF {k + 1} ({np.round(np.mean(group["IF_tl"][:, k, :, :]))} Hz)')
    plt.ylabel('parcel')
    plt.xlabel('Time (s)')
    fig.colorbar(c, ax=ax)

plt.plot(t, np.transpose(emd.imftools.ip_circular_mean(group["IF_tl"], axis=3)[0, :, :]))

##% phase alignment
tmp = scipy.io.loadmat('/Users/matsvanes/Data/YunzheData/Mats/Study2/replay_idx_mask1.mat')
topidx = np.squeeze(tmp['topidx'])
q = []
# Iterate over a sequence of numbers from 0 to 4
for i in range(22):
    # In each iteration, add an empty list to the main list
    q.append([])
    for j in range(38):
        q[i].append([])

for nimf in np.arange(2):
    waveformshape = {'pa_if_noreplay_mean': np.empty((38, 48, nsubs)), 'pa_if_noreplay_std': np.empty((38, 48, nsubs)),
                     'pa_if_replay_mean': np.empty((38, 48, nsubs)), 'pa_if_replay_std': np.empty((38, 48, nsubs)),
                     'pa_if_1stlvl': np.empty((38, 48, nsubs)), 'pa_r_mean': np.empty((38, 48, nsubs)),
                     'pa_r_std': np.empty((38, 48, nsubs)), 'df': q.copy(), 'Comp1_1stlvl': np.empty((38,nsubs)),
                     'Comp2_1stlvl': np.empty((38,nsubs)), 'Comp3_1stlvl': np.empty((38,nsubs))}
    for iSj in np.arange(nsubs):
        print(iSj + 1)
        for iparc in np.arange(38):
            IF = group['IF'][iparc, :, :, iSj]
            IP = group['IP'][iparc, :, :, iSj]
            IA = group['IA'][iparc, :, :, iSj]

            replays = np.squeeze(topidx[int(iSj - 1)])
            replay_ts = np.zeros_like(IF[0, :])
            replay_ts[replays] = 1

            # little bit of smoothing in replay onset timing
            replay_ts = ndimage.maximum_filter(replay_ts, 5)

            C = emd.cycles.Cycles(IP[nimf, :])

            # Compute per-cycle metrics
            C.compute_cycle_metric('max_amp', IA[nimf, :], np.max)
            C.compute_cycle_metric('duration', IA[nimf, :], len)
            C.compute_cycle_metric('mean_if', IF[nimf, :], np.mean)
            C.compute_cycle_metric('has_replay', replay_ts, np.max)

            conditions = ['is_good==1', 'max_amp>0.002', 'mean_if<15']
            df = C.get_metric_dataframe(conditions=conditions)

            # Phase align IF
            Citer = C.iterate(conditions=conditions)
            pa_if, phasex = emd.cycles.phase_align(IP[nimf, :], IF[nimf, :], cycles=Citer)
            pa_if = pa_if[:, :-1]  # Extra cycle? probably a bug...

            # metrics regarding phase aligned IF
            norep = pa_if[:, df['has_replay'] == False]
            waveformshape['pa_if_noreplay_mean'][iparc, :, iSj] = norep.mean(axis=1)
            waveformshape['pa_if_noreplay_std'][iparc, :, iSj] = norep.std(axis=1)
            rep = pa_if[:, df['has_replay'] == True]
            waveformshape['pa_if_replay_mean'][iparc, :, iSj] = rep.mean(axis=1) if rep.shape[
                                                                                        1] > 0 else np.nan * np.ones(48)
            waveformshape['pa_if_replay_std'][iparc, :, iSj] = rep.std(axis=1) if rep.shape[
                                                                                      1] > 0 else np.nan * np.ones(48)
            waveformshape['pa_if_1stlvl'][iparc, :, iSj], _ = stats.ttest_ind(rep, norep, axis=1, equal_var=False)

            # Phase align timing of replay within cycles containing replay
            Citer = C.iterate(conditions=['has_replay>0'])
            pa_r, phasex = emd.cycles.phase_align(IP[nimf, :], replay_ts, cycles=Citer)
            waveformshape['pa_r_mean'][iparc, :, iSj] = pa_r.mean(axis=1)
            waveformshape['pa_r_std'][iparc, :, iSj] = pa_r.std(axis=1)

            # Compute PCA across phase aligned cycles
            pca = sails.utils.PCA(pa_if.T, np.min((pa_if.shape[1] - 1, 10)))

            # Add PCA compoent scores to dataframe
            df['Comp1'] = pca.scores[:, 0]  # continuum of sharp peak and wide through to wide peak and sharp through
            df['Comp2'] = pca.scores[:, 1]  # elongated ascending edge and an elongated descending edge
            df['Comp3'] = pca.scores[:, 2] if pca.scores.shape[1] >= 3 else np.nan * np.ones(
                pa_if.shape[1])  # pca.scores[:, 2]  # shapes with a left or right “tilt” around their extrema
            waveformshape['df'][iSj][iparc] = df

            repidx = np.where(waveformshape['df'][iSj][iparc]['has_replay'] == 1)[0]
            norepidx = np.where(waveformshape['df'][iSj][iparc]['has_replay'] == 0)[0]
            waveformshape['Comp1_1stlvl'][iparc, iSj], _ = stats.ttest_ind(
                waveformshape['df'][iSj][iparc]['Comp1'][repidx].values,
                waveformshape['df'][iSj][iparc]['Comp1'][norepidx].values, axis=0, equal_var=False)
            waveformshape['Comp2_1stlvl'][iparc, iSj], _ = stats.ttest_ind(
                waveformshape['df'][iSj][iparc]['Comp2'][repidx].values,
                waveformshape['df'][iSj][iparc]['Comp2'][norepidx].values, axis=0, equal_var=False)
            waveformshape['Comp3_1stlvl'][iparc, iSj], _ = stats.ttest_ind(
                waveformshape['df'][iSj][iparc]['Comp3'][repidx].values,
                waveformshape['df'][iSj][iparc]['Comp3'][norepidx].values, axis=0, equal_var=False)
    waveformshape['phasex'] = phasex


    filehandler = open(
        "".join(('/Users/matsvanes/Data/YunzheData/Mats/Study2/emd/', f"group_cycle_metrics_imf{nimf + 1}.p")), "wb")
    pickle.dump({"waveformshape": waveformshape}, filehandler)
    filehandler.close()

plt.figure()
for i in range(38):
    plt.subplot(4, 10, i + 1)
    tmp1 = waveformshape['pa_if_noreplay_mean'][i, :, :] - waveformshape['pa_if_noreplay_mean'][i, :, :].mean(
        axis=0)  # demean individual's stats
    plt.plot(phasex, tmp1.mean(axis=1))
    plt.fill_between(phasex, tmp1.mean(axis=1) - (tmp1.std(axis=1) / np.sqrt(22)),
                     tmp1.mean(axis=1) + (tmp1.std(axis=1) / np.sqrt(22)), alpha=0.5)
    tmp2 = waveformshape['pa_if_replay_mean'][i, :, :] - waveformshape['pa_if_replay_mean'][i, :, :].mean(axis=0)
    plt.plot(phasex, tmp2.mean(axis=1))  # demean individual's stats
    plt.fill_between(phasex, tmp2.mean(axis=1) - (tmp2.std(axis=1) / np.sqrt(22)),
                     tmp2.mean(axis=1) + (tmp2.std(axis=1) / np.sqrt(22)), alpha=0.5)
    plt.ylim((-0.5, 0.5))

tvals, pvals = stats.ttest_rel(waveformshape['pa_r_mean'],
                               np.tile(np.nanmean(waveformshape['pa_r_mean'], axis=0), (38, 1, 1)), axis=2,
                               nan_policy='omit')
tvals, pvals = stats.ttest_rel(waveformshape['pa_if_1stlvl'], 0 * waveformshape['pa_if_1stlvl'], axis=2)

plt.figure()
for i in range(38):
    plt.subplot(4, 10, i + 1)
    plt.plot(phasex, tvals[i, :])
    try:
        msk = tvals[i, :].copy()
        msk[np.where(pvals[i, :] > 0.05 / 38)] = np.nan
        plt.plot(phasex, msk)
    except:
        print('')

    plt.figure()
    plt.subplot(131)
    plt.plot(noreplay['Comp1'], noreplay['Comp2'], '.')
    plt.plot(replay['Comp1'], replay['Comp2'], 'o')
    plt.legend(['No Replay', 'Replay'])
    plt.xlabel('Shape comp 1')
    plt.ylabel('Shape comp 2')

    plt.subplot(132)
    plt.plot(noreplay['Comp1'], noreplay['Comp3'], '.')
    plt.plot(replay['Comp1'], replay['Comp3'], 'o')
    plt.legend(['No Replay', 'Replay'])
    plt.xlabel('Shape comp 1')
    plt.ylabel('Shape comp 3')

    plt.subplot(133)
    plt.plot(noreplay['mean_if'], noreplay['max_amp'], '.')
    plt.plot(replay['mean_if'], replay['max_amp'], 'o')
    plt.legend(['No Replay', 'Replay'])
    plt.xlabel('Mean IF')
    plt.ylabel('Max Amp')
