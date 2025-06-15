from annotated_fasta import *
import numpy as np
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve, RocCurveDisplay
import matplotlib.pyplot as plt


colors_list = ['red', 'blue', 'forestgreen', 'limegreen', 'brown', 'pink', 'tomato', 'gray', 'darkred', 'darkgreen',
               'orange', 'yellow', 'darkseagreen', 'darkviolet', 'wheat', 'tan', 'cadetblue', 'deepskyblue', 'sienna']


# this add 'index' for what to delete and 'Y'
def aff_mask_ac(af1, tag):
    msk = list(af1[tag].replace('1', '0').replace('-', '1'))
    msk2 = [eval(i) for i in msk]  # list of ints
    af1['index'] = np.nonzero(msk2)  # Return the indices of the elements that are non-zero (needs not to be included).
    af1['Y'] = [eval(i) for i in np.delete(list(af1[tag]), af1['index'])] # replace eval with int


def aff_get_yx_dict(af, tag):
    yx_dict = {}
    prd_list = []
    for ac in af['data']:
        aff_mask_ac(af['data'][ac], tag)
        if 'scores' in af['data'][ac]:
            for prd in af['data'][ac]['scores']:
                if prd not in prd_list:
                    prd_list.append(prd)
    for prd in prd_list:
        # if prd not in yx_dict:
        #     yx_dict[prd] = {}  # 'yy': [], 'sc': [], 'miss_ac': []
        yy = np.array([], dtype='int32')
        sc = np.array([], dtype='float32')
        miss_ac = []
        for ac in af['data']:
            if prd in af['data'][ac]['scores']:
                yy = np.concatenate((yy, af['data'][ac]['Y']))
                sc = np.concatenate((sc, np.delete(af['data'][ac]['scores'][prd], af['data'][ac]['index'])))
            else:
                miss_ac.append(ac)
        yx_dict[prd] = {'yy': yy, 'sc': sc, 'miss_ac': miss_ac}
    return yx_dict


def aff_roc_figure(prd_use, af, tag, title=None, display=True, color_dict=None, min_auc=0.5,
                   sort_auc=True, out_file=None, auc_file=None):
    if title is None:
        title = tag
    auc_out = None
    yx_dict = aff_get_yx_dict(af, tag)
    print(f"af size:\t{len(af['data'])}")
    auc_dict = {}
    plt.rcParams.update({'font.size': 18})
    plt.rcParams['figure.figsize'] = [9.5, 9]
    for prd in prd_use:
        auc_dict[prd] = roc_auc_score(yx_dict[prd]['yy'], yx_dict[prd]['sc'])
    if sort_auc:
        auc_dict = dict(sorted(auc_dict.items(), key=lambda item: item[1], reverse=True))
    if display or out_file is not None:
        if auc_file:
            auc_out = open(auc_file, 'w')
        print("Predictor\tAUC\tmissing_seq\ttotal_AAs", file=auc_out)
        gray_lbl = 'General Protein Binding Tools'
        cc = 0
        for prd in auc_dict:
            print(f"{prd}\t{auc_dict[prd]:1.4f}\t{len(yx_dict[prd]['miss_ac'])}\t{len(yx_dict[prd]['yy'])}",
                  file=auc_out)
            if auc_dict[prd] < min_auc:
                continue
            fpr, tpr, _ = roc_curve(yx_dict[prd]['yy'], yx_dict[prd]['sc'])
            ls = 'solid'
            lbl = f"{prd} ({auc_dict[prd]:1.3f})"
            cc += 1
            clr = 'gray'
            if color_dict is not None:
                if prd in color_dict:
                    clr = color_dict[prd]['color']
                    ls = color_dict[prd]['line']
            else:
                clr = None
            if clr == 'gray':
                plt.plot(fpr, tpr, lw=1, label=gray_lbl, color=clr, linestyle=ls)
                gray_lbl = None
            else:
                plt.plot(fpr, tpr, lw=2, label=lbl, color=clr, linestyle=ls)  # colors_list[cc]

        plt.plot([0, 1], [0, 1], color="navy", lw=1, linestyle="--", label=f"Naive (0.500)")
        plt.legend(loc="lower right", fontsize=15)
        plt.title(title)
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        if out_file is not None:
            plt.savefig(out_file, dpi=300)
        if display:
            plt.show()
    return auc_dict


def aff_filter_for_success(af, tag, cut=5):
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        cnt_0 = af['data'][ac][tag].count('0')
        cnt_1 = af['data'][ac][tag].count('1')
        if cnt_0 < cut or cnt_1 < cut:
            del af['data'][ac]


def aff_success(af, tag, prd_list):
    succ_dict = {}
    for prd in prd_list:
        succ_cnt = 0.0
        for ac in af['data']:
            scores_c01 = [0, 0]
            cnt = [0, 0]
            for ii in range(len(af['data'][ac]['seq'])):
                tg = af['data'][ac][tag][ii]
                if tg in ['0', '1']:
                    scores_c01[int(tg)] += af['data'][ac]['scores'][prd][ii]
                    cnt[int(tg)] += 1
            succ_cnt += int(float(scores_c01[1]) / cnt[1] > float(scores_c01[0]) / cnt[0])
        succ_dict[prd] = float(succ_cnt) / len(af['data'])

    return succ_dict


# def get_high(af, tag):
#     if tag != 'N_bind':
#         return
#     n_max = 0.6
#     for ac in af['data']:
#         sz = len(af['data'][ac]['seq'])
#         for ii in range(sz):
#             if af['data'][ac]['IPA_nucleic'][ii] > n_max and af['data'][ac]['N_bind'][ii] == '1':
#                 n_max = af['data'][ac]['IPA_nucleic'][ii]
#     print(f"n_max:\t{n_max}")  # 0.87638
#     prot_nucleotide = {}
#     for ac in af['data']:
#         sz = len(af['data'][ac]['seq'])
#         for ii in range(sz):
#             if af['data'][ac]['IPA_nucleic'][ii] > n_max and af['data'][ac]['N_bind'][ii] == '0':
#                 if ac not in prot_nucleotide:
#                     prot_nucleotide[ac] = 0
#                 prot_nucleotide[ac] += 1
#                 print(ac, af['data'][ac]['seq'][ii], ii, af['data'][ac]['IPA_nucleic'][ii], af['data'][ac]['N_bind'][ii])
#
#     for ac in prot_nucleotide:
#         print(ac, prot_nucleotide[ac])


def aff_get_aps(yy, sc):
    lst = []
    precision = np.zeros(len(yy), dtype='float32')
    for y, s in zip(yy, sc):
        lst.append({'Y': y, 'S': s})
    lst_sorted = sorted(lst, key=lambda x: x['S'], reverse=True)
    accumulated = lst_sorted[0]['Y']
    precision[0] = accumulated
    for ii in range(1, len(yy)):
        accumulated = accumulated + lst_sorted[ii]['Y']
        precision[ii] = accumulated / (ii + 1)
    return np.average(precision)

def aff_precision_recall(y, sc, steps=400):
    pre_rec = {'precision': np.zeros(steps+1, dtype='float32'), 'recall': np.zeros(steps+1, dtype='float32')}
    mn = np.amin(sc)
    mx = np.amax(sc)
    cut_list = [mn]
    st = (mx - mn) / (steps-1)
    for i in range(steps):
        cut_list.append(cut_list[i] + st)
    for i, cut in enumerate(cut_list):
        yh = (sc > cut) * 1
        tn, fp, fn, tp = confusion_matrix(y, yh).ravel()
        if (tp + fp) > 0:
            pre_rec['precision'][i] = tp / (tp + fp)
        if (tp + fn) > 0:
            pre_rec['recall'][i] = tp / (tp + fn)

    return pre_rec


def aff_pr_figure(prd_use, af, tag, title=None, min_recall=0.05, display=True, color_dict=None, out_file=None):
    # APS: average precision score
    if title is None:
        title = tag
    yx_dict = aff_get_yx_dict(af, tag)

    max_precision = 0.0
    prd = prd_use[0]
    h_line = sum(yx_dict[prd]['yy']) / len(yx_dict[prd]['yy'])
    print(f"-----\t{h_line:1.4}")
    plt.rcParams.update({'font.size': 18})
    plt.rcParams['figure.figsize'] = [9.5, 9]
    plt.title(f"{title}", fontsize=18)
    gray_lbl = 'General Protein Binding Tools'

    for prd in prd_use:
        ls = 'solid'
        pre_rec = aff_precision_recall(yx_dict[prd]['yy'], yx_dict[prd]['sc'], steps=200)
        pre_rec['APS'] = aff_get_aps(yy=yx_dict[prd]['yy'], sc=yx_dict[prd]['sc'])
        print(f"{prd}:\t{pre_rec['APS']:0.4}")
        j0 = 0
        for j in range(len(pre_rec['recall'])):
            if max_precision < pre_rec['precision'][j]:
                max_precision = pre_rec['precision'][j]
            if pre_rec['recall'][j] < min_recall:
                j0 = j - 1
                break

        if (j0 + 1) < len(pre_rec['recall']):
            pre_rec['precision'][j0] = pre_rec['precision'][j0-1]
            pre_rec['recall'][j0] = min_recall
            j0 += 1
        pre_rec['precision'][0] = h_line
        pre_rec['recall'][0] = 0.999
        clr = None
        if color_dict is not None:
            if prd in color_dict:
                clr = color_dict[prd]['color']
                ls = color_dict[prd]['line']
        else:
            clr = None
        if clr == 'gray':
            plt.plot(pre_rec['recall'][:j0], pre_rec['precision'][:j0], color=clr, lw=1, label=gray_lbl)
            gray_lbl = None
        else:
            lbl = f'{prd} ({pre_rec["APS"]:.3})'
            plt.plot(pre_rec['recall'][:j0], pre_rec['precision'][:j0], linestyle=ls, lw=2, label=lbl)
    plt.plot([0, 1], [h_line, h_line], color="navy", lw=1, linestyle="--", label=f"Priors ({h_line:.3})")
    plt.ylim((0, max_precision + 0.05))
    plt.xlim((0, 1.0))
    plt.legend(loc="upper right", fontsize=16)
    plt.ylabel("Precision", fontsize=16)
    plt.xlabel("Recall", fontsize=16)
    if out_file is not None:
        plt.savefig(out_file, dpi=300)
    if display:
        plt.show()
