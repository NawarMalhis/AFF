import copy
from annotated_fasta import *
import numpy as np
from sklearn.metrics import average_precision_score, confusion_matrix, roc_auc_score, roc_curve, precision_recall_curve
import matplotlib.pyplot as plt


colors_list = ['red', 'blue', 'forestgreen', 'limegreen', 'brown', 'pink', 'tomato', 'gray', 'darkred', 'darkgreen',
               'orange', 'yellow', 'darkseagreen', 'darkviolet', 'wheat', 'tan', 'cadetblue', 'deepskyblue', 'sienna']


# this add 'index' for what to delete and 'Y'
def _mask_ac(af1, tag):
    msk = list(af1[tag].replace('1', '0').replace('-', '1'))
    msk2 = [eval(i) for i in msk]  # list of ints
    af1['index'] = np.nonzero(msk2)  # Return the indices of the elements that are non-zero (needs not to be included).
    af1['Y'] = [eval(i) for i in np.delete(list(af1[tag]), af1['index'])] # replace eval with int


def _get_yx_dict(af, tag):
    yx_dict = {}
    prd_list = []
    for ac in af['data']:
        _mask_ac(af['data'][ac], tag)
        if 'scores' in af['data'][ac]:
            for prd in af['data'][ac]['scores']:
                if prd not in prd_list:
                    prd_list.append(prd)
    for prd in prd_list:
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


def _fill_line_format_dict(prd_list, line_format_dict):
    lf_dict = {}
    if line_format_dict is None:
        line_format_dict = {}
    if len(line_format_dict) == 0:
        for prd in prd_list:
            lf_dict[prd] = {'weight': None, 'style': None, 'color': None, 'formated': True}
    else:
        for prd in prd_list:
            if prd not in line_format_dict:
                lf_dict[prd] = {'weight': 1, 'style': 'solid', 'color': 'gray', 'formated': False}
                continue
            lf_dict[prd] = {'formated': True}
            for ff in ['weight', 'style', 'color']:
                if ff in line_format_dict[prd]:
                    lf_dict[prd][ff] = line_format_dict[prd][ff]
                else:
                    lf_dict[prd][ff] = None
    return lf_dict


def _filter_for_success(af, tag, cut=5):
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        cnt_0 = af['data'][ac][tag].count('0')
        cnt_1 = af['data'][ac][tag].count('1')
        if cnt_0 < cut or cnt_1 < cut:
            del af['data'][ac]


def aff_roc(af, tag, prd_list, title=None, min_auc=0.5, display=True, line_format_dict=None,
            figure_file=None, auc_file=None, legend_font_size=12, gray_label='Others', sort_auc=True):  #
    if title is None:
        title = tag
    plotted_list = []

    lf_dict = _fill_line_format_dict(prd_list=prd_list, line_format_dict=line_format_dict)
    auc_out = None
    yx_dict = _get_yx_dict(af, tag)
    auc_dict = {}
    plt.rcParams.update({'font.size': 18})
    plt.rcParams['figure.figsize'] = [10, 9]
    for prd in prd_list:
        auc_dict[prd] = roc_auc_score(yx_dict[prd]['yy'], yx_dict[prd]['sc'])
    if sort_auc:
        auc_dict = dict(sorted(auc_dict.items(), key=lambda item: item[1], reverse=True))

    if display or figure_file is not None:
        if auc_file:
            auc_out = open(auc_file, 'w')
        print("Predictor\tAUC\tmissing_seq\ttotal_AAs", file=auc_out)
        for formated in [True, False]:
            for prd in auc_dict:
                if lf_dict[prd]['formated'] != formated:
                    continue
                print(f"{prd}\t{auc_dict[prd]:1.4f}\t{len(yx_dict[prd]['miss_ac'])}\t{len(yx_dict[prd]['yy'])}",
                      file=auc_out)
                if auc_dict[prd] < min_auc:
                    continue
                plotted_list.append(prd)
                fpr, tpr, _ = roc_curve(yx_dict[prd]['yy'], yx_dict[prd]['sc'])
                if not lf_dict[prd]['formated']:
                    lbl = gray_label
                    gray_label = None
                else:
                    lbl = f"{prd} ({auc_dict[prd]:1.3f})"
                clr = lf_dict[prd]['color']
                ls = lf_dict[prd]['style']
                lw = lf_dict[prd]['weight']

                plt.plot(fpr, tpr, lw=lw, label=lbl, color=clr, linestyle=ls)

        plt.plot([0, 1], [0, 1], color="navy", lw=1, linestyle="--", label=f"Naive (0.500)")
        plt.legend(loc="lower right", fontsize=legend_font_size)
        plt.title(title)
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        if figure_file is not None:
            plt.savefig(figure_file, dpi=300)
        if display:
            plt.show()
    return auc_dict, plotted_list


def aff_precision_recall(af, tag, prd_list, title=None, min_recall=0.05, display=True, line_format_dict=None,
                         figure_file=None, aps_file=None, legend_font_size=12, gray_label='Others', plotted_list=None):
    if title is None:
        title = tag
    lf_dict = _fill_line_format_dict(prd_list, line_format_dict)
    yx_dict = _get_yx_dict(af, tag)
    aps_dict = {}
    plt.rcParams.update({'font.size': 18})
    plt.rcParams['figure.figsize'] = [10, 9]
    aps_out = None
    if aps_file:
        aps_out = open(aps_file, 'w')
    print("Predictor\tAPS", file=aps_out)
    for prd in prd_list:
        aps_dict[prd] = average_precision_score(yx_dict[prd]['yy'], yx_dict[prd]['sc'])
        print(f"{prd}\t{aps_dict[prd]:1.3}", file=aps_out)

    max_precision = 0.0
    data_priors = sum(yx_dict[prd_list[0]]['yy']) / len(yx_dict[prd_list[0]]['yy'])
    if display or figure_file is not None:
        for formated in [True, False]:
            for prd in plotted_list:
                if lf_dict[prd]['formated'] != formated:
                    continue
                precision, recall, thresholds = precision_recall_curve(yx_dict[prd]['yy'], yx_dict[prd]['sc'])
                # if prd in line_format_dict:
                clr = lf_dict[prd]['color']
                ls = lf_dict[prd]['style']
                lw = lf_dict[prd]['weight']
                if not lf_dict[prd]['formated']:
                    lbl = gray_label
                    gray_label = None
                else:
                    lbl = f"{prd} ({aps_dict[prd]:1.3f})"

                j0 = 0
                for j in range(len(recall)):
                    if max_precision < precision[j]:
                        max_precision = precision[j]
                    if recall[j] < min_recall:
                        j0 = j - 1
                        break

                if (j0 + 1) < len(recall):
                    precision[j0] = precision[j0-1]
                    recall[j0] = min_recall
                    j0 += 1
                precision[0] = data_priors
                recall[0] = 0.999
                plt.plot(recall[:j0], precision[:j0], color=clr, lw=lw, linestyle=ls, label=lbl)
        plt.plot([0, 1], [data_priors, data_priors], color="navy", lw=1, linestyle="--",
                 label=f"Priors ({data_priors:.3})")
        plt.ylim((0, max_precision + 0.05))
        plt.xlim((0, 1.0))
        plt.legend(loc="upper right", fontsize=legend_font_size)
        plt.ylabel("Precision", fontsize=16)
        plt.xlabel("Recall", fontsize=16)
        plt.title(f"{title}", fontsize=18)
        if figure_file is not None:
            plt.savefig(figure_file, dpi=300)
        if display:
            plt.show()
    return aps_dict


def aff_success_rate(af, tag, prd_list, success_rate_file=None, success_data_file=None):
    af_copy = copy.deepcopy(af)
    aff_remove_missing_scores(af_copy)
    _filter_for_success(af_copy, tag=tag, cut=5)
    if success_data_file is not None:
        aff_save(af=af_copy, f_name=success_data_file)
    success_rate_dict = {}
    for prd in prd_list:
        success_cnt = 0.0
        for ac in af_copy['data']:
            scores_c01 = [0, 0]
            cnt = [0, 0]
            for ii in range(len(af_copy['data'][ac]['seq'])):
                tg = af_copy['data'][ac][tag][ii]
                if tg in ['0', '1']:
                    scores_c01[int(tg)] += af_copy['data'][ac]['scores'][prd][ii]
                    cnt[int(tg)] += 1
            success_cnt += int(float(scores_c01[1]) / cnt[1] > float(scores_c01[0]) / cnt[0])
        success_rate_dict[prd] = float(success_cnt) / len(af_copy['data'])
    sr_out = None
    if success_rate_file:
        sr_out = open(success_rate_file, 'w')
    print("Predictor\tSuccess Rate", file=sr_out)
    for prd in prd_list:
        print(f"{prd}\t{success_rate_dict[prd]:1.4f}", file=sr_out)
    return success_rate_dict

