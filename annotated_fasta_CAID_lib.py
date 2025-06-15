from annotated_fasta import *
import os


def aff_load_prd_caid_scores(af, sc_caid_file, prd):
    with open(sc_caid_file, 'r') as fin:
        ac = ''
        for line in fin:
            line = line.strip()
            if len(line) < 2:
                continue
            if line[0] == '>':
                ac = line[1:]
                if ac in af['data']:
                    if 'scores' not in af['data'][ac]:
                        af['data'][ac]['scores'] = {}
                    af['data'][ac]['scores'][prd] = []
                else:
                    ac = ''
                continue
            if len(ac) < 2:
                continue
            sc = float(line.split()[2])
            af['data'][ac]['scores'][prd].append(sc)
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        # foe each predictor
        if prd in af['data'][ac]['scores']:
            if len(af['data'][ac]['scores'][prd]) != len(af['data'][ac]['seq']):
                # scores do not match this seq
                del af['data'][ac]['scores'][prd]


def aff_load_caid_scores_single_file(af, scores_path, prd_list, remove_missing_scores=False):
    path_files = os.listdir(scores_path)
    # prd_cnt = {}
    for prd in prd_list:
        if f"{prd}.caid" not in path_files:
            print(f"{prd} not found")
            continue
        # prd_cnt[prd] = 0
        aff_load_prd_caid_scores(af, sc_caid_file=f'{scores_path}{prd}.caid', prd=prd)

    used_prd_set = set()
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        for prd in af['data'][ac]['scores']:
            # prd_cnt[prd] += 1
            used_prd_set.add(prd)

    if remove_missing_scores:
        aff_remove_missing_scores(af)
    return used_prd_set


def aff_remove_missing_scores(af):
    used_prd_set = set()
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        for prd in af['data'][ac]['scores']:
            used_prd_set.add(prd)
    for ac in ac_list:
        for prd in used_prd_set:
            if prd not in af['data'][ac]['scores']:
                del af['data'][ac]
                break
