from annotated_fasta import *
import os
from lib_prd import get_ac_list


def load_cross_disprot_caid_af(dis_prot_file, caid_file, verbose=False):
    dp_af = annotated_fasta_load(dis_prot_file)
    caid_af = annotated_fasta_load(caid_file)

    ac_list = list(caid_af['data'])
    for ac in ac_list:
        if ac in dp_af['data']:
            caid_af['data'][ac] = dp_af['data'][ac]
        else:
            lst = get_ac_list(seq=caid_af['data'][ac]['seq'])
            if verbose:
                print(ac, len(caid_af['data'][ac]['seq']), caid_af['data'][ac]['seq'])
                print(lst, '\n')
            del caid_af['data'][ac]
    caid_af['metadata']['tags'] = dp_af['metadata']['tags']
    return caid_af


def annotated_fasta_add_scores(af):
    for ac in af['data']:
        if 'scores' not in af['data'][ac]:
            af['data'][ac]['scores'] = {}


def annotated_fasta_load_prd_caid_scores(af, sc_caid_file, prd, del_no_scores=True):
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
        if 'scores' not in af['data'][ac]:
            if del_no_scores:
                del af['data'][ac]
            continue
        if prd in af['data'][ac]['scores']:
            if len(af['data'][ac]['scores'][prd]) != len(af['data'][ac]['seq']):
                print(prd, ac, 'Bad scores')
                del af['data'][ac]['scores'][prd]



def annotated_fasta_load_caid_scores(af, scores_path, prd_list):
    path_files = os.listdir(scores_path)
    annotated_fasta_add_scores(af)
    for prd_f in prd_list:
        if prd_f not in path_files:
            print(f"{prd_f} not found")
            continue
        annotated_fasta_load_prd_caid_scores(af, sc_caid_file=f'{scores_path}{prd_f}', prd=prd_f)

