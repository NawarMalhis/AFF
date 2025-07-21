from annotated_fasta import *


def aff_mdb_fasta_to_af(mdb_file_list: list=None, tag_list: list=None, quality_list: list=None):
    if mdb_file_list is None:
        print(f"Bad MobiDB fasta files:\t{mdb_file_list}", flush=True)
        return None
    if tag_list is None:
        tag_list = ['disorder', 'lip', 'binding_mode_disorder_to_disorder']
    if quality_list is None:
        quality_list = ['curated', 'derived']
    af = annotated_fasta(data_name="MobiDB-LIP")
    for fl in mdb_file_list:
        with open(fl, 'r') as fin:
            # ac = ''
            seq_next = False
            seq = ''
            tag = ''
            ac_seq = ''
            for line in fin:
                line = line.strip()
                if len(line) == 0:
                    continue
                if line[0] == '>':  # ==================================================
                    lst = line[1:].split('|')
                    if len(lst) >= 2:
                        if lst[1] == 'sequence':
                            if lst[0] not in af['data']:
                                ac_seq = lst[0]
                                af['data'][ac_seq] = {'seq': '', 'tags': {}, 'databases': {}, 'scores': {}}
                                seq_next = True
                            else:
                                seq_next = False
                            continue
                        else:
                            seq_next = False
                            t_lst = lst[1].split('-')
                            if t_lst[0] in quality_list and t_lst[1] in tag_list:
                                tag = lst[1]
                                continue
                            else:
                                tag = ''
                else:
                    if seq_next:
                        af['data'][ac_seq]['seq'] = line
                        continue
                    if len(tag) > 0 and len(ac_seq) > 0:
                        af['data'][ac_seq]['tags'][tag] = line
                        if tag not in af['metadata']['tags_dict']:
                            af['metadata']['tags_dict'][tag] = 'MobiDB tag'
                            af['metadata']['tags_list'].append(tag)
    for ac in af['data']:
        for tg in af['metadata']['tags_list']:
            if tg not in af['data'][ac]['tags']:
                af['data'][ac]['tags'][tg] = '-' * len(af['data'][ac]['seq'])

    af['metadata']['tags_list'].sort()
    # print(af['metadata']['tags_list'], flush=True)
    return af