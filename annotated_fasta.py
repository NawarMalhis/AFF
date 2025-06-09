import copy


def annotated_fasta(tags=None, name_tags=None):
    if tags is None:
        tags = []
    if name_tags is None:
        name_tags = []
    return {'data': {}, 'metadata': {'tags': tags, 'name_tags': name_tags, 'statistics': None}}


def annotated_fasta_load(in_file: str):  # , _mark=None
    af_sequences = {}
    tags = []
    name_tags = []
    _more_tags = True
    with open(in_file, 'r') as fin:
        ac = ''
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                if _more_tags:
                    _lst = line.split()
                    if len(_lst) > 1:
                        if _lst[1] == 'TAG':
                            tags.append(_lst[2])
                continue
            _more_tags = False
            if len(tags) == 0:
                tags.append('mask')
            if line[0] == '>':
                ac_lst = line[1:].split('|')
                ac = ac_lst[0]
                af_sequences[ac] = {'seq': ''}
                for extra in ac_lst[1:]:
                    ex_lst = extra.split('=')
                    if len(ex_lst) < 2:
                        continue
                    af_sequences[ac][ex_lst[0]] = ex_lst[1]
                    if ex_lst[0] not in name_tags:
                        name_tags.append(ex_lst[0])
                #     if len(ac_lst) == 2:
                #         ox_lst = ac_lst[1].split('=')
                #         if len(ox_lst) > 1:
                #             if ox_lst[0] == 'OX':
                #                 ox = ox_lst[1]
                # af_sequences[ac] = {'mark': _mark, 'OX': ox}
                continue
            af_sz = len(af_sequences[ac])
            if af_sequences[ac]['seq'] == '':
                af_sequences[ac]['seq'] = line
                tg_i0 = af_sz  # len(af_sequences[ac])
            else:
                af_sequences[ac][tags[af_sz - tg_i0]] = line.replace('x', '-')
            continue
    af = {'data': af_sequences, 'metadata': {'tags': tags, 'name_tags': name_tags, 'statistics': None}}
    # print(af['metadata']['name_tags'], flush=True)
    return af


def annotated_fasta_gen_statistics(af):
    # counts = {'seq': 0, 'seg': 0}
    if len(af['data']) == 0:
        af['metadata']['statistics'] = None
        return
    af['metadata']['statistics'] = {}
    for tg in af['metadata']['tags']:
        af['metadata']['statistics'][tg] = {}
        af['metadata']['statistics'][tg]['seq'] = 0
        af['metadata']['statistics'][tg]['seg'] = 0
        for ac in af['data']:
            mask = str(af['data'][ac][tg])
            mask = mask.replace('-', '0')
            cnt = len([xx for xx in mask.split('0') if xx])
            if cnt > 0:
                af['metadata']['statistics'][tg]['seq'] += 1
                af['metadata']['statistics'][tg]['seg'] += cnt
        for cc in ['0', '1', '-']:
            af['metadata']['statistics'][tg][cc] = 0
            for ac in af['data']:
                for i in range(len(af['data'][ac][tg])):
                    if af['data'][ac][tg][i] == cc:
                        af['metadata']['statistics'][tg][cc] += 1


def get_string_stat(af):
    if not af['metadata']['statistics']:
        annotated_fasta_gen_statistics(af)
    _msg = "# Statistics:\n#\t---\ttag\tSeq#\tSeg#\t'0'\t'1'\t'-'"
    for tg in af['metadata']['tags']:
        _msg = _msg + f"\n#\tTAG\t{tg}"
        # print(af['metadata']['statistics'].keys(), flush=True)
        for cc in af['metadata']['statistics'][tg]:
            _msg = _msg + f"\t{af['metadata']['statistics'][tg][cc]:,}"
    return _msg


def annotated_fasta_save_fasta(af, f_name):
    with open(f_name, 'w') as fout:
        for ac in af['data']:
            ac_o = ac
            # for tg in af['data'][ac]:
            #     if tg not in af['metadata']['tags']:
            #         ac_o = f"{ac_o}|{tg}={af['data'][ac][tg]}"
            print(f">{ac_o}\n{af['data'][ac]['seq']}", file=fout)


def annotated_fasta_save(af, f_name, data_name=None,  header_top=None, header_bottom=None):
    with open(f_name, 'w') as fout:
        if data_name:
            print(f"# dataset: {data_name}\n#", file=fout)
        if header_top:
            print(header_top, file=fout)
        print(f"# Sequences:\t{len(af['data']):,}", file=fout)
        print("#", file=fout)
        print("# Format:", file=fout)
        print("#\t>accession", file=fout)
        print("#\tAmino acid sequence", file=fout)
        for tg in af['metadata']['tags']:
            print(f"#\t{tg} annotation", file=fout)
        print("#", file=fout)
        print(get_string_stat(af), file=fout)
        if header_bottom:
            print('#', file=fout)
            print(header_bottom, file=fout)
        print('#', file=fout)
        for ac in af['data']:
            ac_o = ac
            for tg in af['data'][ac]:
                if tg in af['metadata']['name_tags']:
                    ac_o = f"{ac_o}|{tg}="
                    for nm in af['data'][ac][tg]:
                        if ac_o[-1] != '=':
                            ac_o = f"{ac_o};"
                        ac_o = f"{ac_o}{nm}"
            print(f">{ac_o}\n{af['data'][ac]['seq']}", file=fout)
            for tg in af['metadata']['tags']:
                print(f"{af['data'][ac][tg]}", file=fout)
    return


def annotated_fasta_remove_tags_list(af, tags_list_out):
    tag_list = []
    for tg in af['metadata']['tags']:
        if tg not in tags_list_out:
            tag_list.append(tg)
    for ac in af['data']:
        for tg in af['metadata']['tags']:
            if tg in tags_list_out:
                del af['data'][ac][tg]
    af['metadata']['tags'] = tag_list


def annotated_fasta_rename_tag(af, old_tag, new_tag):
    for ii in range(len(af['metadata']['tags'])):
        if af['metadata']['tags'][ii] == old_tag:
            af['metadata']['tags'][ii] = new_tag
            break
    for ac in af['data']:
        af['data'][ac][new_tag] = af['data'][ac][old_tag]  # copy.deepcopy
        del af['data'][ac][old_tag]


def annotated_fasta_merge_simple(af_to, af_from):
    for ac in af_from['data']:
        if ac not in af_to['data']:
            af_to['data'][ac] = af_from['data'][ac]


def annotated_fasta_merge2(af1, af2):
    # print(af2['metadata']['tags'])
    merged2 = {'data': copy.deepcopy(af1['data']), 'metadata': {'tags': af1['metadata']['tags'], 'statistics': None}}
    for ac in af2['data']:
        if ac not in merged2['data']:
            print(f"AC {ac} not in merged2", flush=True)
            merged2['data'][ac] = copy.deepcopy(af2['data'][ac])
        else:
            if len(af2['data'][ac]['seq']) != len(merged2['data'][ac]['seq']):
                print(f"{ac} DELETED: seq size is different among the two input af data", flush=True)
                del merged2['data'][ac]
                continue
            if af2['data'][ac]['seq'] != merged2['data'][ac]['seq']:
                w_msg = f"{ac} sequences are not identical among the two input af data, same size."
                print(f"WARNING: {w_msg}", flush=True)
                if 'warnings' not in merged2:
                    merged2['warnings'] = {}
                merged2['warnings'][ac] = {'message': w_msg, 'alternative': af2['data'][ac]['seq']}
            for tg in af2['metadata']['tags']:
                if tg == 'seq':
                    continue
                if tg not in merged2['data'][ac].keys():
                    merged2['data'][ac][tg] = af2['data'][ac][tg]
                else:
                    m_tag = list(merged2['data'][ac][tg])
                    for i in range(len(af2['data'][ac][tg])):
                        if af2['data'][ac][tg][i] == '1':
                            m_tag[i] = '1'
                        elif af2['data'][ac][tg][i] == '0':
                            if m_tag[i] == '-':
                                m_tag[i] = '0'
                    merged2['data'][ac][tg] = ''.join(m_tag)
    return merged2


def annotated_fasta_merge_list(af_lst):
    if len(af_lst) == 0:
        return None
    elif len(af_lst) == 1:
        return af_lst[0]
    merged = copy.deepcopy(af_lst[0])
    for af in af_lst[1:]:
        merged = annotated_fasta_merge2(merged, af)
    return merged


def annotated_fasta_remove_no_1_tag(af, tag):
    if tag not in af['metadata']['tags']:
        return
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        if '1' not in af['data'][ac][tag]:
            # print(ac)
            del af['data'][ac]


def annotated_fasta_remove_no_info_tag(af, tag):
    if tag not in af['metadata']['tags']:
        return
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        if af['data'][ac][tag].count('-') == len(af['data'][ac][tag]):
            del af['data'][ac]


def annotated_fasta_remove_no_1_all(af):
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        rmv = True
        for tg in af['metadata']['tags']:
            if '1' in af['data'][ac][tg]:
                rmv = False
                break
        if rmv:
            del af['data'][ac]


def annotated_fasta_load_fasta(in_file):
    af = annotated_fasta()
    with open(in_file, 'r') as fin:
        ac = ''
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            if line[0] == '>':
                ac = line[1:]
                af['data'][ac] = {'seq': ''}
            else:
                af['data'][ac]['seq'] = af['data'][ac]['seq'] + line
    return af


def annotated_fasta_merge_annotations(ann1, ann2):
    if len(ann1) != len(ann2):
        return None
    sz = len(ann1)
    if sz <= 3:
        return None
    lst = ['-'] * sz
    for ii in range(sz):
        if ann1[ii] == '1' or ann2[ii] == '1':
            lst[ii] = '1'
        elif ann1[ii] == '0' or ann2[ii] == '0':
            lst[ii] = '0'
    return ''.join(lst)
