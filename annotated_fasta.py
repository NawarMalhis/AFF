import copy
from crc64iso.crc64iso import crc64
import requests


def annotated_fasta(data_name: str='Data has no Name', names_list: list=None, accession: str=None):
    if names_list is None:
        names_list = []
    if accession is not None:
        if accession not in names_list:
            names_list.append(accession)
    return {'data': {}, 'metadata': {'tags_dict': {}, 'names_list': names_list, 'counts': None,
                                     'accession': accession, 'data_name': data_name}}


def aff_load0(in_file: str, data_name: str='Data has no Name', accession: str=None):  # , _mark=None
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

            if line[0] == '>':
                ac_lst = line[1:].split('|')
                ac = ac_lst[0]
                af_sequences[ac] = {'seq': '', 'scores': {}}
                for extra in ac_lst[1:]:
                    ex_lst = extra.split('=')
                    if len(ex_lst) < 2:
                        continue
                    af_sequences[ac][ex_lst[0]] = ex_lst[1]
                    if ex_lst[0] not in name_tags:
                        name_tags.append(ex_lst[0])
                continue
            af_sz = len(af_sequences[ac])
            if af_sequences[ac]['seq'] == '':
                af_sequences[ac]['seq'] = line
                tg_i0 = af_sz  # len(af_sequences[ac])
            else:
                af_sequences[ac][tags[af_sz - tg_i0]] = line.replace('x', '-')
            continue
    af0 = {'data': af_sequences, 'metadata': {'tags': tags, 'name_tags': name_tags, 'statistics': None}}
    # -------------------------------------------------------------------------------------------------------
    # for ac in af0['data']:
    #     for ntg in af0['metadata']['name_tags']:
    #
    af = annotated_fasta(data_name=data_name)
    for tg in af0['metadata']['tags']:
        af['metadata']['tags_dict'][tg] = 'xxx'
    if accession is None:
        accession = 'Fasta'
    af['metadata']['accession'] = accession
    af['metadata']['names_list'].append(accession)
    for ntg in af0['metadata']['name_tags']:
        af['metadata']['names_list'].append(ntg.split(';'))
    for ac in af0['data']:
        af['data'][ac] = af0['data'][ac]
        af['data'][ac][accession] = [ac]
        for ntg in af0['metadata']['name_tags']:
            af['data'][ac][ntg] = af0['data'][ac][ntg].split(';')
    # print(len(af['data']), flush=True)
    # for ac in af['data']:
    #     print(ac, list(af['data'][ac].keys()))
    _gen_name_counts(af)
    return af


def aff_load2(in_file: str, data_name: str='Data has no name'):  # , _mark=None
    af_sequences = {}
    tags_dict = {}
    tags_list = []
    names_list = []
    _more_tags = False
    _id_counts = False
    accession = None
    with open(in_file, 'r') as fin:
        ac = ''
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                if 'Data Name:' in line:
                    data_name = line.split('\t')[1]
                if 'Amino acid sequence' in line:
                    _more_tags = True
                    continue
                if _more_tags:
                    _lst = line.split('\t')
                    if len(_lst) > 1:
                        if len(_lst) > 1:
                            info = _lst[2]
                        else:
                            info = 'Annotation'
                        tags_dict[_lst[1]] = info
                        tags_list.append(_lst[1])
                    if 'IDs Counts:' in line:
                        _more_tags = False
                        _id_counts = True
                    continue
                if _id_counts:
                    _lst = line.split()
                    if 'Unique#' in line:
                        continue
                    if 'Tags Counts:' in line:
                        _id_counts = False
                        continue
                    if len(_lst) < 4:
                        continue
                    if len(_lst) == 5:
                        accession = _lst[1]
                    names_list.append(_lst[1])
                continue
            if len(tags_dict) == 0:
                print(f"Error: No tags are found in {in_file}")
                return None

            if line[0] == '>':
                ac_lst = line[1:].split('|')
                ac = ac_lst[0]
                af_sequences[ac] = {'seq': '', 'scores': {}}
                for extra in ac_lst[1:]:
                    ex_lst = extra.split('=')
                    if len(ex_lst) != 2:
                        continue
                    lst2 = ex_lst[1].split(';')
                    af_sequences[ac][ex_lst[0]] = lst2
                    if ex_lst[0] not in names_list:
                        names_list.append(ex_lst[0])
                continue
            af_sz = len(af_sequences[ac])
            if af_sequences[ac]['seq'] == '':
                af_sequences[ac]['seq'] = line
                tg_i0 = af_sz
            else:
                af_sequences[ac][tags_list[af_sz - tg_i0]] = line.replace('x', '-')
            continue
    af = {'data': af_sequences, 'metadata': {'tags_dict': tags_dict, 'names_list': names_list, 'counts': None,
                                             'accession': accession, 'data_name': data_name}}
    for ac in af['data']:
        for ntg in af['metadata']['names_list']:
            if ntg not in af['data'][ac]:
                af['data'][ac][ntg] = ''
    _gen_name_counts(af)
    return af


def aff_load_simple(in_file: str, data_name: str='Data has no name', tag: str= 'ANN',
                    tag_description: str= 'Source annotation'):
    af_sequences = {}
    tags_dict = {tag: tag_description}
    names_list = ['Fasta']
    _more_tags = False
    _id_counts = False
    accession = 'Fasta'
    l_num = 0
    with open(in_file, 'r') as fin:
        ac = ''
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            if line[0] == '>':
                l_num = 1
                hd_lst = line[1:].split('|')
                ac = hd_lst[0]
                af_sequences[ac] = {accession: [ac]}
                for hdn in hd_lst:
                    two_parts = hdn.replace(':', '=').split('=')
                    if len(two_parts) == 2:
                        af_sequences[ac][two_parts[0]] = two_parts[1].split(';')
                        if two_parts[0] not in names_list:
                            names_list.append(two_parts[0])
                continue
            if l_num == 1:
                af_sequences[ac]['seq'] = line
                l_num = 2
                continue
            if l_num == 2:
                af_sequences[ac][tag] = line
                l_num = 0
                continue
            # else:
            #     print(f"ERROR\t{line}", flush=True)
    af = {'data': af_sequences, 'metadata': {'tags_dict': tags_dict, 'names_list': names_list, 'counts': None,
                                             'accession': accession, 'data_name': data_name}}
    for ac in af['data']:
        for ntg in af['metadata']['names_list']:
            if ntg not in af['data'][ac]:
                af['data'][ac][ntg] = ''
    _gen_name_counts(af)
    return af


def aff_load_fasta(in_file: str, data_name: str='Data has no name'):
    af = annotated_fasta(data_name=data_name)
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
                af['data'][ac] = {'seq': '', 'scores': {}}
            else:
                af['data'][ac]['seq'] = af['data'][ac]['seq'] + line
    return af

def _validate_header_extra(he: str=None):
    _ret = True
    if he is None:
        return True
    lst = he.split('\n')
    for line in lst:
        if len(line) == 0:
            return False
        if line[0] != '#':
            return False
    return True


def aff_save2(af, f_name: str, header_top: str =None, header_bottom: str =None):
    with open(f_name, 'w') as fout:
        print(f"# Data Name:\t{af['metadata']['data_name']}\n#", file=fout)
        if header_top:
            if _validate_header_extra(header_top):
                print(header_top, file=fout)
            else:
                print(f"BAD header_top:\t{header_top}", flush=True)
        print(f"# Sequences:\t{len(af['data']):,}", file=fout)
        print("#", file=fout)
        print("# Format:", file=fout)
        print("#\t>accession", file=fout, end='')
        for n_tg in af['metadata']['names_list']:
            print(f"|{n_tg}={n_tg}_ID", file=fout, end='')
        print(file=fout)
        print("#\tAmino acid sequence", file=fout)
        for tg in af['metadata']['tags_dict']:
            print(f"#\t{tg}\t{af['metadata']['tags_dict'][tg]}", file=fout)
        print("#", file=fout)
        aff_gen_counts(af)
        str_counts = _get_string_counts(af)
        print(str_counts, file=fout)
        if header_bottom:
            if _validate_header_extra(header_bottom):
                print('#', file=fout)
                print(header_bottom, file=fout)
            else:
                print(f"BAD header_bottom:\t{header_bottom}", flush=True)
        print('#', file=fout)
        for ac in af['data']:
            ac_o = ac
            for tg in af['data'][ac]:
                if tg in af['metadata']['names_list']:
                    ac_o = f"{ac_o}|{tg}={af['data'][ac][tg][0]}"
                    if len(af['data'][ac][tg]) > 1:
                        for xx in af['data'][ac][tg][1:]:
                            ac_o = f"{ac_o};{xx}"
            print(f">{ac_o}\n{af['data'][ac]['seq']}", file=fout)
            for tg in af['metadata']['tags_dict']:
                print(f"{af['data'][ac][tg]}", file=fout)
    return


def aff_save_fasta(af, f_name: str):
    with open(f_name, 'w') as fout:
        for ac in af['data']:
            ac_o = ac
            print(f">{ac_o}\n{af['data'][ac]['seq']}", file=fout)


def aff_merge_simple(af_to, af_from):
    for ac in af_from['data']:
        if ac not in af_to['data']:
            af_to['data'][ac] = af_from['data'][ac]


def merge_annotations(ann1: str, ann2: str):
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


# needs validation
def aff_remove_tags_list(af, tags_list_out: list):
    tags_dict = af['metadata']['tags_dict']
    for tg in af['metadata']['tags_dict']:
        if tg in tags_list_out:
            del tags_dict[tg]
    for ac in af['data']:
        for tg in af['metadata']['tags_dict']:
            if tg in tags_list_out:
                del af['data'][ac][tg]
    af['metadata']['tags_dict'] = tags_dict


def aff_rename_tag(af, old_tag: str, new_tag: str, new_info: str):
    tag_used = False
    for tg in range(len(af['metadata']['tags_dict'])):
        if tg == old_tag:
            af['metadata']['tags_dict'][new_tag] = new_info
            del af['metadata']['tags_dict'][tg]
            tag_used = True
            break
    if tag_used:
        for ac in af['data']:
            af['data'][ac][new_tag] = af['data'][ac][old_tag]
            del af['data'][ac][old_tag]


def aff_remove_missing_scores(af):
    used_scores_set = set()
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        for prd_sc in af['data'][ac]['scores']:
            used_scores_set.add(prd_sc)
    for ac in ac_list:
        for prd_sc in used_scores_set:
            if prd_sc not in af['data'][ac]['scores']:
                del af['data'][ac]
                break


def aff_remove_no_info_tag(af, tag: str):
    if tag not in af['metadata']['tags_dict']:
        return
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        if af['data'][ac][tag].count('-') == len(af['data'][ac][tag]):
            del af['data'][ac]


def aff_remove_no_class_tag(af, tag: str, cl: str):
    if tag not in af['metadata']['tags_dict']:
        return
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        if cl not in af['data'][ac][tag]:
            del af['data'][ac]


def aff_remove_no_class_any(af, cl: str):
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        rmv = True
        for tg in af['metadata']['tags_dict']:
            if cl in af['data'][ac][tg]:
                rmv = False
                break
        if rmv:
            del af['data'][ac]


def aff_add_tag(af, tag: str):
    if tag not in af['metadata']['tags_dict']:
        af['metadata']['tags_dict'].append(tag)
        for ac in af['data']:
            af['data'][ac][tag] = '-' * len(af['data'][ac]['seq'])
    else:
        print(f"{tag} exist in af", flush=True)


def _get_string_counts(af):
    _msg = ''
    if af['metadata']['counts'] is None:
        return _msg
    if len(af['metadata']['names_list']) > 0:
        _msg = _msg + f"# IDs Counts:\n#\tID \tAll#\tUnique#"
        for ntg in af['metadata']['names_list']:
            _msg = _msg + f"\n#\t{ntg}\t{af['metadata']['counts']['names_dict'][ntg]['total']:,}"
            _msg = _msg + f"\t{af['metadata']['counts']['names_dict'][ntg]['unique']:,}"
            if ntg == af['metadata']['accession']:
                _msg = _msg + "\tAC"
        _msg = _msg + "\n#\n"
    _msg = _msg + "# Tags Counts:\n#\ttag\tSeq#\tSeg#\t'0'\t'1'\t'-'"
    for tg in af['metadata']['tags_dict']:
        _msg = _msg + f"\n#\t{tg}"
        for cc in ['seq', 'seg', '0', '1', '-']:  # af['metadata']['counts']['tags_dict'][tg]:
            _msg = _msg + f"\t{af['metadata']['counts']['tags_dict'][tg][cc]:,}"
    return _msg


def aff_gen_counts(af):
    _gen_tag_counts(af)
    _gen_name_counts(af)


def _gen_tag_counts(af):
    if af['metadata']['counts'] is None:
        af['metadata']['counts'] = {}
    af['metadata']['counts']['tags_dict'] = {}
    for tg in af['metadata']['tags_dict']:
        af['metadata']['counts']['tags_dict'][tg] = {'seq': 0, 'seg': 0, '0': 0, '1': 0, '-': 0}
        for ac in af['data']:
            mask = str(af['data'][ac][tg])
            mask = mask.replace('-', '0')
            cnt = len([xx for xx in mask.split('0') if xx])
            if cnt > 0:
                af['metadata']['counts']['tags_dict'][tg]['seq'] += 1
                af['metadata']['counts']['tags_dict'][tg]['seg'] += cnt
            for cc in ['0', '1', '-']:
                af['metadata']['counts']['tags_dict'][tg][cc] += af['data'][ac][tg].count(cc)


def _gen_name_counts(af):
    if af['metadata']['counts'] is None:
        af['metadata']['counts'] = {}
    af['metadata']['counts']['names_dict'] = {}
    ntg_set_dict = {}
    for ntg in af['metadata']['names_list']:
        af['metadata']['counts']['names_dict'][ntg] = {'total': 0, 'unique': 0}
        ntg_set_dict[ntg] = set()
        for ac in af['data']:
            if ntg in af['data'][ac]:
                # print(ac, af['data'][ac][ntg], flush=True)
                if len(af['data'][ac][ntg]) > 0:
                    af['metadata']['counts']['names_dict'][ntg]['total'] += 1
                    ntg_set_dict[ntg].add(af['data'][ac][ntg][0])
    for ntg in af['metadata']['names_list']:
        af['metadata']['counts']['names_dict'][ntg]['unique'] = len(ntg_set_dict[ntg])


def _get_url(url, **kwargs):
    response = requests.get(url, **kwargs)
    if not response.ok:
        print(response.text)
        return None
        # response.raise_for_status()
        # sys.exit()
    return response


def get_ac_list(seq):
    print(seq)
    ac_list = []
    checksum = crc64(seq)
    r = _get_url(f"https://rest.uniprot.org/uniparc/search?query=checksum: {checksum}")
    if r is not None:
        data = r.json()
        # for xx in data["results"][0]['uniParcCrossReferences']:
            # if 'UniProtKB' in xx['database']:
            #     ac_list.append(xx['id'])
        for ac in data["results"][0]['uniProtKBAccessions']:
            response = requests.get(f'https://rest.uniprot.org/uniprotkb/{ac}.fasta')
            if response.ok:
                r_seq = ''
                for ss in response.text.split('\n')[1:]:
                    if len(ss) > 0:
                        r_seq = r_seq + ss
                if r_seq == seq:
                    ac_list.append(ac)
    return ac_list



