from annotated_fasta import *
import json

release = '2024_12'
dp_path = '/home/nmalhis/data/DisProt/'
_in_path = f'{dp_path}{release}/'
_out_path = f"{_in_path}af/"
# tags_list = ['IDR', 'DtO', 'Linker', 'Binding', 'P_bind', 'N_bind', 'I_bind', 'L_bind', 'SM_bind']
tags_names_dict = {'IDR': 'Protein disordered region', 'DtO': 'Disordered to ordered transition',
                   'Linker': 'Linker regions', 'binding': 'IDR binding in general',
                   'binding_protein': 'IDR-binding to proteins', 'binding_nucleic': 'IDR-binding to nucleic',
                   'binding_ion': 'IDR-binding to ions', 'binding_lipid': 'IDR-binding to lipids',
                   'binding_SM': 'IDR binding to small molecule'}


def _load_tags_sets_dict(dp_path, tags_list):
    tags_dict = {}
    _p = f'{dp_path}GO/Ontologies_Sets/'
    for tag in tags_list:
        with open(f"{_p}{tag}.set", 'r') as fin:
            for line in fin:
                line = line.strip()
                if line not in tags_dict:
                    tags_dict[line] = []
                tags_dict[line].append(tag)
    return tags_dict


if __name__ == '__main__':
    tags_list = list(tags_names_dict.keys())

    # example of tags_sets_dict: GO:1905576 ['binding', 'binding_ion', 'binding_lipid']
    tags_sets_dict = _load_tags_sets_dict(dp_path=dp_path, tags_list=tags_list)
    json_file = f"{_in_path}DisProt release_{release}.json"
    with open(json_file, 'r') as fin:
        data = json.load(fin)
    print(f"json clean size:\t{len(data['data']):,}")

    names_list = ['DisProt', 'UniProt', 'OX']
    disprot = annotated_fasta(data_name=f'DisProt release_{release}, annotations include descendants',
                              names_list=names_list, accession='DisProt')
    disprot['metadata']['tags_dict'] = tags_names_dict

    for dd in data['data']:
        ac = dd['disprot_id']
        disprot['data'][ac] = {'DisProt': ac, 'UniProt': dd['acc'], 'OX': f"{dd['ncbi_taxon_id']}",
                               'seq': dd['sequence']}

        for tid in tags_list:
            disprot['data'][ac][f"lst_{tid}"] = ['0'] * len(disprot['data'][ac]['seq'])
        for rg in dd['regions']:
            tid = rg['term_id']
            if tid in tags_sets_dict:
                st = int(rg['start']) - 1
                ed = int(rg['end'])
                for tag in tags_sets_dict[tid]:
                    for i in range(st, ed):
                        disprot['data'][ac][f"lst_{tag}"][i] = '1'

    json_file = f"{_in_path}DisProt release_{release} with_ambiguous_evidences.json"
    with open(json_file, 'r') as fin:
        data = json.load(fin)
    print(f"json ambiguous size:\t{len(data['data']):,}")
    for dd in data['data']:
        ac = dd['disprot_id']
        if ac not in disprot['data']:
            continue
        for rg in dd['regions']:
            tid = rg['term_id']
            if tid in tags_sets_dict:
                st = int(rg['start']) - 1
                ed = int(rg['end'])
                for tag in tags_sets_dict[tid]:
                    for i in range(st, ed):
                        if disprot['data'][ac][f"lst_{tag}"][i] == '0':
                            disprot['data'][ac][f"lst_{tag}"][i] = '-'

        for tid in tags_list:
            disprot['data'][ac][tid] = ''.join(disprot['data'][ac][f"lst_{tid}"])

    aff_remove_no_class_tag(af=disprot, tag='IDR', cl='1')

    out_file = f"DisProt_{release}.af"
    aff_save2(f_name=out_file, af=disprot)
