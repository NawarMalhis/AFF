from annotated_fasta import *
import json

# The .obo file format.
# obo is short for (Open Biomedical Ontologies)
# go-basic.obo source: https://current.geneontology.org/ontology/go-basic.obo
# IDPO_v0.3.0.obo source is DisProt
def aff_process_go(dp_path):
    # go_path = '/home/nmalhis/data/DisProt/GO/'
    mf_base_dict = {'GO:0005488': 'binding', 'GO:0005515': 'binding_protein', 'GO:0003676': 'binding_nucleic',
                    'GO:0003723': 'binding_RNA', 'GO:0003677': 'binding_DNA', 'GO:0008289': 'binding_lipid',
                    'GO:0043167': 'binding_ion', ' GO:0036094': 'binding_SM'}  # SM is for 'small molecule'
    n_space_dict = {'molecular_function': {'in_obo': 'go-basic.obo', 'base_dict': mf_base_dict},
                    'structural_state': {'in_obo': 'IDPO_v0.3.0.obo', 'base_dict': {'IDPO:00076': 'IDR'}},
                    'disorder_function': {'in_obo': 'IDPO_v0.3.0.obo', 'base_dict': {'IDPO:00502': 'Linker'}},
                    'structural_transition': {'in_obo': 'IDPO_v0.3.0.obo', 'base_dict': {'IDPO:00050': 'DtO'}}}

    for n_space in n_space_dict:
        base_dict = n_space_dict[n_space]['base_dict']
        in_obo = n_space_dict[n_space]['in_obo']
        ancestors_dict = {}
        with open(f"{dp_path}GO/OBO/{in_obo}", 'r') as fin:
            in_term = False
            namespace = ''
            go_id = ''
            for line in fin:
                line = line.strip()
                if '[Term]' in line:
                    in_term = True
                if not in_term:
                    continue
                if line[:3] == 'id:':
                    go_id = line.split()[1]
                if 'namespace:' == line[:10]:
                    namespace = line.split()[1]
                    if namespace == n_space:
                        ancestors_dict[go_id] = []
                if namespace != n_space:
                    continue
                if 'is_a:' == line[:5]:
                    ancestors_dict[go_id].append(line.split()[1])

        # remove non molecular_function ancestors
        id_list = list(ancestors_dict.keys())
        for go_id in id_list:
            is_a = []
            for aa in ancestors_dict[go_id]:
                if aa in id_list:
                    is_a.append(aa)
            ancestors_dict[go_id] = is_a

        # base_dicts are the targets
        for go_base in base_dict:
            descendants_set = set()
            descendants_set.add(go_base)
            sz = len(descendants_set)
            while True:
                for go_id in ancestors_dict:
                    for isa in ancestors_dict[go_id]:
                        if isa in descendants_set:
                            descendants_set.add(go_id)

                if sz == len(descendants_set):
                    print(f'{base_dict[go_base]}\t{go_base} set size:', len(descendants_set))
                    break
                sz = len(descendants_set)

            with open(f"{dp_path}GO/{base_dict[go_base]}.set", 'w') as fout:
                for go_id in descendants_set:
                    print(go_id, file=fout)


def aff_process_disprot(release, dp_path):
    tags_names_dict = {'IDR': 'Protein disordered region', 'DtO': 'Disordered to ordered transition',
                       'Linker': 'Linker regions', 'binding': 'IDR binding in general',
                       'binding_protein': 'IDR-binding to proteins', 'binding_nucleic': 'IDR-binding to nucleic',
                       'binding_ion': 'IDR-binding to ions', 'binding_lipid': 'IDR-binding to lipids',
                       'binding_SM': 'IDR binding to small molecule'}

    tags_list = list(tags_names_dict.keys())

    # example of tags_sets_dict: GO:1905576 ['binding', 'binding_ion', 'binding_lipid']
    tags_sets_dict = _load_tags_sets_dict(dp_path=dp_path, tags_list=tags_list)
    json_file = f"{dp_path}JSON/DisProt release_{release}.json"
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

    json_file = f"{dp_path}JSON/DisProt release_{release} with_ambiguous_evidences.json"
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

    out_file = f"{dp_path}DisProt_{release}.af"
    aff_save2(f_name=out_file, af=disprot)


def _load_tags_sets_dict(dp_path, tags_list):
    tags_dict = {}
    _p = f'{dp_path}GO/'
    for tag in tags_list:
        with open(f"{_p}{tag}.set", 'r') as fin:
            for line in fin:
                line = line.strip()
                if line not in tags_dict:
                    tags_dict[line] = []
                tags_dict[line].append(tag)
    return tags_dict

