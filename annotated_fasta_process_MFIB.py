import xml.etree.ElementTree as et
from annotated_fasta import *
from miscellaneous import is_float


def _get_tags(sz: int=0):
    return {'IDR': '', 'binding_protein': '',
            'list': {'IDR': ['-'] * sz,
                     'binding_protein': ['0'] * sz}
            }


def get_general(en):
    general = en.find('general')
    pdb = general.find('pdb_id').text.strip()
    method = general.find('exp_method').text.strip()
    return pdb, method


def aff_mfib_to_af(in_file, af):
    try:
        tree = et.parse(in_file)
        root = tree.getroot()
    except FileNotFoundError:
        print(f"File not found:\t{in_file}")
        exit(1)
    except et.ParseError:
        print(f"Invalid XML format:\t{in_file}")
        exit(1)
    acc = root.find('accession').text
    # pdb, method = get_general(en=entry)
    # ================= macromolecules for chains
    macromolecules = root.find('macromolecules')
    for chain_o in macromolecules.findall('chain'):
        c_id = chain_o.find('id').text  # key
        c_up_o = chain_o.find('uniprot')
        c_up_ac = c_up_o.find('id').text.strip()  # 3
        c_up_st = int(c_up_o.find('start').text.strip())  # 4
        c_up_ed = int(c_up_o.find('end').text.strip())  # 5
        c_up_seq = c_up_o.find('sequence').text.strip()  # 2
        sz = len(c_up_seq)
        if c_up_ac not in af['data']:
            af['data'][c_up_ac] = {'seq': c_up_seq, 'tags': _get_tags(sz=sz), 'scores': {},
                                   'databases': {'srcUniProt': [c_up_ac]}}
        for ii in range(c_up_st-1, c_up_ed):
            af['data'][c_up_ac]['tags']['list']['IDR'][ii] = '1'
            af['data'][c_up_ac]['tags']['list']['binding_protein'][ii] = '1'
        # regions_o = chain_o.find('regions')
        # for rg in regions_o:
        #     r_type = rg.find('region_type').text.strip()
        #     r_st = int(rg.find('region_start').text.strip())
        #     r_ed = int(rg.find('region_end').text.strip())
        #     if r_type != 'secondary structure':
        #         continue
        #     for ii in range(r_st - 1, r_ed):
        #         print(acc, c_up_ac, sz, ii, flush=True)
        #         af['data'][c_up_ac]['tags']['list']['DtoO'][ii] = '1'
        for tg in af['metadata']['tags_list']:
            af['data'][c_up_ac]['tags'][tg] = ''.join(af['data'][c_up_ac]['tags']['list'][tg])
    return
