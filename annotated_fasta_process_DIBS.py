import xml.etree.ElementTree as et
from annotated_fasta import *
from miscellaneous import is_float, get_go_term_lineage, get_uniprot_seq
# from miscellaneous import get_url_response

go_obsolete_dict = {'GO:0097159': '', 'GO:0005070': '', 'GO:1990521': '', 'GO:0005057': '', 'GO:0004871': '',
                    'GO:0008022': 'GO:0005515', 'GO:0018024': 'GO:0140938', 'GO:0030374': 'GO:0003713',
                    'GO:0047485': 'GO:0005515', 'GO:0035064': 'GO:0140566', 'GO:0001190': 'GO:0003676',
                    'GO:0001083': 'GO:0003676', 'GO:0001135': '	GO:0003712', 'GO:0001076': 'GO:0003676',
                    'GO:0000989': 'GO:0008134', 'GO:0043621': 'GO:0005515'}

def _get_tags(sz: int=0):
    return {'IDR': '', 'IDR_Partner': '', 'DtoO': '', 'binding': '', 'binding_protein': '', 'binding_nucleic': '',
            'binding_lipid': '', 'binding_SM': '',
            'list': {'IDR': ['0'] * sz,
                     'IDR_Partner': ['-'] * sz,
                     'DtoO': ['0'] * sz,
                     'binding': ['0'] * sz,
                     'binding_protein': ['0'] * sz,
                     'binding_nucleic': ['0'] * sz,
                     'binding_lipid': ['0'] * sz,
                     'binding_SM': ['0'] * sz},
            }


def get_general(en):
    kd = None
    general = en.find('general')
    kd_o = general.find('kd')
    if kd_o is not None:
        kd = kd_o.find('value').text.strip()
    pdb = general.find('pdb_id').text.strip()
    method = general.find('exp_method').text.strip()
    d_statue = general.find('disorder_status').text.strip()
    return pdb, method, d_statue, kd


def get_functions(en):
    tag_set = {'binding'}
    func = en.find('function')
    mf = func.find('molecular_function')
    if mf is None:
        return tag_set
    for go in mf.findall('go'):
        go_acc_o = go.find('accession')
        if go_acc_o is None:
            continue
        go_acc = go_acc_o.text.strip()
        if go_acc in go_obsolete_dict:
            go_acc = go_obsolete_dict[go_acc]
            if len(go_acc) < 3:
                continue
        ancestors = get_go_term_lineage(go_acc)
        if ancestors is None:
            print(f"==========\tBAD: {go_acc}")
            continue
        if 'GO:0003674' not in ancestors:
            continue
        # if 'GO:0005488' in ancestors:
        # tag_set.add('binding')
        if 'GO:0003676' in ancestors:
            tag_set.add('binding_nucleic')
        if 'GO:0005515' in ancestors:
            tag_set.add('binding_protein')
        if 'GO:0008289' in ancestors:
            tag_set.add('binding_lipid')
        if 'GO:0036094' in ancestors:
            tag_set.add('binding_SM')
    return tag_set


def aff_dibs_to_af(in_file: str=None, q_use_list: list=None, partners_file: str=None):
    fout = None
    if partners_file:
        fout = open(partners_file, 'w')
    if q_use_list is None:
        q_use_list = ['Confirmed', 'Inferred from motif', 'Inferred from homology']
    database_list = ['srcUniProt', 'OX', 'UniProt']
    tags_dict = {'IDR': 'disorder', 'IDR_Partner': 'IDR partner', 'DtoO': 'Disorder to Order', 'binding': 'Binding',
                 'binding_protein': 'protein bind', 'binding_nucleic': 'nucleic bind', 'binding_lipid': 'lipid binding',
                 'binding_SM': ''}
    af = annotated_fasta(database_list=database_list, tags_dict=tags_dict)
    root = None
    try:
        tree = et.parse(in_file)
        root = tree.getroot()
    except FileNotFoundError:
        print("File not found.")
        exit(1)
    except et.ParseError:
        print("Invalid XML format.")
        exit(1)

    cnt = 0
    e_set = set()
    for entry in root.findall('entry'):
        # entry contain: accession, general, function, macromolecules, evidence
        cnt += 1
        # ================= accession
        acc = entry.find('accession').text
        if acc not in e_set:
            e_set.add(acc)
        else:
            print(acc)
        # ================= general
        pdb, method, d_statue, kd = get_general(en=entry)
        # print(pdb, method, d_statue, kd, flush=True)
        # ================= function
        # function for tags based on molecular_function
        tag_set = get_functions(en=entry)
        # print(list(tag_set), flush=True)
        # ================= macromolecules for chains
        macromolecules = entry.find('macromolecules')
        chain_dict = {'Disordered': [], 'Ordered': []}
        for chain_o in macromolecules.findall('chain'):
            c_id = chain_o.find('id').text  # key
            c_type = chain_o.find('type').text.strip()  # 1
            if 'Ordered' in c_type:
                c_type = 'Ordered'
            if c_type not in ['Ordered', 'Disordered']:
                print('BAD c_type', acc, c_type)
            # c_seq = chain_o.find('sequence').text.strip()  # 2
            c_up_o = chain_o.find('uniprot')
            c_up_ac = c_up_o.find('id').text.strip()   # 3
            c_up_st = int(c_up_o.find('start').text.strip())  # 4
            c_up_ed = int(c_up_o.find('end').text.strip())  # 5
            if c_up_ac not in af['data']:
                seq, ox = get_uniprot_seq(c_up_ac)
                if seq is not None:
                    ox_lst = []
                    if ox is not None:
                        ox_lst = [ox]
                    af['data'][c_up_ac] = {'seq': seq, 'tags': _get_tags(sz=len(seq)),
                                           'databases': {'srcUniProt': [c_up_ac], 'OX': ox_lst, 'UniProt': [c_up_ac]},
                                           'scores': {}}
                else:
                    print(acc, pdb, c_id, c_up_ac)
                    continue  # next chain
            if c_up_ed > len(af['data'][c_up_ac]['seq']):
                print(f"ERROR: Bad region in protein {c_up_ed}, PDB {pdb} chain {c_id}", flush=True)
                continue  # next chain
            chain_dict[c_type].append({'UP': c_up_ac, 'start': c_up_st, 'end': c_up_ed})
            if c_type == 'Disordered':
                tag_lst = ['IDR'] + list(tag_set)
            else:
                tag_lst = ['IDR_Partner']

            for tag in tag_lst:
                # print(c_up_ac, tag, c_up_st-1, c_up_ed, flush=True)
                for ii in range(c_up_st-1, c_up_ed):
                    af['data'][c_up_ac]['tags']['list'][tag][ii] = '1'

            regions_o = chain_o.find('regions')
            for rg in regions_o:
                r_type = rg.find('region_type').text.strip()
                r_st = int(rg.find('region_start').text.strip())
                r_ed = int(rg.find('region_end').text.strip())
                if r_type != 'secondary structure':
                    continue
                if c_type == 'Ordered':
                    for ii in range(r_st-1, r_ed):
                        af['data'][c_up_ac]['tags']['list']['IDR_Partner'][ii] = '0'
                elif c_type == 'Disordered':
                    for ii in range(r_st-1, r_ed):
                        if d_statue in q_use_list:
                            af['data'][c_up_ac]['tags']['list']['DtoO'][ii] = '1'
                        else:
                            if af['data'][c_up_ac]['tags']['list']['DtoO'][ii] != '1':
                                af['data'][c_up_ac]['tags']['list']['DtoO'][ii] = '-'
        if fout:
            if len(chain_dict['Disordered']) == 1:
                dis_up = chain_dict['Disordered'][0]['UP']
                dis_st = chain_dict['Disordered'][0]['start']
                dis_ed = chain_dict['Disordered'][0]['end']
                for ch_or in chain_dict['Ordered']:
                    print(f"{acc}\t{dis_up}\t{dis_st}\t{dis_ed}", end='\t', file=fout)
                    print(f"{ch_or['UP']}\t{ch_or['start']}\t{ch_or['end']}\t{kd}", file=fout, flush=True)
            else:
                print(f"{acc}\tlen(chain_dict['Disordered'])\t{len(chain_dict['Disordered'])}", flush=True)

    if fout:
        fout.close()
    for ac in af['data']:
        for tg in af['metadata']['tags_list']:
            af['data'][ac]['tags'][tg] = ''.join(af['data'][ac]['tags']['list'][tg])
        del af['data'][ac]['tags']['list']
    return af

#
