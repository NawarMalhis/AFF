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
    return {'IDR': '', 'DtoO': '', 'binding': '', 'binding_protein': '', 'binding_nucleic': '',
            'binding_lipid': '', 'binding_SM': '',
            'list': {'IDR': ['-'] * sz,
                     'DtoO': ['-'] * sz,
                     'binding': ['-'] * sz,
                     'binding_protein': ['-'] * sz,
                     'binding_nucleic': ['-'] * sz,
                     'binding_lipid': ['-'] * sz,
                     'binding_SM': ['-'] * sz},
            }


def aff_dibs_to_af(in_file: str=None):
    database_list = ['srcUniProt', 'OX', 'UniProt']
    af = annotated_fasta(database_list=database_list)
    af['metadata']['tags_dict'] = {'IDR': 'DIBS disorder', 'P_Bind': 'DIBS protein bind'}
    af['metadata']['tags_list'] = ['IDR', 'P_Bind']
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
    # m_func_set = set()

    for entry in root.findall('entry'):
        # entry contain: accession, general, function, macromolecules, evidence
        tag_set = set()
        cnt += 1
        # ================= accession
        acc = entry.find('accession').text

        # ================= general
        general = entry.find('general')
        pdb = general.find('pdb_id').text
        method = general.find('exp_method').text
        d_statue = general.find('disorder_status').text
        # ================= function
        # function for tags based on molecular_function
        func = entry.find('function')
        mf = func.find('molecular_function')
        if mf is None:
            continue  # next entry
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
            if 'GO:0005488' in ancestors:
                tag_set.add('binding')
            if 'GO:0003676' in ancestors:
                tag_set.add('binding_nucleic')
            if 'GO:0005515' in ancestors:
                tag_set.add('binding_protein')
            if 'GO:0008289' in ancestors:
                tag_set.add('binding_lipid')
            if 'GO:0036094' in ancestors:
                tag_set.add('binding_SM')
        # ================= macromolecules
        macromolecules = entry.find('macromolecules')
        mm_general = macromolecules.find('general')
        nr_of_chains_o = mm_general.find('nr_of_chains')
        if nr_of_chains_o is None:
            continue  # next entry
        nr_of_chains = int(nr_of_chains_o.text)
        chain_dict = {}
        for chain_o in macromolecules.findall('chain'):
            c_id = chain_o.find('id').text
            c_type = chain_o.find('type').text.strip()
            c_seq = chain_o.find('sequence').text.strip()
            c_up_o = chain_o.find('uniprot')
            c_up_ac = c_up_o.find('id')
            c_up_st = int(c_up_o.find('start').text.strip())
            c_up_ed = int(c_up_o.find('end').text.strip())
            chain_dict[c_id] = {'up_ac': c_up_ac, 'type': c_type, 'seq': c_seq, 'start': c_up_st, 'end': c_up_ed}
            if c_up_ac not in af['data']:
                seq, ox = get_uniprot_seq(c_up_ac)
                if seq is not None:
                    ox_lst = []
                    if ox is not None:
                        ox_lst = [ox]
                    af['data'][c_up_ac] = {'seq': seq, 'tags': _get_tags(sz=len(seq)),
                                           'databases': {'srcUniProt': [c_up_ac], 'OX': ox_lst, 'UniProt': []},
                                           'scores': {}}

            regions_o = chain_o.find('regions')
            for rg in regions_o:

                pass



        print(cnt, acc, pdb, method, d_statue, list(tag_set))

    # seq, ox = get_uniprot_seq(ac)

    # for ii, mf in enumerate(m_func_set):
    #     print(ii, mf)
