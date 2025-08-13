import xml.etree.ElementTree as et
from annotated_fasta import *
from miscellaneous import *  # is_float, get_xml_root, dbs_database_list, dbs_tags_dict

com_set = {'This region is ordered when it binds upon DNA.',
           'This region is ordered in the complex with DNA.',
           'This region is ordered in the DNA complex. The disorder-to-order transition locks pol k around the DNA.',
           'This region interacts with DNA',
           'This region is ordered when it binds upon natural IR-1 element (DNA).',
           'This region binds to RNA.',
           'This region binds upon DNA.',
           'This region is ordered in the complex with ssDNA.',
           'This region is ordered in the complex with RNA.',
           'This region is ordered in the DNA complex.'
           }


def get_ideal_condition_set(en, temp_range, ph_range, idp_id):
    condition_id_set = set()
    for cnd in en.findall('Condition'):
        cid = cnd.find('condition_id').text
        cmd_o = cnd.find('method')
        if cmd_o is None:
            continue
        cmd = cmd_o.text
        _ok = False
        if cmd == 'X-RAY':
            _ok = True
            cc = cnd.find('crystallization_condition')
            if cc is None:
                _ok = False
                continue
            tmp_o = cc.find('temperature')
            ph_o = cc.find('pH')
            if tmp_o is not None:
                tmp_t = tmp_o.text.split('K')[0]
                if not is_float(tmp_t):
                    _ok = False
                    # print(f"temperature:\t{idp_id}\t{tmp_o.text}", flush=True)
                else:
                    tmp = float(tmp_t)
                    if tmp < temp_range[0] or tmp > temp_range[1]:
                        _ok = False

            if ph_o is not None:
                ph_t = ph_o.text.strip().split('-')[0]
                ph_t = ph_t.split(')')[0].split('.')[0]
                if not is_float(ph_t):
                    _ok = False
                else:
                    ph = float(ph_t)
                    if ph < ph_range[0] or ph > ph_range[1]:
                        _ok = False
        elif cmd == 'NMR':
            _ok = True
        if _ok:
            condition_id_set.add(cid)
    return condition_id_set


def get_ideal_pros(af, idp_id, en, typ_dict, lq_list):
    for fp in en.findall('Function_pros'):
        phq = True
        pt_o = fp.find('pros_type')
        if pt_o is None:
            print(f"{idp_id} pros_typ is None", flush=True)
            continue
        pt = pt_o.text.strip()
        if pt not in typ_dict:
            typ_dict[pt] = {idp_id}
        else:
            typ_dict[pt].add(idp_id)
        if pt in lq_list:
            phq = False
        or_o = fp.find('order_location')
        if or_o is None:
            print(f"{idp_id}\tNone order_location", flush=True)
            continue

        st_o = or_o.find('order_region_start')
        ed_o = or_o.find('order_region_end')
        com_o = fp.find('comment')
        if st_o is None or ed_o is None:
            print(f"{idp_id}\tNone order_region_start/end", flush=True)
            continue

        if com_o is not None:
            com_t = com_o.text.strip()
            if com_t in com_set:
                print(f"{idp_id}\tComment:\t{com_t}", flush=True)
                continue
        rid_list = or_o.find('order_region_id').text.split(',')
        # print(idp_id, rid_list, flush=True)
        # for rid in rid_list:
        #     if rid in regions_dict[idp_id]:
        #         regions_dict[idp_id][rid]['pros'] = pt
        st = int(st_o.text) - 1
        ed = int(ed_o.text)
        pros = {'start': st, 'end': ed}
        # if idp_id not in regions_dict:
        #     regions_dict[idp_id] = [pros]
        # else:
        #     regions_dict[idp_id].append(pros)
        for ii in range(st, ed):
            if phq:
                af['data'][idp_id]['tags']['list']['binding_protein'][ii] = '1'
                af['data'][idp_id]['tags']['list']['DtoO'][ii] = '1'
                # af['data'][idp_id]['tags']['list']['IDR'][ii] = '1'
            elif af['data'][idp_id]['tags']['list']['binding_protein'][ii] != '1':
                af['data'][idp_id]['tags']['list']['binding_protein'][ii] = '-'
                af['data'][idp_id]['tags']['list']['DtoO'][ii] = '-'
                # if af['data'][idp_id]['tags']['list']['IDR'][ii] not in '1':
                #     af['data'][idp_id]['tags']['list']['IDR'][ii] = '-'

    return


# e_set = {'IID00004', 'IID00010', 'IID00028', 'IID00067', 'IID00077', 'IID00099', 'IID00107', 'IID00116', 'IID00121'}
def process_ideal_idr_all(af, root, temp_range, ph_range):
    # regions_dict = {}
    # bad_pdb_set = set()
    for entry in root.findall('IDEAL_entry'):
        idp_id = get_general(af=af, en=entry)

        # to fill IDR (originally with '-'s) with '0's and '1's
        get_ideal_idr(af=af, idp_id=idp_id, en=entry, temp_range=temp_range, ph_range=ph_range)

    return


def process_pros_all(af, root, lq_list, verbose=False):
    typ_dict = {}
    for entry in root.findall('IDEAL_entry'):
        # to fill DtoO and binding_protein (originally with '0's) with '1' and '-' and also to populate pros_dict
        idp_id = entry.find('idp_id').text
        get_ideal_pros(af=af, idp_id=idp_id, en=entry, typ_dict=typ_dict, lq_list=lq_list)
    if verbose:
        for typ in typ_dict:
            print(f"{typ}:\t{len(typ_dict[typ])}")



def get_ideal_idr(af, idp_id, en, temp_range, ph_range):
    condition_id_set = get_ideal_condition_set(en=en, temp_range=temp_range, ph_range=ph_range, idp_id=idp_id)

    seq = af['data'][idp_id]['seq']
    for rg in en.findall('Region'):
        od = rg.find('order_disorder').text
        if od not in ['order', 'disorder']:
            continue
        cond_id = rg.find('condition_id').text
        if cond_id not in condition_id_set:
            continue
        st = int(rg.find('region_start').text)
        ed = int(rg.find('region_end').text)

        if st < 1 or st > len(seq) or ed < 1 or ed > len(seq):
            print(f"{idp_id}\tSize {len(seq)}\t{st}\t{ed}", flush=True)
            exit(1)

        # -----------------------------------------------------------------------------------------------
        # key is region_id
        chain_id_o = rg.find('chain_id')
        if chain_id_o is None:
            chain_id = None
        else:
            chain_id = chain_id_o.text.strip()
        rri = rg.find('region_id').text.strip()
        # regions_dict[idp_id][rri] = {'chain_id': chain_id, 'start': st, 'end': ed, 'type': od, 'pros': '-------'}
        # ----------------------------------------------------------------------------------------------
        if od == 'disorder':
            for ii in range(st - 1, ed):
                af['data'][idp_id]['tags']['list']['IDR'][ii] = '1'
        elif od == 'order':
            for ii in range(st - 1, ed):
                if af['data'][idp_id]['tags']['list']['IDR'][ii] != '1':
                    af['data'][idp_id]['tags']['list']['IDR'][ii] = '0'
    return


def get_general(af, en):
    idp_id = en.find('idp_id').text
    general = en.find('General')
    uniprot_list = []
    for up in general.findall('uniprot'):
        uniprot_list.append(up.text)
    seq = general.find('sequence').text
    sz = len(seq)
    af['data'][idp_id] = {'seq': seq, 'tags': get_dbs_ac_tags(sz), 'databases': {'srcUniProt': uniprot_list},
                          'scores': {}}
    return idp_id


def process_ne_proc_all(af, root):
    for entry in root.findall('IDEAL_entry'):
        ok = True
        # to fill DtoO and binding_protein (originally with '0's) with '1' and '-' and also to populate pros_dict
        idp_id = entry.find('idp_id').text
        if idp_id not in af['data']:
            continue
        lst = ['-'] * len(af['data'][idp_id]['seq'])
        np = entry.find('NeProc')
        order_o = np.find('order')
        if order_o is not None:
            order_pair_list = order_o.text.split(',')
            for ord_pair in order_pair_list:
                ss, ee = ord_pair.split('-')
                st = int(ss)
                ed = int(ee)
                # print(idp_id, 'order', len(lst), st, ed, flush=True)
                if st < 1 or ed > len(lst):
                    ok = False
                    continue
                for ii in range(st-1, ed):
                    lst[ii] = '0'

        disorder_o = np.find('disorder')
        if disorder_o is not None:
            disorder_pair_list = disorder_o.text.split(',')
            for dis_pair in disorder_pair_list:
                ss, ee = dis_pair.split('-')
                st = int(ss)
                ed = int(ee)
                if st < 1 or ed > len(lst):
                    ok = False
                    continue
                for ii in range(st-1, ed):
                    lst[ii] = '1'
        if ok:
            af['data'][idp_id]['tags']['NP'] = ''.join(lst)
    return


def get_interaction_partners(af, root):
    interactions_dict = {}
    for ac in af['data']:
        seq = af['data'][ac]['seq']
        for ii in range(len(seq)):
            if af['data'][ac]['tags']['list']['IDR'][ii] == '1' or \
                    af['data'][ac]['tags']['list']['DtoO'][ii] == '1':
                af['data'][ac]['tags']['list']['IDR_partner'][ii] = '1'
            elif af['data'][ac]['tags']['list']['IDR'][ii] == '0' or \
                    af['data'][ac]['tags']['list']['DtoO'][ii] == '0' :
                af['data'][ac]['tags']['list']['IDR_partner'][ii] = '0'

    ac_list = list(af['data'].keys())
    for ac in ac_list:
        if af['data'][ac]['tags']['list']['IDR_partner'].count('1') == 0:
            del af['data'][ac]
    ac_list = list(af['data'].keys())

    for interaction in root.findall('IDEAL_interaction'):
        int_id = interaction.find('interaction_id').text.strip()
        id1_o = interaction.find('IDEAL_entry_1')
        id1 = id1_o.find('idp_id').text.strip()
        id2_o = interaction.find('IDEAL_entry_2')
        id2 = id2_o.find('idp_id').text.strip()
        if id1 not in ac_list or id2 not in ac_list:
            continue
        rg1_st_list = []
        rg1_ed_list = []
        rg2_st_list = []
        rg2_ed_list = []
        for rg_o in id1_o.findall('Region'):
            st = rg_o.find('region_start').text
            ed = rg_o.find('region_end').text
            rg1_st_list.append(int(st))
            rg1_ed_list.append(int(ed))

        for rg_o in id2_o.findall('Region'):
            st = rg_o.find('region_start').text
            ed = rg_o.find('region_end').text
            rg2_st_list.append(int(st))
            rg2_ed_list.append(int(ed))
        lst1 = ['-'] * len(af['data'][id1]['seq'])
        lst2 = ['-'] * len(af['data'][id2]['seq'])

        for iii in range(len(rg1_st_list)):
            st = rg1_st_list[iii]
            ed = rg1_ed_list[iii]
            for ii in range(st-1, ed):
                lst1[ii] = '1'
        st1 = ''.join(lst1)

        for iii in range(len(rg2_st_list)):
            st = rg2_st_list[iii]
            ed = rg2_ed_list[iii]
            for ii in range(st-1, ed):
                lst2[ii] = '1'
        st2 = ''.join(lst2)
        interactions_dict[int_id] = {id1: st1, id2: st2}
    return interactions_dict

def aff_ideal_to_af(in_file: str = None, temp_range: list = None, ph_range: list = None, lq_list: list = None):
    if lq_list is None:
        lq_list = ['predicted']
    af = annotated_fasta(database_list=dbs_database_list, tags_dict=get_dbs_tags_dict(source='IDEAL'))
    root = get_xml_root(xml_file=in_file)
    # {'version': 1, 'IDEAL_entry': 1110, 'IDEAL_BiologicalProcess': 177, 'IDEAL_interaction': 1089}

    process_ideal_idr_all(af, root, temp_range=temp_range, ph_range=ph_range)
    process_pros_all(af, root, lq_list=lq_list, verbose=True)
    process_ne_proc_all(af, root)  # IDR prediction NeProc
    interactions_dict = get_interaction_partners(af, root)

    for tg in af['metadata']['tags_list']:
        for ac in af['data']:
            af['data'][ac]['tags'][tg] = ''.join(af['data'][ac]['tags']['list'][tg])

    for ac in af['data']:
        if 'NP' in af['data'][ac]['tags']:
            if af['data'][ac]['tags']['NP'][-1] == '-':
                af['data'][ac]['tags']['NP'] = '-' * len(af['data'][ac]['seq'])
        else:
            af['data'][ac]['tags']['NP'] = '-' * len(af['data'][ac]['seq'])

    return af, interactions_dict
