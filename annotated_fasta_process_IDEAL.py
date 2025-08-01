import xml.etree.ElementTree as et
from annotated_fasta import *
from miscellaneous import is_float


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
                    # print(f"pH:\t{idp_id}\t{ph_o.text}", flush=True)
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


def get_ideal_pros(idp_id, seq, en, typ_dict, pros_set):
    d_to_o = ['-'] * len(seq)
    for fp in en.findall('Function_pros'):
        pt_o = fp.find('pros_type')
        pros_set.add(idp_id)
        if pt_o is None:
            continue
        pt = pt_o.text.strip()
        if pt not in typ_dict:
            typ_dict[pt] = {idp_id}
        else:
            typ_dict[pt].add(idp_id)
        if pt != 'verified':
            continue
        or_o = fp.find('order_location')
        if or_o is None:
            continue

        st_o = or_o.find('order_region_start')
        ed_o = or_o.find('order_region_end')
        com_o = fp.find('comment')
        if st_o is None or ed_o is None:
            continue
        if com_o is not None:
            com_t = com_o.text.strip()
            if com_t in com_set:
                continue

        st = int(st_o.text) - 1
        ed = int(ed_o.text) - 1
        for ii in range(st, ed):
            d_to_o[ii] = '1'
    return ''.join(d_to_o)


def get_ideal_idr(idp_id, seq, en, temp_range, ph_range):
    condition_id_set = get_ideal_condition_set(en=en, temp_range=temp_range, ph_range=ph_range, idp_id=idp_id)
    idr_list = ['-'] * len(seq)
    for rg in en.findall('Region'):
        st = int(rg.find('region_start').text) - 1
        ed = int(rg.find('region_end').text) - 1
        od = rg.find('order_disorder').text
        cid = rg.find('condition_id').text
        if st < 0 or st >= len(seq) or ed < 0 or ed >= len(seq):
            print(f"{idp_id}\tSize {len(seq)}\t{st}\t{ed}", flush=True)
            exit(0)
        if cid not in condition_id_set:
            continue
        if od == 'disorder':
            for ii in range(st, ed):
                idr_list[ii] = '1'
        elif od == 'order':
            for ii in range(st, ed):
                if idr_list[ii] != '1':
                    idr_list[ii] = '0'
    return ''.join(idr_list)


def aff_ideal_to_af(in_file: str=None, temp_range: list=None, ph_range: list=None):
    pros_set = set()
    typ_dict = {}
    database_list = ['srcUniProt']
    af = annotated_fasta(database_list=database_list)
    af['metadata']['tags_dict'] = {'IDR': 'IDEAL disorder', 'P_Bind': 'IDEAL protein bind'}
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

    # {'version': 1, 'IDEAL_entry': 1110, 'IDEAL_BiologicalProcess': 177, 'IDEAL_interaction': 1089}
    for entry in root.findall('IDEAL_entry'):
        idp_id = entry.find('idp_id').text
        general = entry.find('General')
        uniprot_list = []
        for up in general.findall('uniprot'):
            uniprot_list.append(up.text)
        seq = general.find('sequence').text
        af['data'][idp_id] = {'seq': seq, 'tags': {}, 'databases': {'srcUniProt': uniprot_list}, 'scores': {}}
        af['data'][idp_id]['tags']['IDR'] = get_ideal_idr(idp_id=idp_id, seq=seq, en=entry, temp_range=temp_range,
                                                          ph_range=ph_range)
        af['data'][idp_id]['tags']['P_Bind'] = get_ideal_pros(idp_id=idp_id, seq=seq, en=entry, pros_set=pros_set,
                                                              typ_dict=typ_dict)
    print(f"len(pros_set):\t{len(pros_set)}\t{len(typ_dict['verified'].union(typ_dict['possible']))}")
    for typ in typ_dict:
        print(f"{typ}:\t{len(typ_dict[typ])}")
    return af

