import xml.etree.ElementTree as et
from annotated_fasta import *
from miscellaneous import is_float, get_xml_root


def _get_tags(sz: int=0):
    return {'IDR': '', 'binding_protein': '',
            'list': {'IDR': ['-'] * sz,
                     'binding_protein': ['0'] * sz}
            }


def aff_fuzdb_to_af(in_file):
    database_list = ['srcUniProt']
    tags_dict = {'IDR': 'disorder', 'binding_protein': 'protein bind'}
    af = annotated_fasta(database_list=database_list, tags_dict=tags_dict)
    root = get_xml_root(xml_file=in_file)

    cnt = 0
    for fdb in root.findall('fuzdb'):
        cnt += 1
        fz_id = fdb.find('entry_id').text.strip()
        ac = fdb.find('uniprot_acc').text.strip()
        seq = fdb.find('sequence').text.strip()
        sz = len(seq)
        if ac not in af['data']:
            af['data'][ac] = {'seq': seq, 'tags': _get_tags(sz=sz), 'databases': {'srcUniProt': [ac]}, 'scores': {}}
        for frg_o in fdb.findall('fuzzy_region'):
            st_text = frg_o.find('start').text
            ed_text = frg_o.find('end').text
            if st_text == 'null' or ed_text == 'null':
                continue
            print(fz_id, ac, st_text, ed_text, flush=True)
            fr_st = int(st_text)
            fr_ed = int(ed_text)
            for ii in range(fr_st - 1, fr_ed):
                af['data'][ac]['tags']['list']['IDR'][ii] = '1'
                af['data'][ac]['tags']['list']['binding_protein'][ii] = '1'
    for ac in af['data']:
        for tg in af['data'][ac]['tags']['list']:
            af['data'][ac]['tags'][tg] = ''.join(af['data'][ac]['tags']['list'][tg])
    return af