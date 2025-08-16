import requests
import xml.etree.ElementTree as et
import json
import os


dbs_database_list = ['srcUniProt']

def get_dbs_tags_dict(source):
    return {'IDR': f'{source} disorder', 'DtoO': f'{source} disorder to order transition',
            'binding_protein': f'{source} protein bind', 'binding_nucleic': f'{source} nucleic bind',
            'IDR_partner': f'{source} partner disorder', 'binding_partner': f'{source} partner protein bind'
            }

def get_dbs_ac_tags(sz):
    return {'IDR': '-' * sz, 'DtoO': '-' * sz, 'binding_protein': '-' * sz, 'binding_nucleic': '-' * sz,
            'IDR_partner': '-' * sz, 'binding_partner': '-' * sz,
            'list': {'IDR': ['-'] * sz, 'DtoO': ['0'] * sz, 'binding_protein': ['0'] * sz,
                     'binding_nucleic': ['0'] * sz, 'IDR_partner': ['-'] * sz, 'binding_partner': ['0'] * sz
                     }
            }


def get_url_response(url, **kwargs):
    response = None
    for attempt in range(1, 11):
        try:
            response = requests.get(url, **kwargs)
            if response is not None:
                break
            else: print(f"Attempt {attempt} failed, response is None:\t{url}")
        except Exception as e:
            print(f"Attempt {attempt} failed: {url}")

    if not response.ok:
        print(response.text)
        return None
    return response


def is_float(text):
    try:
        float(text)
        return True
    except ValueError:
        return False


def get_go_term_lineage(go_id, verbose=False):
    url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_id}/ancestors?relations=is_a,part_of"
    response = get_url_response(url)  # requests.get(url, headers={"Accept": "application/json"})
    if response is None:
        return None
    # get_url_response
    if response.status_code == 200:
        data = response.json()
        if 'ancestors' not in data['results'][0]:
            print("KeyError: 'ancestors'")
            return None
        ancestors = data['results'][0]['ancestors']
        if verbose:
            print(f"GO Term: {go_id}")
            print("Ancestors:", ancestors)
        return ancestors
    else:
        print(f"Failed to retrieve data for {go_id}")
    return None


def get_xml_root(xml_file):
    root = None
    try:
        tree = et.parse(xml_file)
        root = tree.getroot()
    except FileNotFoundError:
        print(f"{xml_file}\tFile not found.")
        exit(1)
    except et.ParseError:
        print(f"{xml_file}\tInvalid XML format.")
        exit(1)
    return root

def load_json(file_path):
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
        return data
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON: {e}")
        return None
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None


def get_uniprot_seq(ac):
    seq = None
    ox = None
    url = f'https://rest.uniprot.org/uniprotkb/{ac}.json'
    response = get_url_response(url)
    if response is None:
        print("_process_uniprot_list, response is None")
        return seq, ox
    data = response.json()

    ss = data.get('sequence')
    if not ss:
        return ss, ox
    seq = ss['value']

    _taxa = data.get('organism')
    if not _taxa:
        return ss, ox
    if 'taxonId' not in _taxa:
        return ss, ox

    ox = str(_taxa['taxonId'])

    return seq, ox


def run_iup(fasta_path, iup_path, out_path):
    f_list = os.listdir(fasta_path)
    for ff in f_list:
        cmd = f"python3 {iup_path}iupred3.py -a {fasta_path}{ff} short > {out_path}{ff.split('.')[0]}.iup3"
        os.system(cmd)


def split_fasta(in_file, out_path):
    from annotated_fasta import aff_load_fasta
    af = aff_load_fasta(in_file)
    for ac in af['data']:
        o_file = f"{out_path}{ac}.fasta"
        with open(o_file, 'w') as fout:
            print(f">{ac}\n{af['data'][ac]['seq']}", file=fout)

def split_morf_chibi(in_file, out_path):
    prd_list = ['MCW', 'MCL', 'MC']
    # fout = {'MCW': None, 'MCL': None, 'MC': None}  # , 'IDP': None}
    ac = ''
    fout = {}
    for prd in prd_list:
        fout[prd] = None
        if not os.path.isdir(f"{out_path}/{prd}"):
            os.system(f"mkdir {out_path}/{prd}")
    with open(in_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            if line[0] == '>':
                ac = line[1:]
                for prd in fout:
                    if fout[prd] is not None:
                        if not fout[prd].closed:
                            fout[prd].close()
                    fout[prd] = open(f"{out_path}/{prd}/{ac}.caid", 'w')
                    print(f">{ac}", file = fout[prd])
                continue
            lst = line.split()
            for ii, prd in enumerate(prd_list):
                print(f"{lst[0]}\t{lst[1]}\t{lst[2 + ii]}", file = fout[prd])
        for prd in prd_list:
            if not fout[prd].closed:
                fout[prd].close()
