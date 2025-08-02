import requests


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