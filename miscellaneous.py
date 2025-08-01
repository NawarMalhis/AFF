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
    response = requests.get(url, headers={"Accept": "application/json"})

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

