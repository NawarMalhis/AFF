



def aff_load_cdhit_clusters(af, cdhit_clstr_file):
    cluster_dict = {}
    with open(cdhit_clstr_file, 'r') as fin:
        cluster = ''
        p_identity = ''
        for line in fin:
            line = line.strip()
            if len(line) < 2:
                continue
            lst = line.split()
            if lst[0][0] == '>':
                cluster = lst[1]
                continue
            ac = lst[2].split('|')[0][1:]
            if lst[3] == '*':
                cluster_dict[cluster] = ac
                p_identity = '100.00%'
            elif lst[3] == 'at':
                p_identity = lst[4]
            else:
                print(f"ERROR ===============================================\t{line}", flush=True)
            if ac not in af['data']:
                print(ac, flush=True)
                continue
            af['data'][ac]['cluster'] = int(cluster)
            af['data'][ac]['p_identity'] = p_identity
    ac_list = list(af['data'].keys())
    # print("----", len(ac_list), flush=True)
    # exit(0)
    for ac in ac_list:
        if 'cluster' not in af['data'][ac]:
            del af['data'][ac]
        # af['data'][ac]['c_center'] = cluster_dict[str(af['data'][ac]['cluster'])]
    af['metadata']['name_tags'].append('cluster')
    af['metadata']['name_tags'].append('p_identity')
    af['metadata']['name_tags'].append('c_center')
    pass
