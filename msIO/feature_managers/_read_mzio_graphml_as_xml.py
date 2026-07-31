import xml.etree.ElementTree as ET

import pandas as pd

def get_node_and_edge_data(file):
    tree = ET.parse(file)
    root = tree.getroot()

    # the first few entries define columns of the dataframe
    columns = [c.attrib for c in root]
    columns.pop()  # last element is graph
    column_types = pd.DataFrame.from_records(columns)

    graph = root[-1]

    nodes = [n for n in graph if n.tag.endswith('node')]
    edges = [e for e in graph if e.tag.endswith('edge')]

    print(f'#nodes: {len(nodes)}, #edges: {len(edges)}')
    rename_dict = {row['id']: row['attr.name'] for _, row in column_types.iterrows()}
    node_data = pd.DataFrame([{e.attrib['key']: e.text for e in n} for n in nodes]).rename(columns=rename_dict)
    edge_data = pd.DataFrame([n.attrib | {e.attrib['key']: e.text for e in n} for n in edges]).rename(columns=rename_dict)
    edge_data.index = [e.tag for e in edges]

    return node_data, edge_data


def get_modified_cosine_similarities(edge_data):
    mask_ms2_edge = edge_data.type == 'MS2 (modified) cosine'
    return (
        edge_data
        .loc[mask_ms2_edge, ['id1', 'id2', 'score']]
        .rename(columns={'id1': 'feature_id1', 'id2': 'feature_id2', 'score': 'cosine_similarity'})
    )


def collapse_ms2_similarities_to_compounds(comp_to_f_ids: dict[int, list[int]], ms2_sims: pd.DataFrame):
    """add compound id to each feature, group by compound id pairs and take maximums"""
    def find_cp_id(f_id):
        for cp_id, f_ids in comp_to_f_ids.items():
            if f_id in f_ids:
                return cp_id
        else:
            raise ValueError(f'feature id {f_id} not in compound ids')

    ms2_sims_cp = ms2_sims.copy()
    ms2_sims_cp.loc[:, 'compound_id1'] = ms2_sims_cp.feature_id1.apply(find_cp_id)
    ms2_sims_cp.loc[:, 'compound_id2'] = ms2_sims_cp.feature_id2.apply(find_cp_id)
    return (
        ms2_sims_cp
        .loc[:, [['compound_id1', 'compound_id2', 'cosine_similarity']]]
        .groupby(by=['compound_id1', 'compound_id2'])
        .max()
    )


if __name__ == '__main__':
    file = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas new method height recursive\mzmine\guaymas_mzmine_networking_new_iimn.graphml"

    node_data, edge_data = get_node_and_edge_data(file)

    sims_long = get_modified_cosine_similarities(edge_data)

    df = sims_long
