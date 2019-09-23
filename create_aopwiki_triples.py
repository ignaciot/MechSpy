import pandas as pd
import itertools

# The tsv files can be obtained from the Downloads section of the AOPwiki:
# https://aopwiki.org/downloads

VALID_TYPES = [
                'GO', 
                'CHEBI', 
                'CL', 
                'UBERON', 
                'PR'
                ]
IGNORE_PROCESSES = [
                    'GO:0010467',
                    'GO:0046903',
                    'GO:0002790',
                    'GO:0023052',
                    'GO:0007588',
                    'GO:0010468',
                    'GO:0046983',
                    'GO:0006810',
                    'GO:0035624',
                    'GO:0005488',
                    ]

IGNORE_OBJECTS = [
                    'CHEBI:36080',
                    ]




# Create concept URI for a given type and ID
#################################################### 
def create_uri(concept_type, concept_id):
    uri = ''
    if concept_type == 'GO':
        uri = '<http://purl.obolibrary.org/obo/GO_%s>' % concept_id.split(':')[1]
    elif concept_type == 'CHEBI':
        uri = '<http://purl.obolibrary.org/obo/CHEBI_%s>' % concept_id.split(':')[1]
    elif concept_type == 'CL':
        uri = '<http://purl.obolibrary.org/obo/CL_%s>' % concept_id.split(':')[1]
    elif concept_type == 'UBERON':
        uri = '<http://purl.obolibrary.org/obo/UBERON_%s>' % concept_id.split(':')[1]
    elif concept_type == 'PR':
        uri = '<http://purl.obolibrary.org/obo/PR_%s>' % concept_id.split(':')[1]
    return uri


# Get usable nodes (if any) from the given dataframe row
############################################################# 
def get_usable_nodes(event_row):
    # if only the object or the process are given, use that
    
    # if both are given, use the process, EXCEPT if the process is:
    #   GO:0010467	gene expression
    #   GO:0046903	secretion
    #   MESH	D009154	mutation
    #   GO:0002790	peptide secretion
    #   GO:0023052	signaling
    #   GO:0007588	excretion
    #   GO:0010468	regulation of gene expression
    #   GO:0046983	protein dimerization activity
    #   GO:0006810	transport
    #   MI	MI:0195	covalent binding
    #   GO:0035624	receptor transactivation
    #   GO:0005488	binding
    #   
    # ... then use the object, since it is more relevant
    nodes = []
    if event_row.obj_type == '':
        if event_row.proc_type in VALID_TYPES and event_row.proc_id not in IGNORE_PROCESSES:
            nodes.append(create_uri(event_row.proc_type, event_row.proc_id))
        else:
            #print('skipping %s' % event_row.proc_id)
            next

    elif event_row.proc_type == '':
        if event_row.obj_type in VALID_TYPES and event_row.obj_id not in IGNORE_OBJECTS:
            nodes.append(create_uri(event_row.obj_type, event_row.obj_id))
        else:
            #print('skipping %s' % event_row.obj_id)
            next
    
    else:
        if event_row.proc_type in VALID_TYPES and event_row.proc_id not in IGNORE_PROCESSES:
            nodes.append(create_uri(event_row.proc_type, event_row.proc_id))
        elif event_row.obj_type in VALID_TYPES and event_row.obj_id not in IGNORE_OBJECTS:
                nodes.append(create_uri(event_row.obj_type, event_row.obj_id))
    return nodes



# Main stuff
####################################

# Relations dataframe
rel_df = pd.read_csv('aopwiki/aop_ke_ker.tsv', header=None, sep="\t", 
                     na_filter=False, usecols=[1, 2, 4],
                     names=['event1', 'event2', 'rel_type'])

# key event dataframe
ke_df = pd.read_csv('aopwiki/aop_ke_ec.tsv', header=None, sep="\t", 
                    na_filter=False, usecols=[1, 3, 4, 6, 7],
                    names=['event', 'obj_type', 'obj_id', 'proc_type', 'proc_id'])

with open('new_aopwiki_triples.txt', 'w') as new_fp:
    for rel in rel_df.itertuples():
        events1 = ke_df[ke_df['event'] == rel.event1].drop_duplicates()
        node_1s = []
        if len(events1) > 0:
            for idx in list(range(len(events1))):
                event_df_row = events1.iloc[idx]
                node_1s = get_usable_nodes(event_df_row)

        if len(node_1s) > 0:
            # Only bother checking concept 2 if concept 1 is something 
            # we can use to create a new edge
            events2 = ke_df[ke_df['event'] == rel.event2].drop_duplicates()
            if len(events2) > 0:
                for idx in list(range(len(events2))):
                    event_df_row = events2.iloc[idx]
                    node_2s = get_usable_nodes(event_df_row)

            if len(node_2s) > 0:
                # Only create a new edge if both related concepts are usable
                for some_node_1 in node_1s:
                    for some_node_2 in node_2s:
                        new_fp.write('%s <AOPwiki_upstream_of> %s\n' % (some_node_1, some_node_2))





