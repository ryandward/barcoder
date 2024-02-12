from elasticsearch import Elasticsearch

# Create a connection to Elasticsearch
es = Elasticsearch([{'host': 'localhost', 'port': 9200}])

# Index some example data
es.index(index='test_seqs', body={
    'sense_sequence': 'ATCG',
    'antisense_sequence': 'CGAT',
})

es.index(index='test_seqs', body={
    'sense_sequence': 'GCTA',
    'antisense_sequence': 'TAGC',
})
# Refresh the index to make sure all data is searchable
es.indices.refresh(index='test_seqs')

# The sequence to search for
query_sequence = 'ATCG'

# Perform a multi-match query
response = es.search(index='test_seqs', body={
    'query': {
        'multi_match': {
            'query': query_sequence,
            'fields': ['sense_sequence', 'antisense_sequence'],
        },
    },
})
if response:
    # Print the search results
    for hit in response['hits']['hits']:
        print(f"Found match in document {hit['_id']}:")
        print(hit['_source'])