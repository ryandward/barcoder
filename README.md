A suite of tools for pooled barcoded genomic DNA experiments
# Barcoder Tools
## Development note:
* The structure of the locus map powers most of the functions. That structure is as follows:
```json
{
    "<id>": {
        "<position>": [
            {
                "locus_tag": "<locus_tag>",
                "gene_name": "<gene_name>",
                "start": "<start>",
                "end": "<end>",
                "strand": "<strand>"
            },
            ...
        ],
        ...
    },
    ...
}
```
