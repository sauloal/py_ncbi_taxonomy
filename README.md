# NCBI Taxonomy to Parquet Converter

## Introduction
Reads NCBI dump files and convert to intermediate parquet files and
single merged parquet.

## Original Source
https://github.com/sauloal/py_ncbi_taxonomy

## Prerequisites
See ```INSTALL.md```

## How to run
```python3 parse_pq.py```

You must be in the same folder as the ```dmp``` files.

## Input files
```
citations.dmp
delnodes.dmp
division.dmp
gencode.dmp
merged.dmp
names.dmp
nodes.dmp
```

## Output Files
* ```pq_division.pq```
```
division (12, 4)
division_id             int64
division_cde         category
division_name        category
division_comments      object
```


* ```pq_gencode.pq```
```
gencode (28, 5)
genetic_code_id                 int64
genetic_code_abbreviation    category
genetic_code_name            category
genetic_code_cde             category
genetic_code_starts          category
```


* ```pq_nodes.pq```
```
nodes (2174492, 13)
tax_id                              int64
parent_tax_id                       int64
rank                             category
embl_code                        category
division_id                         int64
inherited_div_flag                   bool
genetic_code_id                     int64
inherited_GC_flag                    bool
mitochondrial_genetic_code_id       int64
inherited_MGC_flag                   bool
GenBank_hidden_flag                  bool
hidden_subtree_root_flag             bool
comments                           object
```


* ```pq_names.pq```: Duplicated names merged by ```name_tax_id``` into list encoded as JSON.
```
names (2174492, 4)
name_tax_id       int64
name_txt       category
name_unique    category
name_class     category
```

With ```name_txt```, ```name_unique``` and ```name_class``` being json strings
containing the list of names.


* ```pq_all.pq```: All tables merged.
```
all_data (2174492, 27)
tax_id                              int64
parent_tax_id                       int64
rank                             category
embl_code                        category
division_id                         int64
inherited_div_flag                   bool
genetic_code_id                     int64
inherited_GC_flag                    bool
mitochondrial_genetic_code_id       int64
inherited_MGC_flag                   bool
GenBank_hidden_flag                  bool
hidden_subtree_root_flag             bool
comments                           object
name_tax_id                         int64
name_txt                         category
name_unique                      category
name_class                       category
genetic_code_abbreviation        category
genetic_code_name                category
genetic_code_cde                 category
genetic_code_starts              category
division_cde                     category
division_name                    category
division_comments                  object
is_parent                            bool
rank_id                          category
asc                                object
```

With:
- ```is_parent```: Whether it is present in ```parent_tax_id```

- ```rank_id```: Index of ```rank``` in the global order
of Ranks (dumped in ```pq_tree_ranks.json```).

- ```asc```: JSON string containing list of list of
```rank_id``` and ```parent_tax_id``` for the whole ascendence
of the tax not including itself. ```NULL``` if rank 0 (```no rank```)
or if it is a parent node (according to ```is_parent```.


### 4. Trees
#### 4.1. Ranks
```pq_tree_ranks.json``` contains the Rank names in their orders.

```Rank ID``` is the index of the rank in this list.

```json
[
 "no rank",
 "realm",
 "subrealmdomain",
 "superkingdom",
 "kingdom",
 "subkingdom",
 "infrakingdom",
 "branch",
 "superphylum",
 "phylum",
 "subphylum",
 "infraphylum",
 "microphylum",
 "superclass",
 "class",
 "subclass",
 "infraclass",
 "subterclass",
 "parvclass",
 "legion",
 "cohort",
 "subcohort",
 "magnorder",
 "superorder",
 "order",
 "suborder",
 "infraorder",
 "parvorder",
 "superfamily",
 "family",
 "subfamily",
 "infrafamily",
 "supertribe",
 "tribe",
 "subtribe",
 "infratribe",
 "genus",
 "subgenus",
 "section",
 "subsection",
 "series",
 "species group",
 "species subgroup",
 "species",
 "subspecies",
 "varietas",
 "variety",
 "forma",
 "form"
]
```

#### 4.2. Tree
Nested Dict of nodes

pq_tree.json
```json
{
 "2": {
  "tid": 2,
  "rid": 3,
  "par": null,
  "chl": {
   "1224": {
    "tid": 1224,
    "rid": 9,
    "par": 2,
    "chl": {}
   }
  }
 }
}
```
With:
* ```TID```: Tax ID
* ```RID```: Rank ID
* ```PAR```: Parent Tax ID
* ```CHL```: Children

```pq_tree.pkl``` as Pickled of the same data.


#### 4.3. Trees/
Newick trees for each root tax.

### 5. Ancestor Matrix
Instance of sc.sparse.dok_matrix containing the ancestor matrix.

```
pq_matrix.pickle
```
