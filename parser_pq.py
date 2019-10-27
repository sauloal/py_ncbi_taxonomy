#!/usr/bin/env python3

import os
import sys
import re
import csv
import pickle

from collections import OrderedDict

import simplejson as json

import pandas as pd
import numpy  as np
import scipy  as sc
import fastparquet
from   tqdm.auto       import tqdm
from   IPython.display import display

from tax import RANKS

tqdm.pandas()

__version__ = '191027.2130'

SEP = r'|'

print("parser PQ")


def gen_parents(nodes):
	print("generating parents")
	sys.stdout.flush()

	parents   = OrderedDict()
	
	print(" nodes")
	display_df(nodes)

	parents = nodes["parent_tax_id"].unique().tolist()
	ids     = nodes["tax_id"       ].unique().tolist()

	no_parents = list(set(ids) - set(parents))
	
	print("  total      {:12,d}".format(len(nodes     )))
	print("  parents    {:12,d}".format(len(parents   )))
	print("  no parents {:12,d}".format(len(no_parents)))

	sys.stdout.flush()

	return no_parents, parents


def trim_last_col(x):
	return x.strip("|").strip("\t").strip().strip()


def trimmer(x):
	return x.strip("|\t ").strip("|\t ")


def dnaf(x):
	r = tuple(x.dropna())
	if len(r) == 0:
		return None
	else:
		return r


def display_df(df, window=2):
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        if len(df) <= (2*window):
            display(df)
        else:
            display(pd.concat([df.head(window), df.tail(window)]))


def gen_parents(all_data):
	print("generating parents")
	sys.stdout.flush()

	parents    = all_data["parent_tax_id"].unique().tolist()
	ids        = all_data.index.values.tolist()

	no_parents = list(set(ids) - set(parents))

	print("  total       {:12,d}".format(len(all_data  )))
	print("  parents     {:12,d}".format(len(parents   )))
	print("  not parents {:12,d}".format(len(no_parents)))
	sys.stdout.flush()

	all_data['is_parent'] = all_data.index.isin(parents)

	print(sum(all_data['is_parent']))

	return all_data


def parse_dump(name, ifn, ofn, cols, converters, index_col=0, sep=SEP, post=None, save=True):
	print("parsing", name)
	sys.stdout.flush()

	if os.path.exists(ofn):
		print(" loading", ofn)
		sys.stdout.flush()
		
		df = pd.read_parquet(ofn, engine="fastparquet")

	else:
		print(" reading", ifn)
		sys.stdout.flush()

		df = pd.read_csv(
			ifn,
			engine='c',
			sep=SEP,
			header=None,
			names=cols.keys(),
			index_col=False,
			# index_col=index_col,
			# dtype=cols,
			skipinitialspace=True,
			converters=converters
		)

		for col, dtype in cols.items():
			df[col] = df[col].astype(dtype, copy=False)

		if post is not None:
			print(" running post processing")
			sys.stdout.flush()

			df = post(df)

		if save:
			print(" saving", ofn)
			sys.stdout.flush()

			df.to_parquet(ofn, engine="fastparquet", compression="snappy", index=True)#, partition_cols=)
	
	print("  ", name, df.shape)
	print("    dtypes\n", df.dtypes)
	print("    info\n", df.info())
	display_df(df)
	sys.stdout.flush()

	return df


def parse_division(save=True):
	name = "division"
	ifn  = "division.dmp"
	ofn  = "pq_division.pq"

	cols = OrderedDict((
		("division_id"      , np.int    ), # -- taxonomy database division id
		("division_cde"     , 'category'), # -- GenBank division code (three characters)
		("division_name"    , 'category'), # -- e.g. BCT, PLN, VRT, MAM, PRI...
		("division_comments", np.str    )  #
	))

	converters = {k: trimmer for k in cols.keys()}

	df = parse_dump(name, ifn, ofn, cols, converters, save=save)

	return df


def parse_gencode(save=True):
	name = "gencode"
	ifn  = "gencode.dmp"
	ofn  = "pq_gencode.pq"

	cols = OrderedDict((
		("genetic_code_id"          , np.int    ), # -- GenBank genetic code id
		("genetic_code_abbreviation", 'category'), # -- genetic code name abbreviation
		("genetic_code_name"        , 'category'), # -- genetic code name
		("genetic_code_cde"         , 'category'), # -- translation table for this genetic code
		("genetic_code_starts"      , 'category')  # - start codons for this genetic code
	))

	converters = {k: trimmer for k in cols.keys()}

	df = parse_dump(name, ifn, ofn, cols, converters, save=save)

	return df


def parse_nodes(save=True):
	name = "nodes"
	ifn  = "nodes.dmp"
	ofn  = "pq_nodes.pq"

	cols = OrderedDict((
		("tax_id"                       , np.int    ), #          -- node id in GenBank taxonomy database
		("parent_tax_id"                , np.int    ), #          -- parent node id in GenBank taxonomy database
		("rank"                         , 'category'), #          -- rank of this node (superkingdom, kingdom, ...) 
		("embl_code"                    , 'category'), #          -- locus-name prefix; not unique
		("division_id"                  , np.int    ), #          -- see division.dmp file
		("inherited_div_flag"           , np.bool   ), # (1 or 0) -- 1 if node inherits division from parent
		("genetic_code_id"              , np.int    ), #          -- see gencode.dmp file
		("inherited_GC_flag"            , np.bool   ), # (1 or 0) -- 1 if node inherits genetic code from parent
		("mitochondrial_genetic_code_id", np.int    ), #          -- see gencode.dmp file
		("inherited_MGC_flag"           , np.bool   ), # (1 or 0) -- 1 if node inherits mitochondrial gencode from parent
		("GenBank_hidden_flag"          , np.bool   ), # (1 or 0) -- 1 if name is suppressed in GenBank entry lineage
		("hidden_subtree_root_flag"     , np.bool   ), # (1 or 0) -- 1 if this subtree has no sequence data yet
		("comments"                     , np.str    ), #          -- free-text comments and citations
	))

	converters = {k: trimmer for k in cols.keys()}

	df = parse_dump(name, ifn, ofn, cols, converters, save=save)

	return df


def parse_names(save=True):
	name = "names"
	ifn  = "names.dmp"
	ofn  = "pq_names.pq"

	cols = OrderedDict((
		("name_tax_id", np.int    ), # -- the id of node associated with this name
		("name_txt"   , 'category'), # -- name itself
		("name_unique", 'category'), # -- the unique variant of this name if name not unique
		("name_class" , 'category')  # -- (synonym, common name, ...)
	))

	cols_to_group = list(cols.keys())[1:]


	if os.path.exists(ofn):
		print(" loading", ofn)
		sys.stdout.flush()
		
		df = pd.read_parquet(ofn, engine="fastparquet")

	else:
		print(" loading", ifn)
		sys.stdout.flush()

		grp = OrderedDict()

		with open(ifn, "rt", newline='') as csvfile:
			spamreader = csv.reader(csvfile, delimiter='|')
			for r, row in enumerate(spamreader):
				if (r+1) % 100000 == 0:
					print("{:12,d}".format(r+1))
					sys.stdout.flush()

				row    = [trimmer(r) for r in row]
				row[0] = int(row[0])
				group  = row[0]

				if group not in grp:
					grp[group] = [[] for _ in cols_to_group]

				for p in range(len(cols_to_group)):
					v = row[p + 1]
					if len(v) > 0:
						grp[group][p].append(str(v))

		print(" cleaning empty")
		sys.stdout.flush()
		for data in grp.values():
			for p in range(len(cols_to_group)):
				if len(data[p]) == 0:
					data[p] = None

		print(" converting to dataframe")
		sys.stdout.flush()

		df = pd.DataFrame.from_dict(grp, orient='index', columns=cols_to_group)

		print(" converting type")
		sys.stdout.flush()

		for k, v in list(cols.items())[1:]:
			if v == 'category':
				df[k] = df[k].apply(lambda x: json.dumps(x)).astype('category', copy=False)
			else:
				df[k] = df[k].astype(v, copy=False)

		df.index.name = list(cols.keys())[0]
		df = df.reset_index()

		print(df)
		print(df.dtypes)
		print(df.info())

		if save:
			print(" saving", ofn)
			sys.stdout.flush()

			df.to_parquet(ofn, engine="fastparquet", compression="snappy", index=True)#, partition_cols=)

	print("  ", name, df.shape)
	print("    dtypes\n", df.dtypes)
	print("    info\n", df.info())
	display_df(df)

	return df


def get_all(save=True, make_asc=True):
	ofn = "pq_all.pq"

	all_data = None

	if os.path.exists(ofn):
		print(" reading", ofn)
		sys.stdout.flush()

		all_data = pd.read_parquet(ofn, engine="fastparquet")
		
	else:
		print(" reading intermediate databases")
		sys.stdout.flush()

		divisions = parse_division(save=save)
		gencode   = parse_gencode(save=save)
		nodes     = parse_nodes(save=save)
		names     = parse_names(save=save)



		print(" merging nodes and names")
		sys.stdout.flush()

		nodes_names = nodes.merge(names, how='left', left_on='tax_id', right_on='name_tax_id', copy=False)
		# nodes_names.index.name = 'tax_id'
		display_df(nodes_names)

		print("  shape  nodes_names")
		print(nodes_names.shape)

		print("  dtypes nodes_names")
		print(nodes_names.dtypes)

		print(" info   nodes_names")
		print(nodes_names.info())



		print(" merging nodes+names and gencode")
		sys.stdout.flush()

		nodes_names_gencode = nodes_names.merge(gencode, how='left', left_on='genetic_code_id', right_on='genetic_code_id', copy=False)
		# nodes_names_gencode.index.name = 'tax_id'
		display_df(nodes_names_gencode)

		print("  shape  nodes_names_gencode")
		print(nodes_names_gencode.shape)

		print("  dtypes nodes_names_gencode")
		print(nodes_names_gencode.dtypes)

		print("  info   nodes_names_gencode")
		print(nodes_names_gencode.info())



		print(" merging nodes+names+gencode and division")
		sys.stdout.flush()

		nodes_names_gencode_division = nodes_names_gencode.merge(divisions, how='left', left_on='division_id', right_on='division_id', copy=False)
		# nodes_names_gencode_division.index.name = 'tax_id'
		display_df(nodes_names_gencode_division)

		print("shape  nodes_names_gencode_division")
		print(nodes_names_gencode_division.shape)

		print("dtypes nodes_names_gencode_division")
		print(nodes_names_gencode_division.dtypes)

		print("info   nodes_names_gencode_division")
		print(nodes_names_gencode_division.info())



		print(" merging nodes+names+gencode+division and parents")
		sys.stdout.flush()

		nodes_names_gencode_division_parents = gen_parents(nodes_names_gencode_division)
		display_df(nodes_names_gencode_division_parents)

		print("  shape  nodes_names_gencode_division_parents")
		print(nodes_names_gencode_division_parents.shape)

		print("  dtypes nodes_names_gencode_division_parents")
		print(nodes_names_gencode_division_parents.dtypes)

		print("  info   nodes_names_gencode_division_parents")
		print(nodes_names_gencode_division_parents.info())


		print("  adding rank id")
		sys.stdout.flush()

		nodes_names_gencode_division_parents["rank_id"] = nodes_names_gencode_division_parents["rank"].apply(lambda x: RANKS.index(x))

		all_data = nodes_names_gencode_division_parents

		if make_asc:
			all_data = gen_asc(all_data)

		if save:
			print(" saving merged data", ofn)
			all_data.to_parquet(ofn, engine="fastparquet", compression="snappy", index=True)#, partition_cols=)

	return all_data


def gen_asc_apply(nodes, node, par_ranks):
	tax_id = node.name
	tids   = str(tax_id)
	
	if tids in par_ranks:
		tax_id, node_rank_id, parent_tax_id = par_ranks[tids]
	else:
		node_rank_id    = node['rank_id'      ]
		parent_tax_id   = node['parent_tax_id']
		par_ranks[tids] = (tax_id, node_rank_id, parent_tax_id)

	asc = [None]*len(RANKS)

	if tax_id == parent_tax_id:
		asc = []
	
	else:
		orig_tax_id       = int(tax_id)
		# asc[node_rank_id] = int(tax_id)
		
		while True:
			if node_rank_id == 0:
				break

			if tax_id == parent_tax_id:
				break

			tax_id = parent_tax_id
			tids   = str(tax_id)
			if tids in par_ranks:
				tax_id, node_rank_id, parent_tax_id = par_ranks[tids]
			else:
				node            = nodes.loc[tax_id]
				node_rank_id    = node['rank_id'      ]
				parent_tax_id   = node['parent_tax_id']
				par_ranks[tids] = (tax_id, node_rank_id, parent_tax_id)

			if node_rank_id == 0:
				break

			if tax_id == parent_tax_id:
				break

			asc_node_rank_id = asc[node_rank_id]
			if asc_node_rank_id is not None:
				if asc_node_rank_id != tax_id:
					asc = []
					break
			else:
				asc[node_rank_id] = int(tax_id)

	if len(asc) == 0:
		asc = None
	else:
		asc = [(l,asc[l]) for l in range(len(asc)) if asc[l] is not None]
		if len(asc) == 0:
			asc = None

	return json.dumps(asc)


def gen_asc(all_data):
	try:
		all_data.reset_index(inplace=True)
	except:
		pass

	if 'asc' in all_data.columns:
		all_data.drop('asc', inplace=True, axis=1)
	if 'level_0' in all_data.columns:
		all_data.drop('level_0', inplace=True, axis=1)
	if 'index' in all_data.columns:
		all_data.drop('index', inplace=True, axis=1)
		
	all_data.set_index('tax_id', verify_integrity=True, inplace=True)

	par_ranks = {}
	all_data['asc'] = all_data[['rank_id', 'parent_tax_id']].progress_apply(lambda node: gen_asc_apply(all_data, node, par_ranks), axis=1)

	all_data.reset_index(inplace=True)

	display_df(all_data)

	print("  shape  all_data")
	print(all_data.shape)

	print("  dtypes all_data")
	print(all_data.dtypes)

	print("  info   all_data")
	print(all_data.info())

	return all_data


def gen_tree_sub(asc, tree):
    #print('asc', asc, 'tree', tree)
    h = tree
    
    for pos, (rank_id, tax_id) in enumerate(asc):
        par = None
        if pos == 0:
            par = None
        else:
            par = asc[pos - 1][1]

        if tax_id not in tree:
            h[tax_id] = OrderedDict()
        
        el = h[tax_id]
        
        if "chl" not in el:
            el["tid"] = tax_id
            el["rid"] = rank_id
            el["par"] = par
            el["chl"] = OrderedDict()
            h         = el["chl"]
        else:
            h         = el["chl"]
            
    return tree


def gen_tree(all_data):
	tree   = OrderedDict()

	not_parent          = ~all_data['is_parent']
	not_null            = all_data['asc'] != 'null'
	not_parent_not_null = not_parent & not_null
	print("all_data           ", all_data.shape)
	print("not_parent         ", not_parent.sum())
	print("not_null           ", not_null.sum())
	print("not_parent_not_null", not_parent_not_null.sum())

	with tqdm(total=not_parent_not_null.sum()) as pbar:
		for row_id, row in all_data[not_parent_not_null].iterrows():
			pbar.update()
		#     print(row_id, row)
			asc = json.loads(row['asc'])

			if asc is not None:
				tree = gen_tree_sub(asc, tree)

		#     if row_id >= 1000:
		#         break

		# print(json.dumps(tree, indent=1))

	pkl       = 'pq_trees.pkl'
	print(" saving tree as pickle", pkl)
	sys.stdout.flush()

	with open(pkl, 'wb') as fhd:
		pickle.dump(tree, fhd, protocol=0)

	print(" saving tree as json")
	sys.stdout.flush()
	with open("pq_tree.json", "wt") as fhd:
		json.dump(tree, fhd, indent=1, for_json=True)

	return tree


def tree_to_newick(tree):
    newick = []
    for tax_id, data in tree.items():
        chl, l = tree_to_newick(data['chl'])
        # print("tax_id", tax_id, " chl", chl, "l", l, "data", data)
        if   l == 0:
            newick.append("s{}".format(tax_id))
        elif l == 1:
            newick.append("({})s{}".format(chl, tax_id))
        else:
            newick.append("({})s{}".format(chl, tax_id))

    if len(newick) == 0:
        return "", 0
    if len(newick) == 1:
        return "{}".format(newick[0]), len(newick)
    else:
        return "({})".format(",".join(newick)), len(newick)


def dump_tree_as_newick(all_data, tree):
	with tqdm(total=len(tree)) as pbar:
		for root_parent_tax_id, child_tree in tree.items():
			pbar.update()
			root_parent   = all_data[all_data['tax_id'] == root_parent_tax_id]
			division_cde  = root_parent['division_cde' ].tolist()[0]
			division_name = root_parent['division_name'].tolist()[0]
			rank          = root_parent['rank'         ].tolist()[0]
			name          = root_parent['name_txt'     ].tolist()[0]

			name = json.loads(name)[0]
		    # print(root_parent)
		    # print(child_tree)
			# print(root_parent_tax_id, division_cde, division_name, name)
			sys.stdout.flush()

			if not os.path.exists('trees'):
				os.makedirs('trees')

			bn = '{:09d}_{}_{}_{}_{}'.format(root_parent_tax_id, division_cde, division_name, rank, name)
			bn = "".join(["_" if b in "\\/()[]{}.-'\" " else b for b in bn])
			bn = bn.replace("__", "_").replace("__", "_").replace("__", "_")
			bn = os.path.join('trees', bn)
			# print(" basename", bn)

			# print(" creating newick")
			sys.stdout.flush()
			newick, _ = tree_to_newick({root_parent_tax_id: child_tree})
			nwk       = '{}.newick'.format(bn)

			# print(" saving newick", nwk)
			sys.stdout.flush()
			with open(nwk, 'wt') as fhd:
				fhd.write("(")
				fhd.write(newick)
				fhd.write(")root;")


def gen_matrix(all_data):
	# https://www.geeksforgeeks.org/construct-tree-from-ancestor-matrix/

	sparse_matrix = sc.sparse.dok_matrix
	len_all_data  = all_data['tax_id'].max()
	print("len_all_data", len_all_data)
	sm            = sparse_matrix((len_all_data,len_all_data), dtype=np.int8)
	print("sm", repr(sm))

	not_parent          = ~all_data['is_parent']
	not_null            = all_data['asc'] != 'null'
	not_parent_not_null = not_parent & not_null

	print("all_data           ", all_data.shape)
	print("not_parent         ", not_parent.sum())
	print("not_null           ", not_null.sum())
	print("not_parent_not_null", not_parent_not_null.sum())

	with tqdm(total=not_parent_not_null.sum()) as pbar:
		for row_id, row in all_data[not_parent_not_null].iterrows():
			pbar.update()
			asc = json.loads(row['asc'])

			if asc is not None and len(asc) > 1:
				# print("asc", asc)
				for a in range(len(asc)-1):
					a_lvl, a_id = asc[a]
					# print("  a_id", a_id)
					for b in range(a+1, len(asc)):
						b_lvl, b_id = asc[b]
						# print("   b_id", b_id)
						sm[a_id, b_id] = 1
	
	print("sm", repr(sm))
	# print("sm nnz", sm.nnz)
	# print("sm dok", sm.todok())

	print("saving pq_matrix.pkl")
	sys.stdout.flush()
	with open("pq_matrix.pkl", "wb") as fhd:
		pickle.dump(sm, fhd, protocol=0)

	return sm


def main():
	FILTER_LEVEL = RANKS.index("genus")
	FILTER_CLASS = "scientific name"
	FILTER_VAL   = "Solanum"

	all_data     = get_all(save=True, make_asc=True)

	print("ALL DATA")
	print(all_data.shape)
	display_df(all_data)

	print("HAS ASC")
	print(all_data[all_data['asc']=='null'].shape)
	display_df(all_data[all_data['asc']=='null'])

	print("NO ASC")
	print(all_data[all_data['asc']!='null'].shape)
	display_df(all_data[all_data['asc']!='null'])

	print("NO RANK")
	print(all_data[all_data['rank']=='no rank'].shape)
	display_df(all_data[all_data['rank']=='no rank'])

	print("REALM")
	print(all_data[all_data['rank']=='realm'].shape)
	display_df(all_data[all_data['rank']=='realm'])

	print("SUB REALM")
	print(all_data[all_data['rank']=='subrealm'].shape)
	display_df(all_data[all_data['rank']=='subrealm'])

	print("DOMAIN")
	print(all_data[all_data['rank']=='domain'].shape)
	display_df(all_data[all_data['rank']=='domain'])

	print("SUPER KINGDOMS")
	print(all_data[all_data['rank']=='superkingdom'].shape)
	display_df(all_data[all_data['rank']=='superkingdom'])

	print("KINGDOMS")
	print(all_data[all_data['rank']=='kingdom'].shape)
	display_df(all_data[all_data['rank']=='kingdom'])

	print(" CREATING TREE")
	sys.stdout.flush()
	tree      = gen_tree(all_data)

	print(" SAVING TREES AS NEWICK")
	sys.stdout.flush()
	dump_tree_as_newick(all_data, tree)

	print(" CREATING ANCESTOR MATRIX")
	sys.stdout.flush()
	ancestor_matrix = gen_matrix(all_data)

	print("SAVING JSON RANKS")
	sys.stdout.flush()
	with open("pq_tree_ranks.json", "wt") as fhd:
		json.dump(RANKS, fhd, indent=1, for_json=True)


if __name__ == "__main__":
	main()
