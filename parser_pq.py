#!/usr/bin/env python3

import os
import sys
import re
import csv

import simplejson as json

from collections import OrderedDict

import pandas as pd
import numpy  as np
import fastparquet
from   tqdm.auto   import tqdm

from tax import RANKS


tqdm.pandas()

SEP = r'|'


def gen_parents(nodes):
	print("generating parents")
	sys.stdout.flush()

	parents   = OrderedDict()
	
	print(" nodes", nodes.head())

	parents = nodes["parent_tax_id"].unique().tolist()
	ids     = nodes["tax_id"       ].unique().tolist()

	no_parents = list(set(ids) - set(parents))
	
	print("  total      {:12,d}".format(len(nodes     )))
	print("  parents    {:12,d}".format(len(parents   )))
	print("  no parents {:12,d}".format(len(no_parents)))

	sys.stdout.flush()

	return no_parents, parents


def trim_last_col( x ):
	return x.strip("|").strip("\t").strip().strip()


def trimmer(x):
	return x.strip("|\t ").strip("|\t ")


def dnaf(x):
	r = tuple(x.dropna())
	if len(r) == 0:
		return None
	else:
		return r


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
	# print([e for e in list(df['comments'])])
	print(pd.concat([df.head(), df.tail()]))
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
		print(pd.concat([nodes_names.head(), nodes_names.tail()]))

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
		print(pd.concat([nodes_names_gencode.head(), nodes_names_gencode.tail()]))

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
		print(pd.concat([nodes_names_gencode_division.head(), nodes_names_gencode_division.tail()]))

		print("shape  nodes_names_gencode_division")
		print(nodes_names_gencode_division.shape)

		print("dtypes nodes_names_gencode_division")
		print(nodes_names_gencode_division.dtypes)

		print("info   nodes_names_gencode_division")
		print(nodes_names_gencode_division.info())



		print(" merging nodes+names+gencode+division and parents")
		sys.stdout.flush()

		nodes_names_gencode_division_parents = gen_parents(nodes_names_gencode_division)
		print(pd.concat([nodes_names_gencode_division_parents.head(), nodes_names_gencode_division_parents.tail()]))

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
	tax_id        = node.name
	tids = str(tax_id)
	
	if tids in par_ranks:
		node_rank_id, parent_tax_id = par_ranks[tids]
	else:
		node_rank_id  = node['rank_id']
		parent_tax_id = node['parent_tax_id']
		par_ranks[tids] = (node_rank_id, parent_tax_id)

	asc               = [None]*len(RANKS)

	if tax_id == parent_tax_id:
		asc = []
	else:
		asc[node_rank_id] = tax_id
		while True:
			if node_rank_id == 0:
				break

			if tax_id == parent_tax_id:
				break

			tax_id = parent_tax_id
			tids   = str(tax_id)
			if tids in par_ranks:
				node_rank_id, parent_tax_id = par_ranks[tids]
			else:
				node          = nodes.loc[tax_id]
				node_rank_id  = node['rank_id']
				parent_tax_id = node['parent_tax_id']
				par_ranks[tids] = (node_rank_id, parent_tax_id)
			
			asc_node_rank_id = asc[node_rank_id]
			if (asc_node_rank_id is not None) and (asc_node_rank_id != tax_id):
				asc = []
				break
			else:
				asc[node_rank_id] = int(tax_id)
	
	if len(asc) == 0:
		asc = None
	else:
		asc = [(l,asc[l]) for l in range(len(asc)) if asc[l] is not None]

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

	return all_data


def gen_tree(all_data, asc, matrix, tree, rank=0):
	if rank >= len(asc):
		return

	tax_id    = asc[rank]
	if tax_id is None:
		gen_tree(all_data, asc, matrix, tree, rank=rank+1)
	
	else:
		if tax_id not in tree:
			tree[tax_id] = OrderedDict()
			# tree[tax_id].update(all_data.get(tax_id).as_dict())

		el = tree[tax_id]

		if "chl" not in el:
			el["tid"] = tax_id
			el["rid"] = rank
			el["chl"] = OrderedDict()

		gen_tree(all_data, asc, matrix, el["chl"], rank=rank+1)


def main():
	FILTER_LEVEL = RANKS.index("genus")
	FILTER_CLASS = "scientific name"
	FILTER_VAL   = "Solanum"

	all_data = get_all()

	# all_data['asc'] = all_data['tax_id'].head().map(lambda tax_id: gen_asc(all_data, tax_id))

	return


	all_data     = None
	out_data     = None
	individuals  = None
	parents      = None
	matrix       = None

	if os.path.exists("db_out_all_data.pickle") and os.path.exists("db_out_matrix.pickle"):
		print("loading db_out_all_data.pickle")
		sys.stdout.flush()
		with open("db_out_all_data.pickle", "rb") as fhd:
			all_data = pickle.load(fhd)

		print("loading db_out_matrix.pickle")
		sys.stdout.flush()
		with open("db_out_matrix.pickle", "rb") as fhd:
			matrix = pickle.load(fhd)

	else:
		print("creating all data")
		sys.stdout.flush()

		all_data    = AllData()
	
		out_data    = OrderedDict()
		individuals = OrderedDict()
		parents     = OrderedDict()

		print("summarizing data")
		sys.stdout.flush()

		par_ids = all_data.search_by_name(FILTER_VAL)
		print(" par_ids", par_ids)

		ascs   = []
		max_id = 0
		for row in all_data:
			# print(row.tax_id)

			if row.asc is None:
				# print( " no asc")
				continue

			fl = row.asc[FILTER_LEVEL]
			if fl is None:
				# print( " no fl", FILTER_LEVEL)
				continue

			if fl not in par_ids:
				continue

			print( row )

			out_data[row.tax_id] = row

			if row.is_parent:
				parents[row.tax_id] = row
			else:
				individuals[row.tax_id] = row
				ascs.append(row.asc)
				max_id = max([row.tax_id, max_id])
		
		print("  out_data    {:12,d}".format(len(out_data   )))
		print("  individuals {:12,d}".format(len(individuals)))
		print("  parents     {:12,d}".format(len(parents    )))
		print("  max_id      {:12,d}".format(max_id          ))
		sys.stdout.flush()

		# #https://www.geeksforgeeks.org/construct-tree-from-ancestor-matrix/
		# matrix = [0 for p in range(max_id * max_id)]

		# for asc in ascs:
		# 	for r in range(len(asc) - 1, 0, -1):
		# 		rv = asc[r]
		# 		if rv is None:
		# 			continue
		# 		for l in range(r - 1):
		# 			lv = asc[l]
		# 			if lv is None:
		# 				continue
		# 			matrix[rv * lv] = 1

		print("saving db_out_all_data.pickle")
		sys.stdout.flush()
		with open("db_out_all_data.pickle", "wb") as fhd:
			pickle.dump(all_data, fhd, protocol=0)

		print("saving db_out_rows.pickle")
		sys.stdout.flush()
		with open("db_out_rows.pickle", "wb") as fhd:
			pickle.dump(out_data, fhd, protocol=0)

		print("saving db_out_rows_parents.pickle")
		sys.stdout.flush()
		with open("db_out_rows_parents.pickle", "wb") as fhd:
			pickle.dump(parents, fhd, protocol=0)

		print("saving db_out_rows_individuals.pickle")
		sys.stdout.flush()
		with open("db_out_rows_individuals.pickle", "wb") as fhd:
			pickle.dump(individuals, fhd, protocol=0)

		print("saving db_out_matrix.pickle")
		sys.stdout.flush()
		with open("db_out_matrix.pickle", "wb") as fhd:
			pickle.dump(matrix, fhd, protocol=0)

			


	print("creating tree")
	sys.stdout.flush()
	tree        = OrderedDict()
	rows        = OrderedDict()
	for row_num, row in enumerate(all_data):
		# print(tax_id, row.tax_id])
		rows[row.tax_id] = row

		if row.asc:
			# print(row.tax_id, row.asc)
			gen_tree(all_data, row.asc, matrix, tree)

		# if row_num == 100:
		# 	break
	
	print("saving json tree")
	sys.stdout.flush()
	# print(json.dumps(tree, indent=1, for_json=True))
	with open("db_tree.json", "wt") as fhd:
		json.dump(tree, fhd, indent=1, for_json=True)

	print("saving json data")
	sys.stdout.flush()
	# print(json.dumps(tree, indent=1, for_json=True))
	with open("db_tree_data.json", "wt") as fhd:
		json.dump(rows, fhd, indent=1, for_json=True)

	print("saving json ranks")
	sys.stdout.flush()
	# print(json.dumps(tree, indent=1, for_json=True))
	with open("db_tree_ranks.json", "wt") as fhd:
		json.dump(RANKS, fhd, indent=1, for_json=True)

RANKS

if __name__ == "__main__":
	main()
