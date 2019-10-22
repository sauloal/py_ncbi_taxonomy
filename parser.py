#!/usr/bin/env python3

import os
import sys
import pickle

from collections import OrderedDict

import simplejson as json
#import json

from tax import *

import pandas as pd
import fastparquet


class AllData(General):
	__slots__      = [
		"divisions",
		"genetic_codes",
		"names",
		"nodes",
		# "no_parents",
		# "parents"
	]
	__attributes__ = __slots__
	__shared__     = []

	def __init__(self):
		self.divisions                = parse_division()
		self.genetic_codes            = parse_gencode()
		self.nodes                    = parse_nodes()
		self.names                    = parse_names()
		self.no_parents, self.parents = gen_parents(self.nodes)

	def __iter__(self):
		return self.next()

	def __next__(self):
		return self.next()

	def get(self, tax_id):
		if tax_id not in self.nodes:
			raise KeyError("No such tax_id")
		
		node = self.nodes[tax_id]

		(
			tax_id,
			parent_tax_id,
			rank_id,
			rank_name,
			embl_code_id,
			embl_code,
			division_id,
			inherited_div_flag,
			genetic_code_id,
			inherited_GC_flag,
			mitochondrial_genetic_code_id,
			inherited_MGC_flag,
			genBank_hidden_flag,
			hidden_subtree_root_flag
		) = node.get()

		parent         = self.nodes[parent_tax_id]
		parent_rank_id = parent.rank_id
		parent_rank    = parent.rank_name

		if rank_id != 0:
			assert rank_id >= parent_rank_id, "rank {} ({}) < than rank {} ({})".format(rank_id, rank_name, parent_rank_id, parent_rank)

		name_tax_id, names         = self.names[tax_id].get()
		division                   = self.divisions[division_id].get()
		genetic_code               = self.genetic_codes[genetic_code_id]
		mitochondrial_genetic_code = self.genetic_codes[mitochondrial_genetic_code_id]

		_, division_cde_id, division_cde, division_name_id, division_name, division_comments = division

		genetic_code_abbreviation              , genetic_code_name              , genetic_code_cde              , genetic_code_starts               = genetic_code
		mitochondrial_genetic_code_abbreviation, mitochondrial_genetic_code_name, mitochondrial_genetic_code_cde, mitochondrial_genetic_code_starts = mitochondrial_genetic_code

		is_parent = tax_id in self.parents
		asc       = None
		if True: #not is_parent:
			asc  = [None] * len(RANKS)
			try:
				gen_asc(self.nodes, asc, tax_id)
			except AssertionError as e:
				print(e)
				print(rank_id, rank_name, division_cde, division_name, names)
				asc = None

		row = Row(
			tax_id,
			parent_tax_id,

			rank_id,
			rank_name,

			embl_code_id,
			embl_code,
			
			division_id,
			division_cde_id,
			division_cde,
			division_name_id,
			division_name,
			division_comments,
			inherited_div_flag,

			genetic_code_id,
			genetic_code_abbreviation,
			genetic_code_name,
			genetic_code_cde,
			genetic_code_starts,
			inherited_GC_flag,

			mitochondrial_genetic_code_id,
			mitochondrial_genetic_code_abbreviation,
			mitochondrial_genetic_code_name,
			mitochondrial_genetic_code_cde,
			mitochondrial_genetic_code_starts,
			inherited_MGC_flag,

			names,

			genBank_hidden_flag,
			hidden_subtree_root_flag,

			is_parent,
			asc
		)

		return row

	def next(self):
		for node_pos, tax_id in enumerate(self.nodes.keys()):
			if (node_pos+1) % 100000 == 0:
				print(" {:12,d} / {:12,d}".format(node_pos+1, len(self.nodes)))
				sys.stdout.flush()

			yield self.get(tax_id)

	def search_by_name(self, FILTER_VAL):
		nids = []
		for nid, name in self.names.items():
			if FILTER_VAL.lower() in [t.lower() for t in name.txts]:
				nids.append(nid)
		return list(set(nids))


def parser_save(fn, cls, data):
	# for e, (k,v) in enumerate(data.items()):
	# 	if e < 5 or e > len(data) - 5:
	# 		print(k, v)

	print("  converting to dataframe")
	sys.stdout.flush()

	df = pd.DataFrame.from_dict({k: v._to_pandas() for k,v in data.items()}, orient="index", columns=cls.__pds__)
	
	# print(pd.concat([df.head(), df.tail()]))
	
	print("  saving to", fn)
	sys.stdout.flush()

	df.to_parquet(fn, engine="fastparquet", compression="snappy", index=True)#, partition_cols=)
	
	print("  done")
	sys.stdout.flush()

	return df


def parser_load(fn, cls):
	print("  loading", fn)
	sys.stdout.flush()

	data = pd.read_parquet(fn, engine="fastparquet")

	print("  converting to class")
	sys.stdout.flush()

	# print(pd.concat([data.head(), data.tail()]))

	data = data.to_dict(orient="index", into=OrderedDict)
	data = OrderedDict((k, cls(*[v[s] for s in cls.__pds__])) for k,v in data.items())

	# for e, (k,v) in enumerate(data.items()):
	# 	if e < 5 or e > len(data) - 5:
	# 		print(k, v)

	print("  done")
	sys.stdout.flush()

	return data


def parse_division():
	print("parsing division")

	if os.path.exists("db_division.pq"):
		print(" loading db_division.pq")
		sys.stdout.flush()
		
		divisions = parser_load(fn="db_division.pq", cls=Division)

		# with open("db_division.pickle", "rb") as fhd:
		# 	divisions = pickle.load(fhd)

	else:
		divisions        = OrderedDict()
		
		for linenum, cols in read_dump("division.dmp"):
			# print(linenum, cols)

			division = Division(*cols)

			assert division.did not in divisions

			divisions[division.did] = division

			# print(division)
		
		print(" saving db_division.pq")
		sys.stdout.flush()
		
		divisions = parser_save(fn="db_division.pq", cls=Division, data=divisions)
		# divisions2 = parser_load(fn="db_division.pq", cls=Division)
		# assert repr(divisions2) == repr(divisions)

		# with open("db_division.pickle", "wb") as fhd:
		# 	pickle.dump(divisions, fhd, protocol=0)
		
	print("  divisions  {:12,d}".format(len(divisions     )))
	sys.stdout.flush()

	return divisions


def parse_gencode():
	print("parsing gencode")

	if os.path.exists("db_gencode.pickle"):
		print(" loading db_gencode.pickle")
		sys.stdout.flush()
		
		with open("db_gencode.pickle", "rb") as fhd:
			gencodes = pickle.load(fhd)

	else:
		gencodes = OrderedDict()
		
		for linenum, cols in read_dump("gencode.dmp"):
			# print(linenum, cols)

			genetic_code_id, abbreviation, name, cde, starts = cols
			genetic_code_id = int(genetic_code_id)

			gencodes[genetic_code_id] = (abbreviation, name, cde, starts)
		
		print(" saving db_gencode.pickle")
		sys.stdout.flush()
		
		with open("db_gencode.pickle", "wb") as fhd:
			pickle.dump(gencodes, fhd, protocol=0)
		
	print("  gencodes   {:12,d}".format(len(gencodes)))
	sys.stdout.flush()

	return gencodes


def parse_nodes():
	print("parsing nodes")

	if os.path.exists("db_nodes.pq"):
		print(" loading db_nodes.pq")
		sys.stdout.flush()
		
		# with open("db_nodes.pickle", "rb") as fhd:
		# 	nodes = pickle.load(fhd)

		nodes = parser_load(fn="db_nodes.pq", cls=Node)

	else:
		nodes   = OrderedDict()
		
		for linenum, cols in read_dump("nodes.dmp"):
			# print(linenum, cols)
			if (linenum+1) % 100000 == 0:
				print(" {:12,d}".format(linenum+1))
				sys.stdout.flush()

			node = Node(*cols[:12])
			assert node.tax_id not in nodes 
			nodes[node.tax_id] = node
			# print(node)
		
		print(" saving db_nodes.pq")
		sys.stdout.flush()
		
		parser_save(         fn="db_nodes.pq", cls=Node, data=nodes)
		# nodes2 = parser_load(fn="db_nodes.pq", cls=Node)
		# assert repr(nodes2) == repr(nodes)

		# with open("db_nodes.pickle", "wb") as fhd:
		# 	pickle.dump(nodes, fhd, protocol=0)
		
	print("  nodes      {:12,d}".format(len(nodes     )))
	sys.stdout.flush()

	return nodes


def parse_names():
	print("parsing names")

	if os.path.exists("db_names.pq"):
		print(" loading db_names.pq")
		sys.stdout.flush()
		
		names = parser_load(fn="db_names.pq", cls=Name)

		# with open("db_names.pickle", "rb") as fhd:
		# 	names = pickle.load(fhd)

	else:
		names          = OrderedDict()
		
		for linenum, cols in read_dump("names.dmp"):
			# print(linenum, cols)
			if (linenum+1) % 100000 == 0:
				print(" {:12,d}".format(linenum+1))
				sys.stdout.flush()

			name = Name(*cols)

			if name.tax_id in names:
				names[name.tax_id].merge(name)
				# assert tax_id not in names, "duplicated tax_id {}. {}".format(tax_id, names[tax_id])
			else:
				names[name.tax_id] = name
			# print(names[tax_id])
			# print()

		print(" saving db_names.pq")
		sys.stdout.flush()

		parser_save(         fn="db_names.pq", cls=Name, data=names)
		# names2 = parser_load(fn="db_names.pq", cls=Name)
		# assert repr(names2) == repr(names)

		# with open("db_names.pickle", "wb") as fhd:
		# 	pickle.dump(names, fhd, protocol=0)

	print("  names      {:12,d}".format(len(names       )))

	return names


def gen_parents(nodes):
       print("generating parents")
       sys.stdout.flush()

       parents   = OrderedDict()

       for pos, (tax_id, node) in enumerate(nodes.items()):
               (
                       tax_id,
                       parent_tax_id,
                       rank_id,
                       rank_name,
                       embl_code_id,
                       embl_code,
                       division_id,
                       inherited_div_flag,
                       genetic_code_id,
                       inherited_GC_flag,
                       mitochondrial_genetic_code_id,
                       inherited_MGC_flag,
                       genBank_hidden_flag,
                       hidden_subtree_root_flag
               ) = node.get()

               if (pos + 1) % 100000 == 0:
                       print(" {:12,d}".format(pos+1))
                       sys.stdout.flush()

               parents[parent_tax_id] = parents.get(parent_tax_id, 0) + 1

       no_parents = list(set(nodes.keys()) - set(parents.keys()))

       print("  total      {:12,d}".format(len(nodes     )))
       print("  parents    {:12,d}".format(len(parents   )))
       print("  no parents {:12,d}".format(len(no_parents)))

       sys.stdout.flush()

       return no_parents, parents


def gen_asc(nodes, asc, tax_id, level=0):
	node          = nodes[tax_id]
	node_rank_id  = node.rank_id
	parent_tax_id = node.parent_tax_id

	assert asc[node_rank_id] is None and asc[node_rank_id] != tax_id, (tax_id, node_rank_id, RANKS[node_rank_id], parent_tax_id, level, asc)

	if node_rank_id != 0:
		asc[node_rank_id] = tax_id
		gen_asc(nodes, asc, parent_tax_id, level=level+1)


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

			# print( row )

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
