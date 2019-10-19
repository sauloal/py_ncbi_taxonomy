#!/usr/bin/env python3

import os
import sys
import re
import pickle

import simplejson as json
#import json

from collections import OrderedDict

import pandas as pd
import fastparquet

"""
Each record consists of one or more fields delimited by "\t|\t"

nodes.dmp
---------

This file represents taxonomy nodes. The description for each node includes 
the following fields:

	tax_id					-- node id in GenBank taxonomy database
 	parent tax_id				-- parent node id in GenBank taxonomy database
 	rank					-- rank of this node (superkingdom, kingdom, ...) 
 	embl code				-- locus-name prefix; not unique
 	division id				-- see division.dmp file
 	inherited div flag  (1 or 0)		-- 1 if node inherits division from parent
 	genetic code id				-- see gencode.dmp file
 	inherited GC  flag  (1 or 0)		-- 1 if node inherits genetic code from parent
 	mitochondrial genetic code id		-- see gencode.dmp file
 	inherited MGC flag  (1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
 	GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
 	hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
 	comments				-- free-text comments and citations

names.dmp
---------
Taxonomy names file has these fields:

	tax_id					-- the id of node associated with this name
	name_txt				-- name itself
	unique name				-- the unique variant of this name if name not unique
	name class				-- (synonym, common name, ...)

division.dmp
------------
Divisions file has these fields:
	division id				-- taxonomy database division id
	division cde				-- GenBank division code (three characters)
	division name				-- e.g. BCT, PLN, VRT, MAM, PRI...
	comments

Genetic codes file (gencode.dmp):
	genetic code id				-- GenBank genetic code id
	abbreviation				-- genetic code name abbreviation
	name					-- genetic code name
	cde					-- translation table for this genetic code
	starts					-- start codons for this genetic code

"""

sep = re.compile(r'\t\|\t')

RANKS = [
	'no rank',

	'realm',#*
	'subrealm' #*

	'domain', #*
	'superkingdom',
	'kingdom',
	'subkingdom',
	'infrakingdom', #*
	'branch',#*

	'superphylum',
	'phylum',
	'subphylum',
	'infraphylum',#*
	'microphylum',#*

	'superclass',
	'class',
	'subclass',
	'infraclass',
	'subterclass',#*
	'parvclass',#*

	'legion',#*
	'cohort',
	'subcohort',

	'magnorder',#*
	'superorder',
	'order',
	'suborder',
	'infraorder',
	'parvorder',

	'superfamily',
	'family', 
	'subfamily',
	'infrafamily',#*

	'supertribe',#*
	'tribe',
	'subtribe',
	'infratribe',#*

	'genus',
	'subgenus',

	'section',
	'subsection',

	'series',

	'species group',
	'species subgroup',
	'species',
	'subspecies',
	'varietas',
	'variety',#*
	'forma',
	'form',#*
]

class General():
	def get(self):
		return [ getattr(self, attr) for attr in self.__attributes__ ]

	def as_dict(self):
		data = OrderedDict()

		for attr in self.__attributes__:
			data[attr] = getattr(self, attr)
		
		return data

	def _to_pandas(self):
		return self.get()

	def __getstate__(self):
		return [getattr(self.__class__, attr) for attr in self.__shared__] + [getattr(self, attr) for attr in self.__slots__]

	def __setstate__(self, args):
		for pos, attr in enumerate(self.__shared__):
			val = getattr(self.__class__, attr)
			if len(val) == 0:
				# print(self.__class__.__name__.upper() + " :: set state", pos, attr, id(self))
				setattr(self.__class__, attr, args[pos])

		for pos, attr in enumerate(self.__slots__):
			setattr(self, attr, args[pos + len(self.__shared__)])

	def __str__(self):
		return "<{} {}>".format(
			self.__class__.__name__.upper(),
			", ".join(["{}: {}".format(k,v) for k,v in self.as_dict().items()])
		)

	def __repr__(self):
		return str(self)

	def toJSON(self):
		return json.dumps(self.as_dict())

	def __json__(self):
		return self.as_dict()
	
	for_json = __json__  # supported by simplejson


class NameAlts(General):
	__slots__      = ["name_txt", "unique_name", "name_class_id"]
	__attributes__ = ["name_txt", "unique_name", "name_class_id", "name_class"]
	__shared__     = ["classes", "classes_i"]
	
	classes        = OrderedDict()
	classes_i      = OrderedDict()

	def __init__(self, name_txt, unique_name, name_class):
		name_class = name_class.lower()

		if name_class not in self.classes_i:
			self.classes_i[  name_class]    = len(self.classes)
			self.classes[len(self.classes)] = name_class

		name_class_id    = self.classes_i[name_class]

		self.name_txt      = name_txt
		self.unique_name   = unique_name
		self.name_class_id = name_class_id

	@property
	def name_class(self):
		return self.classes[self.name_class_id]


class Name(General):
	__slots__      = ["tax_id", "names"]
	__pds__        = ["tax_id", "name_txt", "unique_name", "name_class"]
	__attributes__ = __slots__
	__shared__     = []

	def __init__(self, tax_id, name_txt, unique_name, name_class):
		self.tax_id = int(tax_id)
		if isinstance(name_txt, list):
			assert isinstance(name_txt   , list)
			assert isinstance(unique_name, list)
			assert isinstance(name_class , list)
			assert len(name_txt) == len(unique_name)
			assert len(name_txt) == len(name_class)
			self.names  = [NameAlts(name_txt[p], unique_name[p], name_class[p]) for p in range(len(name_txt))]
		else:
			self.names  = [NameAlts(name_txt, unique_name, name_class)]

	def merge(self, other):
		self.names.extend(other.names)

	@property
	def txts(self):
		return [n.name_txt    for n in self.names]

	@property
	def unis(self):
		return [n.unique_name for n in self.names]

	@property
	def clss(self):
		return [n.name_class  for n in self.names]	

	def _to_pandas(self):
		txt = self.txts
		uni = self.unis
		cls = self.clss

		return [ self.tax_id, txt, uni, cls ]

	def __str__(self):
		tax_id, names = self.get()
		# names_str     = [ "({} - text '{}' unique name '{}' class ({}) '{}')".format(*([n+1]+list(ns))) for n, ns in enumerate(names) ]
		names_str     = [ "{}: {}".format(n, ns) for n, ns in enumerate(names) ]
		return "<NAME tax id {} ({}) [{}]>".format(tax_id, len(self.names), ", ".join(names_str))


class Division(General):
	__slots__      = ["did" , "cde_id",          "name_id",         "comments"]
	__attributes__ = ["did" , "cde_id", "cde"  , "name_id", "name", "comments"]
	__pds__        = ["did" ,           "cde"  ,            "name", "comments"]
	__shared__     = ["cdes", "cdes_i", "names", "names_i"]

	cdes    = OrderedDict()
	cdes_i  = OrderedDict()
	names   = OrderedDict()
	names_i = OrderedDict()

	def __init__(self, did, cde, name, comments):
		did   = int(did)
		cde  = cde.lower()
		name = name.lower()

		if cde  not in self.cdes_i:
			self.cdes_i[       cde]   = len(self.cdes)
			self.cdes[len(self.cdes)] = cde
		
		if name not in self.names_i:
			self.names_i[       name]   = len(self.names)
			self.names[len(self.names)] = name

		cde_id  = self.cdes_i [cde ]
		name_id = self.names_i[name]

		self.did      = did
		self.cde_id   = cde_id
		self.name_id  = name_id
		self.comments = comments

	@property
	def cde(self):
		return self.cdes[self.cde_id]

	@property
	def name(self):
		return self.names[self.name_id]


class Node(General):
	__slots__      = [
		"tax_id",
		"parent_tax_id",
		"rank_id",
		"embl_code_id",
		"division_id",
		"inherited_div_flag",
		"genetic_code_id",
		"inherited_GC_flag",
		"mitochondrial_genetic_code_id",
		"inherited_MGC_flag",
		"genBank_hidden_flag",
		"hidden_subtree_root_flag"
	]
	__attributes__ = [
		"tax_id",
		"parent_tax_id",
		"rank_id",
		"rank_name",
		"embl_code_id",
		"embl_code",
		"division_id",
		"inherited_div_flag",
		"genetic_code_id",
		"inherited_GC_flag",
		"mitochondrial_genetic_code_id",
		"inherited_MGC_flag",
		"genBank_hidden_flag",
		"hidden_subtree_root_flag"
	]
	__pds__        = [
		"tax_id",
		"parent_tax_id",
		"rank_name",
		"embl_code",
		"division_id",
		"inherited_div_flag",
		"genetic_code_id",
		"inherited_GC_flag",
		"mitochondrial_genetic_code_id",
		"inherited_MGC_flag",
		"genBank_hidden_flag",
		"hidden_subtree_root_flag"
	]
	__shared__     = [
		"embls",
		"embls_i"
	]
	
	embls   = OrderedDict()
	embls_i = OrderedDict()

	def __init__(
			self, 
			tax_id,
			parent_tax_id,
			rank_name,
			embl_code,
			division_id,
			inherited_div_flag,
			genetic_code_id,
			inherited_GC_flag,
			mitochondrial_genetic_code_id,
			inherited_MGC_flag,
			genBank_hidden_flag,
			hidden_subtree_root_flag
		):
		self.tax_id                        = int(tax_id)
		self.parent_tax_id                 = int(parent_tax_id)
		self.division_id                   = int(division_id)
		self.genetic_code_id               = int(genetic_code_id)
		self.mitochondrial_genetic_code_id = int(mitochondrial_genetic_code_id)

		self.inherited_div_flag            = True if inherited_div_flag       in ("1", True) else False
		self.inherited_GC_flag             = True if inherited_GC_flag        in ("1", True) else False
		self.inherited_MGC_flag            = True if inherited_MGC_flag       in ("1", True) else False
		self.genBank_hidden_flag           = True if genBank_hidden_flag      in ("1", True) else False
		self.hidden_subtree_root_flag      = True if hidden_subtree_root_flag in ("1", True) else False

		embl_code                        = embl_code.upper()
		if embl_code not in self.embls_i:
			self.embls_i[embl_code]      = len(self.embls)
			self.embls[len(self.embls)]  = embl_code
		
		embl_code_id      = self.embls_i[embl_code]
		self.embl_code_id = embl_code_id

		rank_name      = rank_name.lower()
		
		assert rank_name in RANKS

		self.rank_id   = RANKS.index(rank_name)

	@property
	def rank_name(self):
		return RANKS[self.rank_id]

	@property
	def embl_code(self):
		return self.embls[self.embl_code_id]


class Row(General):
	__slots__      = [
		"tax_id",
		"parent_tax_id",
		
		"rank_id",
		"rank_name",
		
		"embl_code_id",
		"embl_code",
		
		"division_id",
		"division_cde_id",
		"division_cde",
		"division_name_id",
		"division_name",
		"division_comments",
		"inherited_div_flag",

		"genetic_code_id",
		"genetic_code_abbreviation",
		"genetic_code_name",
		"genetic_code_cde",
		"genetic_code_starts",
		"inherited_GC_flag",

		"mitochondrial_genetic_code_id",
		"mitochondrial_genetic_code_abbreviation",
		"mitochondrial_genetic_code_name",
		"mitochondrial_genetic_code_cde",
		"mitochondrial_genetic_code_starts",
		"inherited_MGC_flag",

		"names",

		"genBank_hidden_flag",
		"hidden_subtree_root_flag",

		"is_parent",
		"asc"
	]
	__pds__        = __slots__
	__attributes__ = __slots__
	__shared__     = []

	def __init__(self, *args):
		(
			self.tax_id,
			self.parent_tax_id,

			self.rank_id,
			self.rank_name,

			self.embl_code_id,
			self.embl_code,
			
			self.division_id,
			self.division_cde_id,
			self.division_cde,
			self.division_name_id,
			self.division_name,
			self.division_comments,
			self.inherited_div_flag,

			self.genetic_code_id,
			self.genetic_code_abbreviation,
			self.genetic_code_name,
			self.genetic_code_cde,
			self.genetic_code_starts,
			self.inherited_GC_flag,

			self.mitochondrial_genetic_code_id,
			self.mitochondrial_genetic_code_abbreviation,
			self.mitochondrial_genetic_code_name,
			self.mitochondrial_genetic_code_cde,
			self.mitochondrial_genetic_code_starts,
			self.inherited_MGC_flag,

			self.names,

			self.genBank_hidden_flag,
			self.hidden_subtree_root_flag,

			self.is_parent,
			self.asc
		) = args


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


def parser_load(fn, cls):
	print("  loading", fn)
	sys.stdout.flush()

	df = pd.read_parquet(fn, engine="fastparquet")

	print("  converting to class")
	sys.stdout.flush()

	# print(pd.concat([df.head(), df.tail()]))

	data = df.to_dict(orient="index", into=OrderedDict)
	data = OrderedDict((k, cls(*[v[s] for s in cls.__pds__])) for k,v in data.items())

	# for e, (k,v) in enumerate(data.items()):
	# 	if e < 5 or e > len(data) - 5:
	# 		print(k, v)

	print("  done")
	sys.stdout.flush()

	return data


def read_dump(filename, has_header=False):
	col_names = []
	with open(filename, "rt") as fhd:
		for linenum, line in enumerate(fhd):
			line = line.strip().strip("\t|")
			if len(line) == 0:
				continue
			cols = sep.split(line)
			if linenum == 0 and has_header:
				col_names = cols
			else:
				if has_header:
					yield col_names, linenum, cols
				else:
					yield linenum, cols


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
		
		parser_save(fn="db_division.pq", cls=Division, data=divisions)
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

	assert asc[node_rank_id] is None and asc[node_rank_id] != tax_id, (tax_id, node_rank_id, RANKS[node_rank_id], level, asc)

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
