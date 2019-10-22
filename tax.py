#!/usr/bin/env python3
import sys
import re

from collections import OrderedDict

import simplejson as json

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
	division cde			-- GenBank division code (three characters)
	division name			-- e.g. BCT, PLN, VRT, MAM, PRI...
	comments

gencode.dmp
-----------
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

	def _to_pandas(self):
		return [ getattr(self, attr) for attr in self.__pds__ ]

	def as_dict(self):
		data = OrderedDict()

		for attr in self.__attributes__:
			data[attr] = getattr(self, attr)
		
		return data

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

