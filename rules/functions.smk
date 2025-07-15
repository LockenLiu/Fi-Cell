from snakemake.utils import min_version

import sys
import os
import platform
import re
import csv
import gzip
import warnings

def build_dict_from_key_value_pairs_sc(list_of_dicts):
	'''
	Each dict in the list MUST contain the keys 'id', 'sce_path', and 'cts_score_path'.
	path will be converted to absolute paths.
	Takes the list of dictionaries and makes it into a new dictionary where the keys are the id values from each dictionary and the values are each dictionary
	e.g. [{"id":"a", "value": 1}, {"id":"b","value":2}] ->
	{"a":{"id":"a", "value": 1}, "b":{"id":"b","value":2}}
	'''
	out_dict = {}
	for d in list_of_dicts:
		d['sce_path'] = os.path.abspath(d['sce_path'])
		d['cts_score_path'] = os.path.abspath(d['cts_score_path'])
		out_dict[d['id']] = d
	return(out_dict)

def build_dict_from_key_value_pairs_gwas(list_of_dicts):
	'''
	Each dict in the list MUST contain the keys 'id', and 'sce_path'.
	path will be converted to absolute paths.
	Takes the list of dictionaries and makes it into a new dictionary where the keys are the id values from each dictionary and the values are each dictionary
	e.g. [{"id":"a", "value": 1}, {"id":"b","value":2}] ->
	{"a":{"id":"a", "value": 1}, "b":{"id":"b","value":2}}
	'''
	out_dict = {}
	for d in list_of_dicts:
		d['path'] = os.path.abspath(d['path'])
		out_dict[d['id']] = d
	return(out_dict)

##############################
GWAS_SUMSTATS = build_dict_from_key_value_pairs_gwas(config['GWAS_SUMSTATS_INPUT'])
SC_SUMSTATS = build_dict_from_key_value_pairs_sc(config['SCE_INPUT'])
