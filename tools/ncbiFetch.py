#!/usr/bin/env python
## NCBI FETCH SCRIPT
## Utility to fetch data from NCBI and return in a tabular format.
## Created by Peter Lais on 02/20/24.

# All gene symbols are CASE SENSITIVE.

import argparse
import requests
import urllib
import sys
import pandas as pd
from bs4 import BeautifulSoup
from io import StringIO
from rich.progress import track, Progress

import time

# ISSUE: with many genes, paging limits the genes you can consider. Auto-pull all pages.
def _get_all_response(request_string, curr_json):
    responses = {**curr_json}
    tmp_ = responses
    while "next_page_token" in tmp_.keys():
        tmp_ = requests.get(request_string + "&page_token=%s" % tmp_["next_page_token"]).json()
        responses["reports"].extend(tmp_["reports"])
    responses.pop("next_page_token", None)
    return responses

# Given a symbol string and response JSON, return the subset corresponding to the symbol of interest.
# Decided to add case-insensitivity since we manually specify the taxon anyway.
def _find_symbol_dict(symbol, response_json):

    # Get the main symbol dictionary from the API. If it fails, it fails.
    found = False
    for entry in response_json.get("reports", []):
        gene_info = entry["gene"]
        if gene_info["symbol"].lower() == symbol.lower() or symbol.lower() in [x.lower() for x in gene_info.get("synonyms", [])]:
            found = True
            break
    if found == False: return {} 
    
    # Add the symbol description from the NCBI website.
    gene_info["summary"] = _id_to_description(gene_info.get("gene_id", None))

    return gene_info

# Given a symbol string and response JSON, return the subset corresponding to the symbol of interest.
def _find_id_dict(response_json):

    # Get the complete symbol information for each gene. We should have a gene_id for each id we pass in.
    gene_info = {}
    for entry in track(response_json.get("reports", []), description="Processing gene IDs..."):
        tmp = entry["gene"]
        tmp["summary"] = _id_to_description(tmp["gene_id"])
        gene_info[tmp["gene_id"]] = tmp

    return gene_info

def _id_to_description(gene_id):

    # If nothing comes in, nothing goes out.
    if gene_id is None: return None

    # Second request: get the data.
    soup = BeautifulSoup(requests.get(f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}").content, "html.parser")
    tree = soup.find("dl", id = "summaryDl")

    # Once the child is found, then print the next element that has a non-whitespace only string.
    # If otherwise, return None.
    found_target = False
    for child in tree.children:
        if found_target and child.string is not None and child.string.strip():
            return child.string
        if child.string == "Summary":
            found_target = True
    return None

def symbols_to_info(gene_symbols, taxon):
    return _symbols_to_info_helper(gene_symbols, taxon, True)

# Create a gene information dictionary from a list of symbols.
def _symbols_to_info_helper(gene_symbols, taxon, progress):

    # First try to get the symbols using one request. Then, try with multiple.

    # Get the response json for the complete list of gene symbols.
    with Progress(disable=not progress) as progress:
        id = progress.add_task("Processing gene symbols...", total=len(gene_symbols))
        # Begin the try/catch of trying one/multiple requests.
        try:
            # Format the arguments and fetch the JSON.
            if not hasattr(gene_symbols, "__len__"): gene_symbols = [gene_symbols]
            gene_string = urllib.parse.quote(",".join(gene_symbols), safe="")
            taxon_string = urllib.parse.quote(taxon, safe="")
            req = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/gene/symbol/{gene_string}/taxon/{taxon_string}?table_fields=gene-id"
            response_json = requests.get(req).json()
            response_json = _get_all_response(req, response_json)

            # Make a dictionary containing information for each gene.
            gene_info = {}
            for symbol in progress.track(gene_symbols, task_id=id):
                gene_info[symbol] = _find_symbol_dict(symbol, response_json)
            return gene_info
        
        except requests.exceptions.JSONDecodeError as err:
            progress.reset(id, description="Processing gene symbols using separate requests...")
            if len(gene_symbols) == 1:
                # The one-by-one strategy has been tried, so just return nothing.
                return {gene_symbols[0]: {}}
            else:
                # Try scraping information for the genes one-by-one.
                # Each call will individually fail if it is not recognized.
                gene_info = {}
                for symbol in progress.track(gene_symbols, task_id=id):
                    gene_info[symbol] = _symbols_to_info_helper([symbol], taxon, False)[symbol]
                return gene_info


def ids_to_info(ids):
    # Get the response json for the complete list of gene symbols.
    id_string = urllib.parse.quote(",".join(ids), safe="")
    response_json = requests.get(f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/gene/id/{id_string}").json()

    # Make a dictionary containing information for each gene. Add the missing info as None.
    tmp = _find_id_dict(response_json)
    for leftover_id in set(ids).difference(tmp.keys()):
        tmp[leftover_id] = {}
    return tmp


def create_summary_df(gene_info_dict):
    return pd.DataFrame({"Gene Symbol": [entry.get("symbol") for entry in gene_info_dict.values()],
                         "NCBI ID": [entry.get("gene_id") for entry in gene_info_dict.values()],
                         "Description": [entry.get("description") for entry in gene_info_dict.values()],
                         "Taxonomy Name": [entry.get("taxname") for entry in gene_info_dict.values()],
                         "Gene Type": [entry.get("type") for entry in gene_info_dict.values()],
                         "Orientation": [entry.get("orientation") for entry in gene_info_dict.values()],
                         "NCBI Gene Summary": [entry.get("summary") for entry in gene_info_dict.values()]},
                         index = gene_info_dict.keys())

def process_gene_list(fp):
    with open(fp) as fd:
        genes = fd.read()
    return [gene.strip() for gene in genes.split(",")]

if __name__ == "__main__":

    # Create the argument parser.
    parser = argparse.ArgumentParser(prog = "NCBI Fetch",
                                     description = "Fetch gene data from NCBI using aliases.",
                                     epilog = "Created by Peter Lais on 02/20/2024.")
    parser.add_argument("-t", "--taxon", required=False, type = str)
    parser.add_argument("-s", "--use-symbols", action="store_true", required=False)
    #parser.add_argument("-a", "--use-accessions", action="store_true", required=False)
    parser.add_argument("-i", "--use-ids", action="store_true", required=False)
    parser.add_argument("-f", "--infile", required=False, type = str)
    parser.add_argument("-o", "--outfile", required=False, type = str)
    parser.add_argument("-v", "--csv", action="store_true", required=False)
    parser.add_argument("-x", "--excel", action="store_true", required=False)
    parser.add_argument("-S", "--sep", required=False, type=str, default="\n")
    parser.add_argument("-c", "--column", default=0, type=int, required=False)

    # Parse and check the arguments.
    args = parser.parse_args()
    if args.use_symbols and args.taxon is None:
        raise argparse.ArgumentError(None, message = "If using symbols, must supply a taxon (NCBI taxonomy identifier, common name, or scientific name).\n"
                                                     "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/rest-api/#get-/gene/taxon/-taxon-")
    elif (not args.use_symbols and not args.use_ids) or (args.use_symbols and args.use_ids):
        raise argparse.ArgumentError(None, message = "Please specify either the --use-symbols or --use-ids flag, not both or neither.")
    
    # Get the information in.
    if args.infile is not None:
        with open(args.infile) as fd:
            content = fd.read()
    else:
        content = sys.stdin.read()
    
    # Split the information appropriately.
    # Skip blank entries resulting from errant commas and such.
    if not (args.csv or args.excel):
        elem_list = sorted(set([elem.strip() for elem in content.split(args.sep) if elem.strip()]))
    else:
        buf = StringIO(content)
        tab = pd.read_csv(buf) if args.csv else pd.read_excel(buf)
        elem_list = sorted(set([elem.strip() for elem in tab.iloc[:,args.column] if elem.strip()]))

    # Create the dataframe.
    if args.use_symbols:
        out = create_summary_df(symbols_to_info(elem_list, args.taxon))
    else:
        out = create_summary_df(ids_to_info(elem_list))

    # Get the information out.
    if args.outfile is not None:
        out.to_csv(args.outfile)
    else:
        print(out)
