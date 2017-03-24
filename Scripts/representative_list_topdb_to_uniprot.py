from __future__ import division

from xml.dom.minidom import parse, parseString
import itertools
import numpy

#TOPDB_uniprot mapper

# Importing the IDs of the representative protein entries from another file.
with open("all TOPDB cdhit90%.txt", "r") as f:
    representative_topDB = f.read().splitlines()
f.closed

# Parse the TOPDB database file.
topdbxml = parse('topdb_all.xml')
proteins_in_xml_list = topdbxml.getElementsByTagName('TOPDB')
uniprotcodes_list = []
done_IDs = []
# Go through each XML entry
for index_of_xml, proteins in enumerate(proteins_in_xml_list):

    # Get topdb ID
    topdb_id_of_this_protein = proteins_in_xml_list[int(index_of_xml)].attributes['ID'].value

    # Only analyse those that are on the representative code list.
    if str(topdb_id_of_this_protein) in representative_topDB:
        uniprot_code = proteins.getElementsByTagName('AC')[0].firstChild.nodeValue
        if topdb_id_of_this_protein not in done_IDs:
            print uniprot_code
            done_IDs.append(topdb_id_of_this_protein)
    else:
        pass
