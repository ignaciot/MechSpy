from owlready2 import *

OUTPUT = "labels/all_labels.tsv"

onto_path.append("labels")

with open(OUTPUT, 'w') as new_fd:
    onto = get_ontology("labels/ro-core.owl")
    onto.load()
    for prop in onto.properties():
        if len(prop.label) > 0:
            new_fd.write("<%s>\t%s\n" % (prop.iri, prop.label[0]))

    onto = get_ontology("labels/go.owl")
    onto.load()
    for cl in onto.classes():
        if len(cl.label) > 0:
            new_fd.write("<%s>\t%s\n" % (cl.iri, cl.label[0]))

    with open("labels/Homo_sapiens.gene_info", 'r') as fd:
        for line in fd.readlines():
            chunks = line.split('\t')
            new_fd.write("<http://purl.uniprot.org/geneid/%s>\t%s (%s)\n" % (chunks[1], chunks[2], chunks[8]))


    new_fd.write("<http://purl.obolibrary.org/obo/BFO_0000051>\thas part in\n")
    new_fd.write("<http://purl.obolibrary.org/obo/BFO_0000056>\tparticipates in\n")
    new_fd.write("<http://purl.obolibrary.org/obo/BFO_0002180>\thas component\n")
    new_fd.write("<http://purl.obolibrary.org/obo/RO_0002434>\tinteracts with\n")
    new_fd.write("<http://www.w3.org/2000/01/rdf-schema#subClassOf>\tis a\n")
    
    with open("reactome/UniProt2Reactome_All_Levels.txt", 'r') as fd:
        for line in fd.readlines():
            chunks = line.split('\t')
            new_fd.write("<https://reactome.org/content/detail/%s>\t%s\n" % (chunks[1], chunks[3]))

