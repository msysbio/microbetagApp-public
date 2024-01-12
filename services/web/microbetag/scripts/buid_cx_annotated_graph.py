"""
Aim:
    Gets an edgelist as input along with a df mentioning the NCBI Tax ids and the GTDB ids of the corresponding taxa
    and based on the annotation types asked, it build a cx-format graph that is going to be the final return object of microbetag

Author: 
    Haris Zafeiropoulos
"""
try:
    from .variables import *
    from .utils import *
except ImportError:
    from variables import *
    from utils import *
import pandas as pd
import numpy as np
import ast
import json
import numpy as np
import os


def build_cx_annotated_graph(edgelist_as_df, edgelist_as_a_list_of_dicts, seq_map, cfg, out_dir, children_df=None):
    """
    Builds the annotated network object in the .cx format
    (https://cytoscape.org/cx/cx2/specification/cytoscape-exchange-format-specification-(version-2)/)
    """
    base_cx = []
    init = {}
    init["numberVerification"] = [{"longNumber": 281474976710655}]
    base_cx.append(init)

    metadata = {}
    metadata["metaData"] = [
        {"name": "cyTableColumn", "version": "1.0"},
        {"name": "nodes", "version": "1.0"},
        {"name": "edges", "version": "1.0"},
        {"name": "nodeAttributes", "version": "1.0"},
        {"name": "edgeAttributes", "version": "1.0"},
        {"name": "networkAttributes", "version": "1.0"},
        {"name": "cartesianLayout", "version": "1.0"},
    ]
    base_cx.append(metadata)

    # GET ALL COLUMNS OF ALL TABLES
    table_columns = {}
    table_columns["cyTableColumn"] = []
    cyTableColumns = [
        # for node table mandatory
        {"applies_to": "node_table", "n": "@id", "d": "string"},
        {"applies_to": "node_table", "n": "name", "d": "string"},
        {"applies_to": "node_table", "n": "shared name", "d": "string"},
        {"applies_to": "node_table", "n": "display name", "d": "string"},
        {"applies_to": "node_table", "n": "microbetag::taxon name", "d": "string"},
        {"applies_to": "node_table", "n": "microbetag::namespace"},
        {"applies_to": "node_table", "n": "microbetag::taxonomy", "d": "string"},
        {"applies_to": "node_table", "n": "microbetag::ncbi-tax-id", "d": "string"},
        {"applies_to": "node_table", "n": "microbetag::ncbi-tax-level", "d": "string"},
        {"applies_to": "node_table", "n": "microbetag::gtdb-genomes", "d": "list_of_string"},
        # for edge table mandatory
        {"applies_to": "edge_table", "n": "shared name"},
        {"applies_to": "edge_table", "n": "name"},
        {"applies_to": "edge_table", "n": "interaction type"},
        {"applies_to": "edge_table", "n": "weight::weight", "d": "double"}  # flashweave score
    ]

    """CHILDREN GENOMES AND NCBI IDS"""
    if cfg["get_children"] and children_df is not None:
        cyTableColumns.append({"applies_to": "node_table", "n": "microbetag::children-ncbi-ids", "d": "list_of_string"})
        cyTableColumns.append({"applies_to": "node_table", "n": "microbetag::children-gtdb-genomes", "d": "list_of_string"})

    """Phen traits"""
    if cfg["phenDB"]:
        df_traits_table = pd.read_csv(os.path.join(out_dir, "phen_predictions/phen_traits.tsv"), sep="\t")
        phen_traits = [c for c in df_traits_table.columns.to_list() if "Score" not in c]
        phen_traits.remove("NCBI_ID")
        phen_traits.remove("gtdb_id")
        for term in phen_traits:
            cyTableColumns.append({"applies_to": "node_table", "n": "::".join(["phendb", term]), "d": "boolean"})
            cyTableColumns.append({"applies_to": "node_table", "n": "::".join(["phendbScore", "".join([term, "Score"])]), "d": "double"})
            cyTableColumns.append({"applies_to": "node_table", "n": "::".join(["phendbScore", "".join([term, "Std"])]), "d": "double"})
        df_traits_table["NCBI_ID"] = pd.to_numeric(df_traits_table["NCBI_ID"], errors='coerce')

    """FAPROTAX traits"""
    if cfg["faprotax"]:
        assignments_per_seqId = seqId_faprotax_functions_assignment(os.path.join(out_dir, "faprotax/sub_tables/"))
        faprotax_traits = list(flatten_list(assignments_per_seqId.values()))
        for term in faprotax_traits:
            column = {"applies_to": "node_table", "n": "::".join(["faprotax", term]), "d": "boolean"}
            cyTableColumns.append(column)

    """COMPLEMENTS"""
    if cfg["pathway_complement"]:
        complements_dict = json.load(open(os.path.join(out_dir, "path_compl/complements.json")))

        descrps = pd.read_csv(os.path.join(MAPPINGS, "module_descriptions"), sep="\t", header=None)
        descrps.columns = ["category", "moduleId", "description"]
        column_order = ["moduleId", "description", "category"]
        descrps = descrps[column_order]
        for n_pair, g_pair in complements_dict.items():
            for genomes, compls in g_pair.items():
                pairs_cats = set()  # NOTE use it if needed as last entry for each genome pair so cytoscapeApp can build filtering
                for index, compl in enumerate(compls):
                    triplet = descrps[descrps["moduleId"] == compl[0][0]]
                    try:
                        triplet = triplet.values.tolist()[0]
                        pairs_cats.add(triplet[2])
                        extend = triplet + compl[0][1:]
                        complements_dict[n_pair][genomes][index] = [extend]
                    except ValueError:
                        print(">>>>>", compl)
                        logging.warn("Module with id:", compl[0][0], " not in current the module_descriptions file.")
                        extend = [compl[0][0]] + ["N/A", "N/A",] + compl[0][1:]
                        complements_dict[n_pair][genomes][index] = [extend]

        # ast.literal_eval function to evaluate the key as a Python literal expression and converts it into an object.
        # to safely evaluate a string containing a Python literal or container display (e.g., a dictionary or list) without running arbitrary code.
        # It's a safer alternative to eval when dealing with untrusted input.
        complements_keys = [ast.literal_eval(x) for x in list(complements_dict.keys())]
        complements_dict = dict(zip(complements_keys, complements_dict.values()))

        for edge in edgelist_as_a_list_of_dicts:
            if (
                (seq_map[seq_map["seqId"] == edge["node_a"]]["gtdb_gen_repr"].notna()).any() and
                (seq_map[seq_map["seqId"] == edge["node_b"]]["gtdb_gen_repr"].notna()).any()
            ):
                genomes_a = seq_map[seq_map["seqId"] == edge["node_a"]]["gtdb_gen_repr"].to_list()[0]
                genomes_b = seq_map[seq_map["seqId"] == edge["node_b"]]["gtdb_gen_repr"].to_list()[0]
                for genome_a in genomes_a:
                    for genome_b in genomes_b:
                        edge_col = {"applies_to": "edge_table", "n": "".join(["compl::", genome_a, ":", genome_b]), "d": "list_of_string"}
                        cyTableColumns.append(edge_col)
                        # example of edge_col:  {"applies_to": "edge_table", "n": "compl::GCF_000174815.1:GCF_003253005.1", "d": "list_of_string"}

    """SEED SCORES"""
    if cfg["seed_scores"]:
        seed_scores_dict = json.load(open(os.path.join(out_dir, "seed_scores/seed_scores.json")))
        seed_scores_keys = [ast.literal_eval(x) for x in list(seed_scores_dict.keys())]
        seed_scores_dict = dict(zip(seed_scores_keys, seed_scores_dict.values()))
        cyTableColumns.append({"applies_to": "edge_table", "n": "::".join(["seed", "competition"]), "d": "double"})
        cyTableColumns.append({"applies_to": "edge_table", "n": "::".join(["seed", "competition-std"]), "d": "double"})
        cyTableColumns.append({"applies_to": "edge_table", "n": "::".join(["seed", "cooperation"]), "d": "double"})
        cyTableColumns.append({"applies_to": "edge_table", "n": "::".join(["seed", "cooperation-std"]), "d": "double"})
        for seedFindings in seed_scores_dict.values():
            if len(seedFindings) > 0:
                for scenario in seedFindings.values():
                    cyTableColumns.append(
                        {"applies_to": "edge_table", 
                            "n": "".join(["seedCompl::", scenario["genome-A"], ":", scenario["genome-B"]]), 
                            "d": "list_of_string"
                        }
                    )

    """MANTA CLUSTERS"""
    if cfg["manta"]:
        cartesianLayout = {}; cartesianLayout["cartesianLayout"] = []
        m1 = {"applies_to": "node_table", "n": "::".join(["manta", "cluster"]), "d": "double"}
        m2 = {"applies_to": "node_table", "n": "::".join(["manta", "assignment"]), "d": "string"}
        cyTableColumns.append(m1)
        cyTableColumns.append(m2)
        manta_output_file = "/".join([out_dir, 'manta_annotated.cyjs'])
        manta_net = read_cyjson(manta_output_file)
        clusters = list(manta_net.nodes(data="cluster"))
        assignments = list(manta_net.nodes(data="assignment"))
        positions = list(manta_net.nodes(data="position"))
        manta_annotations = {}
        for pair in clusters:
            manta_annotations[pair[0]] = {}
            manta_annotations[pair[0]]["cluster"] = pair[1]
        for pair in assignments:
            manta_annotations[pair[0]]["assignment"] = pair[1]
        for pair in positions:
            manta_annotations[pair[0]]["position"] = pair[1]

    table_columns["cyTableColumn"] = cyTableColumns
    base_cx.append(table_columns)

    # NETWORK TABLE
    networkAttributes = {}
    networkAttributes["networkAttributes"] = []
    networkAttributes["networkAttributes"].append({"n": "database", "v": "microbetagDB"})
    networkAttributes["networkAttributes"].append({"n": "network type", "v": "annotated network of microbial co-occurrences"})
    networkAttributes["networkAttributes"].append({"n": "name", "v": "microbetag network"})
    networkAttributes["networkAttributes"].append({"n": "uri", "v": "https://hariszaf.github.io/microbetag/"})
    networkAttributes["networkAttributes"].append({"n": "version", "v": "1.0"})
    base_cx.append(networkAttributes)

    # NODES TABLE
    nodes = {}; nodes["nodes"] = []
    nodeAttributes = {}; nodeAttributes["nodeAttributes"] = []

    set_of_nodes = set(edgelist_as_df["node_a"].to_list() + edgelist_as_df["node_b"].to_list())
    node_counter = 1000
    seq_to_nodeID = {}

    for seq in set_of_nodes:
        node_counter += 1
        case = seq_map[seq_map["seqId"] == seq]
        node = {"@id": node_counter, "n": seq}
        seq_to_nodeID[seq] = node_counter
        nodes["nodes"].append(node)

        # NCBI attribute
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "@id", "v": seq, "d": "string"})
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "shared name", "v": seq, "d": "string"})
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "display name", "v": seq, "d": "string"})
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "microbetag::taxon name", "v": case["extendedSpecies"].item(), "d": "string"})
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "microbetag::taxonomy", "v": case["taxonomy"].item(), "d": "string"})

        # Clusters
        if cfg["manta"]:
            if seq in manta_annotations.keys():
                nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "manta::cluster", "v": str(manta_annotations[seq]["cluster"]), "d": "double"})
                nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "manta::assignment", "v": manta_annotations[seq]["assignment"], "d": "string"})
                cartesianLayout["cartesianLayout"].append({"node": node_counter, "x": manta_annotations[seq]["position"]["x"], "y": manta_annotations[seq]["position"]["y"]})

        if isinstance(case["ncbi_tax_id"].item(), pd._libs.missing.NAType):
            continue

        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "microbetag::ncbi-tax-id", "v": str(case["ncbi_tax_id"].item()), "d": "string"})
        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "microbetag::ncbi-tax-level", "v": case["ncbi_tax_level"].item(), "d": "string"})

        if not isinstance(case["gtdb_gen_repr"].item(), type(None)):
                nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "microbetag::gtdb-genomes", "v": case["gtdb_gen_repr"].item(), "d": "list_of_string"})

        if cfg["get_children"] and children_df is not None:
            ch_case = children_df[children_df["parent_ncbi_tax_id"] == case["ncbi_tax_id"].item()].dropna()
            if len(ch_case["child_ncbi_tax_id"].to_list()) > 0:
                children_ncbi_ids = [str(x) for x in list(flatten_list(ch_case["child_ncbi_tax_id"].to_list(), []))]
                nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "microbetag::children-ncbi-ids", "v": children_ncbi_ids, "d": "list_of_string"})
            if len(ch_case["gtdb_gen_repr"].to_list()):
                children_genomes = list(flatten_list(ch_case["gtdb_gen_repr"].to_list(), []))
                nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "microbetag::children-gtdb-genomes", "v": children_genomes, "d": "list_of_string"})

        if cfg["phenDB"]:
            df_traits_table["NCBI_ID_str"] = df_traits_table["NCBI_ID"].astype(str)
            traits_case = df_traits_table[df_traits_table["NCBI_ID_str"] == case["ncbi_tax_id"].item()]

            if traits_case.empty:
                continue
            for trait in phen_traits:
                trait_value = list(set(traits_case[trait].to_list()))
                trait_score = traits_case["".join([trait, "Score"])].mean()
                trait_std = traits_case["".join([trait, "Score"])].std()
                if len(trait_value) == 1:
                    trait_value = True if trait_value[0] == "YES" else False
                    nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "::".join(["phendb", trait]), "v": trait_value, "d": "boolean"})
                    nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "::".join(["phendbScore", "".join([trait, "Score"])]), "v": str(trait_score), "d": "double"})
                    if traits_case.shape[0] == 1:
                        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "::".join(["phendbScore", "".join([trait, "Std"])]), "v": "0.0", "d": "double"})
                    else:
                        nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "::".join(["phendbScore", "".join([trait, "Std"])]), "v": str(trait_std), "d": "double"})

        if cfg["faprotax"]:
            if seq in assignments_per_seqId:
                for term in assignments_per_seqId[seq]:
                    nodeAttributes["nodeAttributes"].append({"po": node_counter, "n": "::".join(["faprotax", term]), "v": True, "d": "boolean"})

    base_cx.append(nodes); base_cx.append(nodeAttributes)
    if cfg["manta"]:
        base_cx.append(cartesianLayout)

    # ADD EDGES AND THEIR ATTRIBUTES
    # [NOTE] Each seqId corresponds to a single node id, however each co-occurrence association from Flashweave or any other tool leads up to 3 edges in
    # the microbetag network; that is because in the case of an mspecies-mspecies association, we build 2 extra edges to describe this pairwised (A-B) association:
    # one where we have A as source and B as target and another that is vice-versa; we do that, so it is clear which complements have A as beneficiary species and B as donor and
    # the other way around. The same applies for the seed scores.
    # [NOTE] With respect to the path complementarities, we need to keep in mind that species A may benefit from species B but not the other way around.
    # Thus, we may have only one of the two extra edges.
    edges = {}; edges["edges"] = []
    edgeAttributes = {}; edgeAttributes["edgeAttributes"] = []
    edge_counter = node_counter + 1000
    for case in edgelist_as_a_list_of_dicts:
        id_a = seq_to_nodeID[case["node_a"]]
        id_b = seq_to_nodeID[case["node_b"]]
        edge = {"@id": edge_counter, "s": id_a, "t": id_b, "i": "cooccurrence/depletion"}
        edges["edges"].append(edge)

        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "weight::weight", "v": str(case["score"]), "d": "double"})
        if float(case["score"]) > 0:
            edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "shared name", "v": " ".join([case["node_a"], "(cooccurss with)", case["node_b"]]), "d": "string"})
            edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "interaction type", "v": "cooccurrence", "d": "string"})
        else:
            edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "shared name", "v": " ".join([case["node_a"], "(depletes)", case["node_b"]]), "d": "string"})
            edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "interaction type", "v": "depletion", "d": "string"})

        ncbi_a = seq_map[seq_map["seqId"] == case["node_a"]]["ncbi_tax_id"].item()
        ncbi_b = seq_map[seq_map["seqId"] == case["node_b"]]["ncbi_tax_id"].item()
        if isinstance(ncbi_a, pd._libs.missing.NAType) or isinstance(ncbi_b, pd._libs.missing.NAType):
            continue

        if ncbi_a == "<NA>" or ncbi_b == "<NA>":
            continue
        ncbi_pair_as_tuple_a_b = (int(ncbi_a), int(ncbi_b))
        ncbi_pair_as_tuple_b_a = (int(ncbi_b), int(ncbi_a))
        edge_counter += 1

        if cfg["pathway_complement"] or cfg["seed_scores"]:

            check = False

            """Edge for A -> B"""
            pot_edge = {"@id": (edge_counter), "s": id_a, "t": id_b, "i": "comp_coop"}

            # Path complements A -> B
            if cfg["pathway_complement"]:
                if ncbi_pair_as_tuple_a_b in complements_keys:
                    edges["edges"].append(pot_edge)
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "shared name", "v": " ".join([case["node_a"], "(completes/competes with)", case["node_b"]]), "d": "string"})
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "interaction type", "v": "completes/competes with", "d": "string"})
                    check = True
                    for genome_combo, complements in complements_dict[ncbi_pair_as_tuple_a_b].items():
                        genome_combo = ast.literal_eval(genome_combo)
                        attr = "".join(["compl::", genome_combo[0], ":", genome_combo[1]])
                        merged_compl = ["^".join(gcompl[0]) for gcompl in complements]
                        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": attr, "v": merged_compl, "d": "list_of_string"})

            # Seed scores A -> B
            if cfg["seed_scores"]:
                if ncbi_pair_as_tuple_a_b in seed_scores_keys:
                    if pot_edge not in edges["edges"]:
                        edges["edges"].append(pot_edge)
                        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "shared name", "v": " ".join([case["node_a"], "(completes/competes with)", case["node_b"]]), "d": "string"})
                        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "interaction type", "v": "completes/competes with", "d": "string"})
                        check = True

                    comp = []; coop = []
                    for index in seed_scores_dict[ncbi_pair_as_tuple_a_b].values():
                        if "competition" in index.keys():
                            comp.append(index["competition"])
                            coop.append(index["cooperation"])
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "::".join(["seed", "competition"]), "v": str(np.mean(comp)), "d": "double"})
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "::".join(["seed", "competition-std"]), "v": str(np.std(comp)), "d": "double"})
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "::".join(["seed", "cooperation"]), "v": str(np.mean(coop)), "d": "double"})
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "::".join(["seed", "cooperation-std"]), "v": str(np.std(coop)), "d": "double"})

                    for scenario in seed_scores_dict[ncbi_pair_as_tuple_a_b].values():
                        if "complements" in scenario.keys():
                            attr = "".join(["seedCompl::", scenario["genome-A"], ":", scenario["genome-B"]])
                            merged_compl = ["^".join(gcompl) for gcompl in scenario["complements"]]
                            edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": attr, "v": merged_compl, "d": "list_of_string"})

            if check:
                edge_counter += 1
                check = False

            """Edge for B -> A"""
            pot_edge = {"@id": (edge_counter), "s": id_b, "t": id_a, "i": "comp_coop"}

            # Path complements B -> A
            if cfg["pathway_complement"]:
                if ncbi_pair_as_tuple_b_a in complements_keys:
                    edges["edges"].append(pot_edge)
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "shared name", "v": " ".join([case["node_b"], "(completes/competes with)", case["node_a"]]), "d": "string"})
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "interaction type", "v": "completes/competes with", "d": "string"})
                    check = True
                    for genome_combo, complements in complements_dict[ncbi_pair_as_tuple_b_a].items():
                        genome_combo = ast.literal_eval(genome_combo)
                        merged_compl = ["^".join(gcompl[0]) for gcompl in complements]
                        attr = "".join(["compl::", genome_combo[0], ":", genome_combo[1]])
                        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": attr, "v": merged_compl, "d": "list_of_string"})

            # Seed scores B -> A
            if cfg["seed_scores"]:
                if ncbi_pair_as_tuple_b_a in seed_scores_keys:
                    if pot_edge not in edges["edges"]:
                        edges["edges"].append(pot_edge)
                        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "shared name", "v": " ".join([case["node_b"], "(completes/competes with)", case["node_a"]]), "d": "string"})
                        edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "interaction type", "v": "completes/competes with", "d": "string"})
                        check = True
                    comp = []; coop = []
                    for index in seed_scores_dict[ncbi_pair_as_tuple_b_a].values():
                        if "competition" in index.keys():
                            comp.append(index["competition"])
                            coop.append(index["cooperation"])
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "::".join(["seed", "competition"]), "v": str(np.mean(comp)), "d": "double"})
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "::".join(["seed", "competition-std"]), "v": str(np.std(comp)), "d": "double"})
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "::".join(["seed", "cooperation"]), "v": str(np.mean(coop)), "d": "double"})
                    edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": "::".join(["seed", "cooperation-std"]), "v": str(np.std(coop)), "d": "double"})

                    for scenario in seed_scores_dict[ncbi_pair_as_tuple_b_a].values():
                        if "complements" in scenario.keys():
                            attr = "".join(["seedCompl::", scenario["genome-B"], ":", scenario["genome-A"]])
                            # NOTE: in case a map has no description or for any reason we have a nan value from the dfs this will fail
                            merged_compl = ["^".join(gcompl) for gcompl in scenario["complements"]]
                            edgeAttributes["edgeAttributes"].append({"po": edge_counter, "n": attr, "v": merged_compl, "d": "list_of_string"})

            if check:
                edge_counter += 1

    base_cx.append(edges); base_cx.append(edgeAttributes)

    # POST-metadata
    post_metadata = {}
    post_metadata["metaData"] = []
    post_metadata["metaData"].append({"name": "nodeAttributes", "elementCount": len(nodeAttributes["nodeAttributes"]), "version": 1.0})
    post_metadata["metaData"].append({"name": "edgeAttributes", "elementCount": len(edgeAttributes["edgeAttributes"]), "version": 1.0})
    post_metadata["metaData"].append({"name": "cyTableColumn", "elementCount": len(table_columns["cyTableColumn"]), "version": 1.0})
    post_metadata["metaData"].append({"name": "edges", "elementCount": len(edges["edges"]), "idCounter": node_counter + 1000, "version": 1.0})
    post_metadata["metaData"].append({"name": "nodes", "elementCount": len(nodes["nodes"]), "idCounter": 1001, "version": 1.0})
    post_metadata["metaData"].append({"name": "networkPropernetworkAttributesties", "elementCount": len(networkAttributes["networkAttributes"]), "version": "1.0"})
    if cfg["manta"]:
        post_metadata["metaData"].append({"name": "cartesianLayout", "elementCount": len(cartesianLayout["cartesianLayout"]), "version": 1.0})

    base_cx.append(post_metadata)

    # Status
    status = {}; status["status"] = []
    status["status"].append({"error": "", "success": True})
    base_cx.append(status)

    return base_cx


def seqId_faprotax_functions_assignment(path_to_subtables):
    """
    Parse the sub tables of the faprotax analysis
    to assign the biological processes related to each sequence id
    """
    seqId_faprotax_assignments = {}
    for subtable_file in os.listdir(path_to_subtables):
        f = os.path.join(path_to_subtables, subtable_file)
        process_name = subtable_file.split(".")[0].replace("_", " ")
        table_file = open(f, "r")
        table_file = table_file.readlines()
        for line in table_file[2:]:
            seqId = line.split("\t")[1]
            if seqId not in seqId_faprotax_assignments:
                seqId_faprotax_assignments[seqId] = [process_name]
            else:
                seqId_faprotax_assignments[seqId].append(process_name)
    return seqId_faprotax_assignments


if __name__ == "__main__":
    """
    Build an annotated .cx using already build objects
    """
    import sys
   

    input_dir = os.path.join(sys.argv[1], "tests/build_annotated_graph")
    edgelist_as_df = pd.read_csv( os.path.join(input_dir, "network.edgelist"), sep="\t", index_col = 0 )
    edgelist_as_a_list_of_dicts = edge_list_of_ncbi_ids(os.path.join(input_dir, "network.edgelist"))
    seq_map = pd.read_csv( os.path.join(input_dir, "seq_map.csv"), sep="\t", index_col = 0 )
    cfg = {
        "taxonomy":"microbetag_prep", "delimiter": ";", "faprotax": False, "phenDB": True, "pathway_complement": True, 
        "seed_scores": False, "manta": False, "get_children": False, "sensitive": True, "heterogeneous": False
    }
    out_dir = input_dir
    annotated_network = build_cx_annotated_graph(
        edgelist_as_df,
        edgelist_as_a_list_of_dicts,
        seq_map,
        cfg,
        out_dir
    )

