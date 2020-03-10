
# Mechanism 1: Triggering of caspase-mediated apoptosis via release of cytochrome C
# -----------
M1_NODE_LABELS = [
            'positive regulation of mitochondrial membrane permeability',
            'positive regulation of release of cytochrome c from mitochondria',
            'caspase activation',
            'apoptotic DNA fragmentation',
            'intrinsic apoptotic signaling pathway in response to DNA damage',
        ]
M1_NODES = [
            "<http://purl.obolibrary.org/obo/GO_0035794>",
            "<http://purl.obolibrary.org/obo/GO_0090200>",
            "<http://purl.obolibrary.org/obo/GO_0006919>",
            "<http://purl.obolibrary.org/obo/GO_0006309>",
            "<http://purl.obolibrary.org/obo/GO_0008630>",
            ]

# Mechanism 2: ATP depletion from calcium homeostasis disruption, resulting in necrosis
# -----------
M2_NODE_LABELS = [
            'regulation of calcium ion transport',
            'positive regulation of cytosolic calcium ion concentration',
            'positive regulation of mitochondrial membrane permeability',
            'negative regulation of ATP biosynthetic process',
            'necrotic cell death',
        ]
M2_NODES = [
            "<http://purl.obolibrary.org/obo/GO_0051924>",
            "<http://purl.obolibrary.org/obo/GO_0007204>",
            "<http://purl.obolibrary.org/obo/GO_0035794>",
            "<http://purl.obolibrary.org/obo/GO_2001170>",
            "<http://purl.obolibrary.org/obo/GO_0070265>",
            ]

# Mechanism 3: Increased cytosolic calcium, resulting in calpain-mediated cytoskeletal damage
# -----------
M3_NODE_LABELS = [
            'regulation of calcium ion transport',
            'positive regulation of cytosolic calcium ion concentration',
            'calcium-dependent cysteine-type endopeptidase (calpain) activity',
            'microtubule severing',
            'necrotic cell death',
            ]
M3_NODES = [
            "<http://purl.obolibrary.org/obo/GO_0051924>",
            "<http://purl.obolibrary.org/obo/GO_0007204>",
            "<http://purl.obolibrary.org/obo/GO_0004198>",
            "<http://purl.obolibrary.org/obo/GO_0051013>",
            "<http://purl.obolibrary.org/obo/GO_0070265>",
            ]

# Mechanism 4: Xenobiotic-induced oxydative stress
# -----------
M4_NODE_LABELS = [
            'membrane lipid catabolic process (peroxidation)',
            'aldehyde oxidase activity',
            'oxidation-reduction process',
            'cellular response to redox state',
            'cellular response to oxidative stress',
            'regulation of oxidative stress-induced cell death',
        ]
M4_NODES = [
            "<http://purl.obolibrary.org/obo/GO_0046466>",
            "<http://purl.obolibrary.org/obo/GO_0004031>",
            "<http://purl.obolibrary.org/obo/GO_0055114>",
            "<http://purl.obolibrary.org/obo/GO_0071461>",
            "<http://purl.obolibrary.org/obo/GO_0034599>",
            "<http://purl.obolibrary.org/obo/GO_1903201>",
            ]

# Mechanism 5: Mitochondria-mediated toxicity by inhibition of electron transport chain
# -----------
M5_NODE_LABELS = [
            'negative regulation of mitochondrial electron transport, NADH to ubiquinone',
            'negative regulation of mitochondrial ATP synthesis coupled proton transport',
            'negative regulation of ATP biosynthetic process',
            'cellular response to reactive oxygen species',
            'mitochondrial DNA repair',
            ]
M5_NODES = [
            "<http://purl.obolibrary.org/obo/GO_1902957>",
            "<http://purl.obolibrary.org/obo/GO_1905707>",
            "<http://purl.obolibrary.org/obo/GO_2001170>",
            "<http://purl.obolibrary.org/obo/GO_0034614>",
            "<http://purl.obolibrary.org/obo/GO_0043504>",
            ]

# Mechanism 6: Inhibition of tissue repair by cell cycle disruption
# -----------
M6_NODE_LABELS = [
            'negative regulation of G0 to G1 transition',
            'negative regulation of cell cycle G2/M phase transition',
            'negative regulation of mitotic cell cycle',
            'positive regulation of cell cycle arrest',
            'positive regulation of apoptotic process'
        ]
M6_NODES = [
            "<http://purl.obolibrary.org/obo/GO_0070317>",
            "<http://purl.obolibrary.org/obo/GO_1902750>",
            "<http://purl.obolibrary.org/obo/GO_0045930>",
            "<http://purl.obolibrary.org/obo/GO_0071158>",
            "<http://purl.obolibrary.org/obo/GO_0043065>",
            ]

# Mechanism 7: Endoplasmic reticulum stress (by chemical or a metabolite covalently bound to proteins)
# -----------
M7_NODE_LABELS = [
            'endoplasmic reticulum unfolded protein response',
            'positive regulation of signal transduction',
            'positive regulation of protein folding',
            'positive regulation of chaperone-mediated protein folding',
            'response to endoplasmic reticulum stress',
            'positive regulation of endoplasmic reticulum stress-induced intrinsic apoptotic signaling pathway',
        ]
M7_NODES = [
            "<http://purl.obolibrary.org/obo/GO_0030968>",
            "<http://purl.obolibrary.org/obo/GO_0009967>",
            "<http://purl.obolibrary.org/obo/GO_1903334>",
            "<http://purl.obolibrary.org/obo/GO_1903646>",
            "<http://purl.obolibrary.org/obo/GO_0034976>",
            "<http://purl.obolibrary.org/obo/GO_1902237>",
            ]


# Mechanism 8: ER-mediated toxicity
# -----------
M8_NODE_LABELS = [
#            'heat shock protein binding',   # ligand to undergo necessary conformational changeu
            'protein homodimerization activity',
            'estrogen response element binding',
            'estrogen receptor activity',
            'intracellular estrogen receptor signaling pathway',
            'cellular response to estrogen stimulus',
        ]
M8_NODES = [
#            "<http://purl.obolibrary.org/obo/GO_0031072",
            "<http://purl.obolibrary.org/obo/GO_0042803>",
            "<http://purl.obolibrary.org/obo/GO_0034056>",
            "<http://purl.obolibrary.org/obo/GO_0030284>",
            "<http://purl.obolibrary.org/obo/GO_0030520>",
            "<http://purl.obolibrary.org/obo/GO_0071391>",
            ]


# Mechanism 9: AHR-mediated toxicity
# -----------
M9_NODE_LABELS = [
            'aryl hydrocarbon receptor binding',
            'protein heterodimerization activity',
            'glutathione transferase activity',
            'glucuronosyltransferase activity',
            'negative regulation of cell cycle phase transition',
        ]
M9_NODES = [
            "<http://purl.obolibrary.org/obo/GO_0017162>",
            "<http://purl.obolibrary.org/obo/GO_0046982>",
            "<http://purl.obolibrary.org/obo/GO_0004364>",
            "<http://purl.obolibrary.org/obo/GO_0015020>",
            "<http://purl.obolibrary.org/obo/GO_1901988>",
            ]


# Mechanism 10: AR-mediated toxicity
# -----------
M10_NODE_LABELS = [
#            'heat shock protein binding',   # ligand to undergo necessary conformational changeu
            'protein dimerization activity',
            'androgen receptor binding',
            'androgen receptor signaling pathway',
        ]
M10_NODES = [
#            "<http://purl.obolibrary.org/obo/GO_0031072",
            "<http://purl.obolibrary.org/obo/GO_0046983>",
            "<http://purl.obolibrary.org/obo/GO_0050681>",
            "<http://purl.obolibrary.org/obo/GO_0030521>",
            ]


# Mechanism 11: PPAR-gamma-mediated alteration of fatty-acid metabolism
# -----------
M11_NODE_LABELS = [
            'peroxisome proliferator activated receptor binding',
            'fatty acid binding',
            'positive regulation of fatty acid biosynthetic process',
            'negative regulation of fatty acid beta-oxidation', # in mitochondria
            'positive regulation of lipid storage',
            'oxidative phosphorylation uncoupler activity',
        ]
M11_NODES = [
            "<http://purl.obolibrary.org/obo/GO_0042975>",
            "<http://purl.obolibrary.org/obo/GO_0005504>",
            "<http://purl.obolibrary.org/obo/GO_0045723>",
            "<http://purl.obolibrary.org/obo/GO_0031999>",
            "<http://purl.obolibrary.org/obo/GO_0010884>",
            "<http://purl.obolibrary.org/obo/GO_0017077>",
            ]


## Mechanism 12, 13, etc...: Add your own here
## -----------
#M12_NODE_LABELS = [
#            'some label for mechanism step 1',
#            'some label for mechanism step 2',
#            'some label for mechanism step 3... etc',
#        ]
#M12_NODES = [
#            "<http://purl.obolibrary.org/obo/GO_xxxxxxx>",
#            "<http://purl.obolibrary.org/obo/GO_xxxxxxx>",
#            "<http://purl.obolibrary.org/obo/GO_xxxxxxx>",
#            ]


MECHANISMS = {
                "M1": M1_NODES,
                "M2": M2_NODES,
                "M3": M3_NODES,
                "M4": M4_NODES,
                "M5": M5_NODES,
                "M6": M6_NODES,
                "M7": M7_NODES,
                "M8": M8_NODES,
                "M9": M9_NODES,
                "M10": M10_NODES,
                "M11": M11_NODES,
#                "M12": M12_NODES,
             }

MECHANISM_DESCRIPTIONS = {
                "M1": M1_NODE_LABELS,
                "M2": M2_NODE_LABELS,
                "M3": M3_NODE_LABELS,
                "M4": M4_NODE_LABELS,
                "M5": M5_NODE_LABELS,
                "M6": M6_NODE_LABELS,
                "M7": M7_NODE_LABELS,
                "M8": M8_NODE_LABELS,
                "M9": M9_NODE_LABELS,
                "M10": M10_NODE_LABELS,
                "M11": M11_NODE_LABELS,
#                "M12": M12_NODE_LABELS,
             }

MECHANISM_TOPS = {
                "M1": "Triggering of caspase-mediated apoptosis via release of cytochrome C",
                "M2": "ATP depletion from calcium homeostasis disruption, resulting in necrosis",
                "M3": "Increased cytosolic calcium, resulting in calpain-mediated cytoskeletal damage",
                "M4": "Xenobiotic-induced oxydative stress",
                "M5": "Mitochondria-mediated toxicity by inhibition of electron transport chain",
                "M6": "Inhibition of tissue repair by cell cycle disruption",
                "M7": "Endoplasmic reticulum stress (by chemical or a metabolite covalently bound to proteins)",
                "M8": "Triggering of estrogen receptor (ER) activity",
                "M9": "Triggering of aryl hydrocarbon receptor (AHR) activity",
                "M10": "Triggering of androgen receptor (AR) activity",
                "M11": "Triggering of peroxisome proliferator-activated receptor gamma (PPAR-gamma) alteration of fatty acid metabolism",
#                "M12": "Some new mechanism",
             }

