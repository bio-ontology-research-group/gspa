"""Fixed, cited biological checks for real GSPA whole-proteome predictions.

A failed biology check is a result to report, never a reason to alter a rule.
All rules are defined from published biological literature and official enzyme
nomenclature before inspecting predictor outcomes.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

V3_DIR = Path(__file__).resolve().parents[1]
if str(V3_DIR) not in sys.path:
    sys.path.insert(0, str(V3_DIR))

from biology.extract import (  # noqa: E402
    BacterialTaxonConstraintIndex,
     canonicalize_ec4,
    ec3_prefix,
    load_bacterial_taxon_constraints,
    load_ec_audit_entries,
    load_price148_targets,
    normalize_ec,
    unpack_ec_record,
)

# ---------------------------------------------------------------------------
# Genome Exemption Policy (Obligate Endosymbionts & Host-Restricted Parasites)
# ---------------------------------------------------------------------------

EXEMPTION_POLICY_CITATION = (
    "McCutcheon JP & Moran NA (2012) Extreme genome reduction in symbiotic "
    "bacteria, Nat. Rev. Microbiol. 10(1):13-26 (doi:10.1038/nrmicro2670); "
    "Moran NA & Bennett GM (2014) The tiniest tiny genomes, Annu. Rev. "
    "Microbiol. 68:195-215 (doi:10.1146/annurev-micro-091213-112901); "
    "Fraser CM et al. (1995) The minimal gene complement of Mycoplasma "
    "genitalium, Science 270(5235):397-403 (doi:10.1126/science.270.5235.397)."
)

EXEMPT_ENDOSYMBIONT_GENERA: Dict[str, str] = {
    "Buchnera": "Obligate aphid nutritional endosymbiont with reduced metabolic repertoire (Moran & Bennett 2014).",
    "Wigglesworthia": "Obligate tsetse fly endosymbiont with reduced genome (McCutcheon & Moran 2012).",
    "Blochmannia": "Obligate carpenter ant endosymbiont with reduced genome (McCutcheon & Moran 2012).",
    "Carsonella": "Extreme psyllid nutritional endosymbiont lacking subsets of aaRS and replication genes (Moran & Bennett 2014).",
    "Tremblaya": "Extreme mealybug endosymbiont lacking aminoacyl-tRNA synthetases (McCutcheon & Moran 2012).",
    "Nasuia": "Extreme leafhopper endosymbiont with 112 kb genome (Moran & Bennett 2014).",
    "Sulcia": "Obligate sap-feeding insect endosymbiont with reduced genome (McCutcheon & Moran 2012).",
    "Hodgkinia": "Extreme cicada endosymbiont with 144 kb genome (McCutcheon & Moran 2012).",
    "Zinderia": "Extreme spittlebug endosymbiont with 208 kb genome (McCutcheon & Moran 2012).",
    "Moranella": "Nested mealybug endosymbiont with reduced genome (McCutcheon & Moran 2012).",
    "Mycoplasma": "Obligate host-associated Mollicutes parasite with reduced metabolic pathways (Fraser et al. 1995).",
    "Mycoplasmoides": "Obligate host-associated Mollicutes parasite with reduced metabolic pathways (Fraser et al. 1995).",
    "Mesomycoplasma": "Obligate host-associated Mollicutes parasite with reduced metabolic pathways (Fraser et al. 1995).",
    "Ureaplasma": "Obligate host-associated Mollicutes lacking canonical F0F1-ATP synthase mode (Fraser et al. 1995).",
    "Phytoplasma": "Obligate phloem-limited plant parasite lacking F0F1-ATP synthase (McCutcheon & Moran 2012).",
    "Chlamydia": "Obligate intracellular bacterium (energy parasite) with reduced central metabolism (McCutcheon & Moran 2012).",
    "Chlamydophila": "Obligate intracellular bacterium with reduced central metabolism (McCutcheon & Moran 2012).",
    "Rickettsia": "Obligate intracellular alphaproteobacterium lacking glycolytic trunk enzymes (McCutcheon & Moran 2012).",
    "Orientia": "Obligate intracellular alphaproteobacterium with reduced central metabolism (McCutcheon & Moran 2012).",
    "Wolbachia": "Obligate intracellular alphaproteobacterium lacking glycolytic trunk enzymes (McCutcheon & Moran 2012).",
    "Ehrlichia": "Obligate intracellular alphaproteobacterium with reduced metabolism (McCutcheon & Moran 2012).",
    "Anaplasma": "Obligate intracellular alphaproteobacterium with reduced metabolism (McCutcheon & Moran 2012).",
}

EXEMPT_ENDOSYMBIONT_TAXIDS: Dict[str, str] = {
    "9": "Buchnera aphidicola (obligate aphid endosymbiont; Moran & Bennett 2014)",
    "51229": "Wigglesworthia glossinidia (obligate tsetse endosymbiont; McCutcheon & Moran 2012)",
    "2093": "Mycoplasma (host-restricted Mollicutes with reduced genome; Fraser et al. 1995)",
    "2097": "Mycoplasmoides genitalium (minimal cellular genome; Fraser et al. 1995)",
    "810": "Chlamydia (obligate intracellular bacterium; McCutcheon & Moran 2012)",
    "780": "Rickettsia (obligate intracellular bacterium; McCutcheon & Moran 2012)",
    "953": "Wolbachia (obligate intracellular endosymbiont; McCutcheon & Moran 2012)",
}

# All 22 Price-148 genomes are free-living heterotrophic Pseudomonadota (Price et al. 2018, Nature 557:503-509).
PRICE22_GENOME_METADATA: Dict[str, Dict[str, Any]] = {
    "GCA_000006965.1": {"taxid": "266834", "organism": "Sinorhizobium meliloti 1021", "lifestyle": "free-living / facultative root-nodule symbiont"},
    "GCA_003443575.1": {"taxid": "221822", "organism": "Phaeobacter inhibens DSM 16374", "lifestyle": "free-living marine heterotroph"},
    "GCA_004402975.1": {"taxid": "1211383", "organism": "Paraburkholderia dipogonis", "lifestyle": "free-living soil/nodule heterotroph"},
    "GCF_000015245.1": {"taxid": "326297", "organism": "Shewanella amazonensis SB2B", "lifestyle": "free-living aquatic heterotroph"},
    "GCF_000016065.1": {"taxid": "323850", "organism": "Shewanella loihica PV-4", "lifestyle": "free-living marine bacterium"},
    "GCF_000018145.1": {"taxid": "398580", "organism": "Dinoroseobacter shibae DFL 12", "lifestyle": "free-living marine heterotroph"},
    "GCF_000020125.1": {"taxid": "398527", "organism": "Paraburkholderia phytofirmans PsJN", "lifestyle": "free-living plant-associated bacterium"},
    "GCF_000146165.2": {"taxid": "211586", "organism": "Shewanella oneidensis MR-1", "lifestyle": "free-living facultative anaerobe"},
    "GCF_000166295.1": {"taxid": "225937", "organism": "Marinobacter adhaerens HP15", "lifestyle": "free-living marine heterotroph"},
    "GCF_000203935.1": {"taxid": "94122", "organism": "Shewanella sp. ANA-3", "lifestyle": "free-living aquatic heterotroph"},
    "GCF_000236665.1": {"taxid": "640081", "organism": "Azospira oryzae PS", "lifestyle": "free-living perchlorate-reducing bacterium"},
    "GCF_000237365.1": {"taxid": "1064539", "organism": "Azospirillum baldaniorum Sp245", "lifestyle": "free-living rhizosphere diazotroph"},
    "GCF_000327065.1": {"taxid": "644801", "organism": "Stutzerimonas stutzeri RCH2", "lifestyle": "free-living denitrifying heterotroph"},
    "GCF_000698265.1": {"taxid": "321846", "organism": "Pseudomonas simiae WCS417", "lifestyle": "free-living rhizosphere heterotroph"},
    "GCF_000832305.1": {"taxid": "68895", "organism": "Cupriavidus basilensis 4G11", "lifestyle": "free-living soil heterotroph"},
    "GCF_001307155.1": {"taxid": "294", "organism": "Pseudomonas fluorescens FW300-N1B4", "lifestyle": "free-living rhizosphere heterotroph"},
    "GCF_001307275.1": {"taxid": "294", "organism": "Pseudomonas fluorescens FW300-N2E2", "lifestyle": "free-living rhizosphere heterotroph"},
    "GCF_001440325.1": {"taxid": "294", "organism": "Pseudomonas fluorescens GW456-L13", "lifestyle": "free-living rhizosphere heterotroph"},
    "GCF_001623525.1": {"taxid": "294", "organism": "Pseudomonas fluorescens FW300-N2C3", "lifestyle": "free-living rhizosphere heterotroph"},
    "GCF_001625455.1": {"taxid": "294", "organism": "Pseudomonas fluorescens FW300-N2E3", "lifestyle": "free-living rhizosphere heterotroph"},
    "GCF_001633105.1": {"taxid": "1813946", "organism": "Acidovorax sp. GW101-3H11", "lifestyle": "free-living groundwater heterotroph"},
    "GCF_002797435.1": {"taxid": "93064", "organism": "Sphingomonas koreensis DSMZ 15582", "lifestyle": "free-living aquatic/soil heterotroph"},
}


# ---------------------------------------------------------------------------
# Check (a): Translation — 20 Aminoacyl-tRNA Charging Routes (EC 6.1.1.1-6.1.1.22)
# ---------------------------------------------------------------------------
# Note: EC 6.1.1.8 is a deleted ENZYME entry; EC 6.1.1.13 is D-alanine--poly(phosphoribitol)
# ligase (teichoic acid biosynthesis, not a proteinogenic tRNA synthetase).
# The remaining 20 entries in 6.1.1.1-6.1.1.22 cover the 20 standard amino acids,
# with GatCAB/GatDE transamidation routes 6.3.5.7 (Gln) and 6.3.5.6 (Asn) accepted
# as valid physiological alternatives where direct GlnRS (6.1.1.18) or AsnRS (6.1.1.22)
# is absent (Sheppard et al. 2008). Non-discriminating AspRS (6.1.1.23) and GluRS (6.1.1.24)
# are also accepted alongside 6.1.1.12 and 6.1.1.17.

TRANSLATION_RULE_CITATION = (
    "Ibba M & Soll D (2000) Aminoacyl-tRNA synthesis, Annu. Rev. Biochem. "
    "69:617-650 (doi:10.1146/annurev.biochem.69.1.617); "
    "Sheppard K et al. (2008) From one amino acid to another: tRNA-dependent "
    "amino acid biosynthesis, Nucleic Acids Res. 36(6):1813-1825 "
    "(doi:10.1093/nar/gkn015); "
    "IUBMB Enzyme Nomenclature / ExPASy ENZYME subclass EC 6.1.1 and EC 6.3.5.6/6.3.5.7."
)

TRANSLATION_CHARGING_RULES: List[Dict[str, Any]] = [
    {
        "amino_acid": "Tyr",
        "name": "tyrosine--tRNA ligase",
        "direct_ecs": ["6.1.1.1"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; IUBMB EC 6.1.1.1",
    },
    {
        "amino_acid": "Trp",
        "name": "tryptophan--tRNA ligase",
        "direct_ecs": ["6.1.1.2"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; IUBMB EC 6.1.1.2",
    },
    {
        "amino_acid": "Thr",
        "name": "threonine--tRNA ligase",
        "direct_ecs": ["6.1.1.3"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; IUBMB EC 6.1.1.3",
    },
    {
        "amino_acid": "Leu",
        "name": "leucine--tRNA ligase",
        "direct_ecs": ["6.1.1.4"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; IUBMB EC 6.1.1.4",
    },
    {
        "amino_acid": "Ile",
        "name": "isoleucine--tRNA ligase",
        "direct_ecs": ["6.1.1.5"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; IUBMB EC 6.1.1.5",
    },
    {
        "amino_acid": "Lys",
        "name": "lysine--tRNA ligase",
        "direct_ecs": ["6.1.1.6"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; IUBMB EC 6.1.1.6",
    },
    {
        "amino_acid": "Ala",
        "name": "alanine--tRNA ligase",
        "direct_ecs": ["6.1.1.7"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; IUBMB EC 6.1.1.7",
    },
    {
        "amino_acid": "Val",
        "name": "valine--tRNA ligase",
        "direct_ecs": ["6.1.1.9"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; IUBMB EC 6.1.1.9",
    },
    {
        "amino_acid": "Met",
        "name": "methionine--tRNA ligase",
        "direct_ecs": ["6.1.1.10"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; IUBMB EC 6.1.1.10",
    },
    {
        "amino_acid": "Ser",
        "name": "serine--tRNA ligase",
        "direct_ecs": ["6.1.1.11"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; IUBMB EC 6.1.1.11",
    },
    {
        "amino_acid": "Asp",
        "name": "aspartate--tRNA ligase (discriminating or non-discriminating)",
        "direct_ecs": ["6.1.1.12", "6.1.1.23"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; Sheppard et al. (2008) NAR 36:1813-1825; IUBMB EC 6.1.1.12 / 6.1.1.23",
    },
    {
        "amino_acid": "Gly",
        "name": "glycine--tRNA ligase",
        "direct_ecs": ["6.1.1.14"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; IUBMB EC 6.1.1.14",
    },
    {
        "amino_acid": "Pro",
        "name": "proline--tRNA ligase",
        "direct_ecs": ["6.1.1.15"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; IUBMB EC 6.1.1.15",
    },
    {
        "amino_acid": "Cys",
        "name": "cysteine--tRNA ligase",
        "direct_ecs": ["6.1.1.16"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; IUBMB EC 6.1.1.16",
    },
    {
        "amino_acid": "Glu",
        "name": "glutamate--tRNA ligase (discriminating or non-discriminating)",
        "direct_ecs": ["6.1.1.17", "6.1.1.24"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; Sheppard et al. (2008) NAR 36:1813-1825; IUBMB EC 6.1.1.17 / 6.1.1.24",
    },
    {
        "amino_acid": "Gln",
        "name": "glutamine--tRNA ligase (direct 6.1.1.18) or GatCAB/GatDE Glu-tRNA(Gln) amidotransferase (indirect 6.3.5.7)",
        "direct_ecs": ["6.1.1.18"],
        "indirect_ecs": ["6.3.5.7"],
        "source": "Sheppard et al. (2008) Nucleic Acids Res. 36(6):1813-1825; IUBMB EC 6.1.1.18 / 6.3.5.7",
    },
    {
        "amino_acid": "Arg",
        "name": "arginine--tRNA ligase",
        "direct_ecs": ["6.1.1.19"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; IUBMB EC 6.1.1.19",
    },
    {
        "amino_acid": "Phe",
        "name": "phenylalanine--tRNA ligase",
        "direct_ecs": ["6.1.1.20"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; IUBMB EC 6.1.1.20",
    },
    {
        "amino_acid": "His",
        "name": "histidine--tRNA ligase",
        "direct_ecs": ["6.1.1.21"],
        "indirect_ecs": [],
        "source": "Ibba & Soll (2000) Annu. Rev. Biochem. 69:617-650; IUBMB EC 6.1.1.21",
    },
    {
        "amino_acid": "Asn",
        "name": "asparagine--tRNA ligase (direct 6.1.1.22) or GatCAB Asp-tRNA(Asn) amidotransferase (indirect 6.3.5.6)",
        "direct_ecs": ["6.1.1.22"],
        "indirect_ecs": ["6.3.5.6"],
        "source": "Sheppard et al. (2008) Nucleic Acids Res. 36(6):1813-1825; IUBMB EC 6.1.1.22 / 6.3.5.6",
    },
]


# ---------------------------------------------------------------------------
# Check (b): Core Information-Processing and Central Metabolism ECs
# ---------------------------------------------------------------------------

CORE_METABOLISM_RULE_CITATION = (
    "Gil R, Silva FJ, Pereto J & Moya A (2004) Determination of the core of a "
    "minimal bacterial gene set, Microbiol. Mol. Biol. Rev. 68(3):518-537 "
    "(doi:10.1128/MMBR.68.3.518-537.2004); "
    "Hutchison CA III et al. (2016) Design and synthesis of a minimal bacterial "
    "genome, Science 351(6280):aad6253 (doi:10.1126/science.aad6253); "
    "Wilkinson A, Day J & Bowater R (2001) Bacterial DNA ligases, Mol. Microbiol. "
    "40(6):1241-1248 (doi:10.1046/j.1365-2958.2001.02479.x); "
    "Romano AH & Conway T (1996) Evolution of carbohydrate metabolic pathways, "
    "Res. Microbiol. 147(6-7):448-455; "
    "Hove-Jensen B et al. (2017) Phosphoribosyl diphosphate (PRPP): biosynthesis, "
    "enzymology, utilization, and metabolic significance, Microbiol. Mol. Biol. "
    "Rev. 81(1):e00040-16; "
    "IUBMB Enzyme Nomenclature / ExPASy ENZYME."
)

CORE_BACTERIAL_EC_RULES: List[Dict[str, Any]] = [
    {
        "role_id": "rna_polymerase",
        "name": "DNA-directed RNA polymerase",
        "category": "information_processing",
        "accepted_ecs": ["2.7.7.6"],
        "source": "Gil et al. (2004) MMBR 68:518-537; Hutchison et al. (2016) Science 351:aad6253; IUBMB EC 2.7.7.6",
    },
    {
        "role_id": "dna_polymerase",
        "name": "DNA-directed DNA polymerase",
        "category": "information_processing",
        "accepted_ecs": ["2.7.7.7"],
        "source": "Gil et al. (2004) MMBR 68:518-537; Hutchison et al. (2016) Science 351:aad6253; IUBMB EC 2.7.7.7",
    },
    {
        "role_id": "dna_topoisomerase_i",
        "name": "DNA topoisomerase type I",
        "category": "information_processing",
        "accepted_ecs": ["5.6.2.1", "5.99.1.2"],
        "source": "Gil et al. (2004) MMBR 68:518-537; Hutchison et al. (2016) Science 351:aad6253; IUBMB EC 5.99.1.2 transferred to 5.6.2.1",
    },
    {
        "role_id": "dna_gyrase_topoisomerase_ii",
        "name": "DNA topoisomerase (ATP-hydrolyzing) / DNA gyrase & topoisomerase IV",
        "category": "information_processing",
        "accepted_ecs": ["5.6.2.2", "5.99.1.3"],
        "source": "Gil et al. (2004) MMBR 68:518-537; Hutchison et al. (2016) Science 351:aad6253; IUBMB EC 5.99.1.3 transferred to 5.6.2.2",
    },
    {
        "role_id": "dna_ligase",
        "name": "DNA ligase (NAD+ or ATP)",
        "category": "information_processing",
        "accepted_ecs": ["6.5.1.2", "6.5.1.1"],
        "source": "Wilkinson et al. (2001) Mol. Microbiol. 40(6):1241-1248; Gil et al. (2004) MMBR 68:518-537; IUBMB EC 6.5.1.2 / 6.5.1.1",
    },
    {
        "role_id": "triose_phosphate_isomerase",
        "name": "triose-phosphate isomerase",
        "category": "central_metabolism",
        "accepted_ecs": ["5.3.1.1"],
        "source": "Gil et al. (2004) MMBR 68:518-537; Romano & Conway (1996) Res. Microbiol. 147:448-455; IUBMB EC 5.3.1.1",
    },
    {
        "role_id": "phosphoglycerate_kinase",
        "name": "phosphoglycerate kinase",
        "category": "central_metabolism",
        "accepted_ecs": ["2.7.2.3"],
        "source": "Gil et al. (2004) MMBR 68:518-537; Hutchison et al. (2016) Science 351:aad6253; IUBMB EC 2.7.2.3",
    },
    {
        "role_id": "enolase",
        "name": "phosphopyruvate hydratase (enolase)",
        "category": "central_metabolism",
        "accepted_ecs": ["4.2.1.11"],
        "source": "Gil et al. (2004) MMBR 68:518-537; Hutchison et al. (2016) Science 351:aad6253; IUBMB EC 4.2.1.11",
    },
    {
        "role_id": "prpp_synthetase",
        "name": "ribose-phosphate diphosphokinase (PRPP synthetase)",
        "category": "central_metabolism",
        "accepted_ecs": ["2.7.6.1"],
        "source": "Gil et al. (2004) MMBR 68:518-537; Hove-Jensen et al. (2017) MMBR 81(1):e00040-16; IUBMB EC 2.7.6.1",
    },
    {
        "role_id": "f0f1_atp_synthase",
        "name": "H(+)-transporting two-sector ATPase (F0F1-ATP synthase)",
        "category": "central_metabolism",
        "accepted_ecs": ["7.1.2.2", "3.6.3.14"],
        "source": "Gil et al. (2004) MMBR 68:518-537; Hutchison et al. (2016) Science 351:aad6253; IUBMB EC 3.6.3.14 transferred to 7.1.2.2",
    },
]

TAXON_CONSISTENCY_CITATION = (
    "Deegan (nee Clark) JI, Dimmer EC & Mungall CJ (2010) Formalization of "
    "taxon-based constraints to detect inconsistencies in annotation and "
    "ontology development, BMC Bioinformatics 11:530 "
    "(doi:10.1186/1471-2105-11-530); bundled GSPA resources "
    "gspa-core/src/main/resources/taxon-constraints/go-taxon-constraints.tsv and "
    "ncbi-taxon-hierarchy.tsv with benchmark/whole_genome_cafa/data/reference/ec2go.txt."
)

TRUTH_AGREEMENT_CITATION = (
    "Price MN et al. (2018) Mutant phenotypes for thousands of bacterial genes "
    "of unknown function, Nature 557(7706):503-509 "
    "(doi:10.1038/s41586-018-0124-0); benchmark/manifests/suite6_price149_manifest.tsv."
)

MUTUAL_EXCLUSIVITY_OMISSION_NOTE = (
    "Skipped per specification (checks c and d take priority): under IUBMB "
    "Enzyme Nomenclature Rule 2 (Webb EC, 1992, Enzyme Nomenclature 1992, "
    "Academic Press; McDonald AG & Tipton KF, 2023, FEBS J. 290:2214-2231), a "
    "single EC number designates both forward and reverse directions of a "
    "reversible reaction, and separate EC numbers for opposing physiological "
    "steps (such as kinase 2.7.1.11 vs. hydrolase 3.1.3.11) represent distinct "
    "chemical reactions that can legitimately co-occur on bifunctional enzymes."
)

STRATA_ORDER = [
    "Stratum_A_Core_Biomass",
    "Stratum_B_Catabolic_Assay_Specific",
    "Stratum_C_Secondary_Peripheral",
]


def check_genome_exemption(
    assembly: str,
    taxon_id: Optional[str] = None,
    organism: Optional[str] = None,
    exempt_assemblies: Optional[Iterable[str]] = None,
    exempt_taxids: Optional[Iterable[str]] = None,
) -> Dict[str, Any]:
    """Determine whether a bacterial genome is exempt from free-living completeness checks."""
    if exempt_assemblies and assembly in set(exempt_assemblies):
        return {
            "exempt": True,
            "reason": f"Assembly {assembly} explicitly marked exempt by caller ({EXEMPTION_POLICY_CITATION})",
        }

    meta = PRICE22_GENOME_METADATA.get(assembly, {})
    eff_taxid = str(taxon_id or meta.get("taxid") or "").removeprefix("NCBITaxon_").strip()
    eff_org = str(organism or meta.get("organism") or "").strip()

    if exempt_taxids:
        norm_exempt = {str(t).removeprefix("NCBITaxon_").strip() for t in exempt_taxids}
        if eff_taxid and eff_taxid in norm_exempt:
            return {
                "exempt": True,
                "reason": f"Taxon {eff_taxid} explicitly marked exempt by caller ({EXEMPTION_POLICY_CITATION})",
            }

    if eff_taxid and eff_taxid in EXEMPT_ENDOSYMBIONT_TAXIDS:
        return {
            "exempt": True,
            "reason": EXEMPT_ENDOSYMBIONT_TAXIDS[eff_taxid],
        }

    if eff_org:
        tokens = eff_org.replace("Candidatus ", "").split()
        genus = tokens[0] if tokens else ""
        if genus in EXEMPT_ENDOSYMBIONT_GENERA:
            return {
                "exempt": True,
                "reason": EXEMPT_ENDOSYMBIONT_GENERA[genus],
            }

    return {
        "exempt": False,
        "reason": None,
        "organism": eff_org or None,
        "taxid": eff_taxid or None,
        "lifestyle": meta.get("lifestyle", "free-living bacterium"),
    }


def _best_among_ecs(
    ec_table: Dict[str, Any],
    candidate_ecs: Iterable[str],
    ec_entries: Optional[Dict[str, Any]] = None,
) -> Optional[Dict[str, Any]]:
    expanded: Set[str] = set()
    for ec in candidate_ecs:
        norm = normalize_ec(ec)
        expanded.add(norm)
        if ec_entries:
            expanded.add(canonicalize_ec4(norm, ec_entries))

    best: Optional[Dict[str, Any]] = None
    for ec in sorted(expanded):
        raw_rec = ec_table.get(ec)
        if raw_rec is None:
            continue
        rec = unpack_ec_record(raw_rec)
        entry = {
            "ec": ec,
            "best_score": rec["best_score"],
            "best_protein": rec["best_protein"],
            "n_rows": rec["n_rows"],
        }
        if best is None:
            best = entry
        else:
            if (
                entry["best_score"] > best["best_score"]
                or (
                    entry["best_score"] == best["best_score"]
                    and entry["best_protein"] < best["best_protein"]
                )
            ):
                best = entry
    return best


def evaluate_translation_arm(
    arm_summary: Dict[str, Any],
    level: str = "called",
    ec_entries: Optional[Dict[str, Any]] = None,
    exemption: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Check all 20 aminoacyl-tRNA charging routes for one arm at 'called' or 'retained' level."""
    if level not in ("called", "retained"):
        raise ValueError(f"Unsupported level: {level!r}")
    ec4_table = arm_summary.get("ec4", {}).get(level, {})
    is_exempt = bool(exemption and exemption.get("exempt"))

    routes: Dict[str, Any] = {}
    satisfied_count = 0
    direct_only_count = 0
    indirect_gat_count = 0
    both_routes_count = 0
    missing_amino_acids: List[str] = []

    for rule in TRANSLATION_CHARGING_RULES:
        aa = rule["amino_acid"]
        direct_hit = _best_among_ecs(ec4_table, rule["direct_ecs"], ec_entries=ec_entries)
        indirect_hit = (
            _best_among_ecs(ec4_table, rule["indirect_ecs"], ec_entries=ec_entries)
            if rule["indirect_ecs"]
            else None
        )

        if direct_hit and indirect_hit:
            route_mode = "both"
            satisfied = True
            both_routes_count += 1
            chosen = (
                direct_hit
                if direct_hit["best_score"] >= indirect_hit["best_score"]
                else indirect_hit
            )
        elif direct_hit:
            route_mode = "direct"
            satisfied = True
            direct_only_count += 1
            chosen = direct_hit
        elif indirect_hit:
            route_mode = "indirect_gat"
            satisfied = True
            indirect_gat_count += 1
            chosen = indirect_hit
        else:
            route_mode = "missing"
            satisfied = False
            chosen = None
            missing_amino_acids.append(aa)

        if satisfied:
            satisfied_count += 1

        routes[aa] = {
            "amino_acid": aa,
            "name": rule["name"],
            "satisfied": satisfied,
            "route_mode": route_mode,
            "best_ec": chosen["ec"] if chosen else None,
            "best_score": chosen["best_score"] if chosen else None,
            "best_protein": chosen["best_protein"] if chosen else None,
            "direct_hit": direct_hit,
            "indirect_hit": indirect_hit,
            "source": rule["source"],
        }

    total_required = len(TRANSLATION_CHARGING_RULES)
    passed = True if is_exempt else (satisfied_count == total_required)

    return {
        "level": level,
        "exempt": is_exempt,
        "exemption_reason": exemption.get("reason") if exemption else None,
        "total_amino_acids": total_required,
        "satisfied_amino_acids": satisfied_count,
        "missing_count": total_required - satisfied_count,
        "missing_amino_acids": missing_amino_acids,
        "direct_only_count": direct_only_count,
        "indirect_gat_count": indirect_gat_count,
        "both_routes_count": both_routes_count,
        "passed": passed,
        "routes": routes,
    }


def evaluate_core_metabolism_arm(
    arm_summary: Dict[str, Any],
    level: str = "called",
    ec_entries: Optional[Dict[str, Any]] = None,
    exemption: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Check presence of the 10 near-universal bacterial information-processing and central-metabolism ECs."""
    if level not in ("called", "retained"):
        raise ValueError(f"Unsupported level: {level!r}")
    ec4_table = arm_summary.get("ec4", {}).get(level, {})
    is_exempt = bool(exemption and exemption.get("exempt"))

    roles: Dict[str, Any] = {}
    present_count = 0
    missing_roles: List[str] = []

    for rule in CORE_BACTERIAL_EC_RULES:
        role_id = rule["role_id"]
        hit = _best_among_ecs(ec4_table, rule["accepted_ecs"], ec_entries=ec_entries)
        present = hit is not None
        if present:
            present_count += 1
        else:
            missing_roles.append(role_id)

        roles[role_id] = {
            "role_id": role_id,
            "name": rule["name"],
            "category": rule["category"],
            "accepted_ecs": list(rule["accepted_ecs"]),
            "present": present,
            "best_ec": hit["ec"] if hit else None,
            "best_score": hit["best_score"] if hit else None,
            "best_protein": hit["best_protein"] if hit else None,
            "source": rule["source"],
        }

    total_roles = len(CORE_BACTERIAL_EC_RULES)
    passed = True if is_exempt else (present_count == total_roles)

    return {
        "level": level,
        "exempt": is_exempt,
        "exemption_reason": exemption.get("reason") if exemption else None,
        "total_roles": total_roles,
        "present_roles": present_count,
        "missing_count": total_roles - present_count,
        "missing_roles": missing_roles,
        "passed": passed,
        "roles": roles,
    }


def evaluate_taxon_consistency_arm(
    arm_summary: Dict[str, Any],
    taxon_constraints: Optional[BacterialTaxonConstraintIndex] = None,
) -> Dict[str, Any]:
    """Evaluate bacterial taxon constraint violations and inspect cell enforcement configuration."""
    tc_block = arm_summary.get("taxon_consistency", {})
    config = arm_summary.get("config", {})

    # Re-verify EC violations against bundled taxon constraints from ec4 tables as well
    ec4_ret = arm_summary.get("ec4", {}).get("retained", {})
    ec4_cal = arm_summary.get("ec4", {}).get("called", {})

    violating_ecs_ret = dict(tc_block.get("violating_ecs_retained", {}))
    violating_ecs_cal = dict(tc_block.get("violating_ecs_called", {}))

    if taxon_constraints is not None:
        for ec, raw_rec in ec4_ret.items():
            fgos = taxon_constraints.forbidden_gos_for_ec(ec)
            if fgos and ec not in violating_ecs_ret:
                rec = unpack_ec_record(raw_rec)
                violating_ecs_ret[ec] = {
                    "best_score": rec["best_score"],
                    "best_protein": rec["best_protein"],
                    "n_rows": rec["n_rows"],
                    "forbidden_go_terms": fgos,
                }
        for ec, raw_rec in ec4_cal.items():
            fgos = taxon_constraints.forbidden_gos_for_ec(ec)
            if fgos and ec not in violating_ecs_cal:
                rec = unpack_ec_record(raw_rec)
                violating_ecs_cal[ec] = {
                    "best_score": rec["best_score"],
                    "best_protein": rec["best_protein"],
                    "n_rows": rec["n_rows"],
                    "forbidden_go_terms": fgos,
                }

    violating_gos_ret = dict(tc_block.get("violating_gos_retained", {}))
    violating_gos_cal = dict(tc_block.get("violating_gos_called", {}))

    rows_ret = int(
        tc_block.get(
            "violation_rows_retained",
            sum(v["n_rows"] for v in violating_ecs_ret.values())
            + sum(v["n_rows"] for v in violating_gos_ret.values()),
        )
    )
    rows_cal = int(
        tc_block.get(
            "violation_rows_called",
            sum(v["n_rows"] for v in violating_ecs_cal.values())
            + sum(v["n_rows"] for v in violating_gos_cal.values()),
        )
    )

    enforce_stage2 = config.get("enforce_stage2")
    symbolic_status = config.get("symbolic_enforcement_status")
    enforcement_on = bool(enforce_stage2 is True and symbolic_status != "skipped")

    return {
        "target_taxon": tc_block.get("target_taxon", "NCBITaxon_2"),
        "genome_taxon": tc_block.get("genome_taxon"),
        "genome_taxon_in_bundled_hierarchy": bool(
            tc_block.get("genome_taxon_in_bundled_hierarchy", False)
        ),
        "enforcement_on": enforcement_on,
        "enforce_stage2": enforce_stage2,
        "enforce_mode": config.get("enforce_mode"),
        "symbolic_enforcement_status": symbolic_status,
        "priors_enabled": config.get("priors_enabled"),
        "priors_fired_rows": int(arm_summary.get("priors_fired_rows", 0)),
        "violation_rows_called": rows_cal,
        "violation_rows_retained": rows_ret,
        "violating_distinct_ecs_called": len(violating_ecs_cal),
        "violating_distinct_ecs_retained": len(violating_ecs_ret),
        "violating_distinct_gos_called": len(violating_gos_cal),
        "violating_distinct_gos_retained": len(violating_gos_ret),
        "violating_ecs_called": violating_ecs_cal,
        "violating_ecs_retained": violating_ecs_ret,
        "violating_gos_called": violating_gos_cal,
        "violating_gos_retained": violating_gos_ret,
    }


def _score_single_truth_protein(
    prot_summary: Dict[str, Any],
    level: str,
    ec_entries: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Score one Price-148 target protein at 4-level and 3-level EC.

    Unconfirmed labels are never treated as negatives: only the experimentally
    assigned Price-148 truth proteins enter this evaluation, and additional
    candidate ECs on other proteins are never counted as false positives.
    """
    gt_ec4s = {normalize_ec(e) for e in prot_summary.get("ground_truth_ecs", []) if e}
    gt_canon_ec4s = {canonicalize_ec4(e, ec_entries) for e in gt_ec4s}
    gt_ec3s = {
        p
        for e in gt_ec4s
        for p in [ec3_prefix(e)]
        if p is not None
    }

    block = prot_summary.get(level, {})
    top_ec4 = block.get("top_ec4")
    top_ec4_tied = [normalize_ec(e) for e in block.get("top_ec4_tied", []) if e]
    top_ec4_score = block.get("top_ec4_score")

    canon_top_ec4 = block.get("canonical_top_ec4")
    canon_tied = [normalize_ec(e) for e in block.get("canonical_top_ec4_tied", []) if e]
    if not canon_tied and top_ec4_tied:
        canon_tied = sorted({canonicalize_ec4(e, ec_entries) for e in top_ec4_tied})
        canon_top_ec4 = canon_tied[0] if canon_tied else None

    top_ec3 = block.get("top_ec3")
    top_ec3_tied = [e for e in block.get("top_ec3_tied", []) if e]
    top_ec3_score = block.get("top_ec3_score")

    # Exact 4-level match (lexical top-1 and fractional tie-aware credit)
    exact_ec4_match = bool(top_ec4 is not None and normalize_ec(top_ec4) in gt_ec4s)
    exact_ec4_tie_hits = sum(1 for e in top_ec4_tied if e in gt_ec4s)
    exact_ec4_fractional = (
        exact_ec4_tie_hits / len(top_ec4_tied) if top_ec4_tied else 0.0
    )

    # Canonical 4-level match (resolving unambiguous single-destination ENZYME transfers)
    canonical_ec4_match = bool(
        canon_top_ec4 is not None and normalize_ec(canon_top_ec4) in gt_canon_ec4s
    )
    canon_ec4_tie_hits = sum(1 for e in canon_tied if e in gt_canon_ec4s)
    canonical_ec4_fractional = (
        canon_ec4_tie_hits / len(canon_tied) if canon_tied else 0.0
    )

    # 3-level prefix match
    ec3_match = bool(top_ec3 is not None and top_ec3 in gt_ec3s)
    ec3_tie_hits = sum(1 for p in top_ec3_tied if p in gt_ec3s)
    ec3_fractional = ec3_tie_hits / len(top_ec3_tied) if top_ec3_tied else 0.0

    return {
        "price_entry": prot_summary.get("price_entry"),
        "proteome_protein_id": prot_summary.get("proteome_protein_id"),
        "metabolic_stratum": prot_summary.get("metabolic_stratum"),
        "ground_truth_ecs": sorted(gt_ec4s),
        "ground_truth_ec3s": sorted(gt_ec3s),
        "has_prediction": top_ec4 is not None,
        "top_ec4": top_ec4,
        "top_ec4_tied": top_ec4_tied,
        "top_ec4_score": top_ec4_score,
        "exact_ec4_match": exact_ec4_match,
        "exact_ec4_fractional": exact_ec4_fractional,
        "canonical_top_ec4": canon_top_ec4,
        "canonical_ec4_match": canonical_ec4_match,
        "canonical_ec4_fractional": canonical_ec4_fractional,
        "top_ec3": top_ec3,
        "top_ec3_tied": top_ec3_tied,
        "top_ec3_score": top_ec3_score,
        "ec3_match": ec3_match,
        "ec3_fractional": ec3_fractional,
    }


def _aggregate_truth_records(records: List[Dict[str, Any]]) -> Dict[str, Any]:
    n_total = len(records)
    n_with_pred = sum(1 for r in records if r["has_prediction"])
    exact_ec4_matches = sum(1 for r in records if r["exact_ec4_match"])
    exact_ec4_fractional = sum(float(r["exact_ec4_fractional"]) for r in records)
    canonical_ec4_matches = sum(1 for r in records if r["canonical_ec4_match"])
    canonical_ec4_fractional = sum(float(r["canonical_ec4_fractional"]) for r in records)
    ec3_matches = sum(1 for r in records if r["ec3_match"])
    ec3_fractional = sum(float(r["ec3_fractional"]) for r in records)

    return {
        "n_truth_proteins": n_total,
        "n_with_prediction": n_with_pred,
        "ec4_exact_matches": exact_ec4_matches,
        "ec4_exact_fractional": exact_ec4_fractional,
        "ec4_exact_rate": (exact_ec4_matches / n_total) if n_total else None,
        "ec4_canonical_matches": canonical_ec4_matches,
        "ec4_canonical_fractional": canonical_ec4_fractional,
        "ec4_canonical_rate": (canonical_ec4_matches / n_total) if n_total else None,
        "ec3_matches": ec3_matches,
        "ec3_fractional": ec3_fractional,
        "ec3_rate": (ec3_matches / n_total) if n_total else None,
    }


def evaluate_truth_agreement_arm(
    arm_summary: Dict[str, Any],
    manifest_targets_for_assembly: Optional[Dict[str, Dict[str, Any]]] = None,
    level: str = "retained",
    ec_entries: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Evaluate Price-148 experimental EC agreement for one arm at 'retained' or 'called' level.

    Proteins not in the Price-148 truth set are never scored as wrong (unconfirmed
    labels are never negatives).
    """
    if level not in ("called", "retained"):
        raise ValueError(f"Unsupported level: {level!r}")

    stored_truth = arm_summary.get("truth_proteins", {})
    if manifest_targets_for_assembly is not None:
        target_ids = sorted(manifest_targets_for_assembly.keys())
    else:
        target_ids = sorted(stored_truth.keys())

    scored_proteins: List[Dict[str, Any]] = []
    for pid in target_ids:
        prot_rec = stored_truth.get(pid)
        if prot_rec is None and manifest_targets_for_assembly is not None:
            meta = manifest_targets_for_assembly[pid]
            prot_rec = {
                "price_entry": meta["price_entry"],
                "proteome_protein_id": pid,
                "ground_truth_ecs": meta["ground_truth_ecs"],
                "ground_truth_ec3s": meta["ground_truth_ec3s"],
                "metabolic_stratum": meta["metabolic_stratum"],
                "retained": {},
                "called": {},
            }
        if prot_rec is not None:
            scored_proteins.append(
                _score_single_truth_protein(prot_rec, level=level, ec_entries=ec_entries)
            )

    by_stratum: Dict[str, Any] = {}
    for stratum in STRATA_ORDER:
        sub = [r for r in scored_proteins if r["metabolic_stratum"] == stratum]
        by_stratum[stratum] = _aggregate_truth_records(sub)

    return {
        "level": level,
        "unconfirmed_never_negative": True,
        "overall": _aggregate_truth_records(scored_proteins),
        "by_stratum": by_stratum,
        "proteins": {r["proteome_protein_id"]: r for r in scored_proteins},
    }


def evaluate_cell_summary(
    cell_summary: Dict[str, Any],
    taxon_constraints: Optional[BacterialTaxonConstraintIndex] = None,
    ec_entries: Optional[Dict[str, Any]] = None,
    price_targets_by_assembly: Optional[Dict[str, Dict[str, Dict[str, Any]]]] = None,
    exempt_assemblies: Optional[Iterable[str]] = None,
    exempt_taxids: Optional[Iterable[str]] = None,
) -> Dict[str, Any]:
    """Run all biological checks (a)-(d) on one cell's stage0 and stage1 summaries."""
    if ec_entries is None:
        ec_entries = load_ec_audit_entries()
    if taxon_constraints is None:
        taxon_constraints = load_bacterial_taxon_constraints(ec_entries=ec_entries)
    if price_targets_by_assembly is None:
        price_targets_by_assembly = load_price148_targets()

    assembly = str(cell_summary["assembly"])
    predictor = str(cell_summary["predictor"])
    arms = cell_summary.get("arms", {})

    stage0_taxon_id = arms.get("stage0", {}).get("config", {}).get("taxon_id")
    exemption = check_genome_exemption(
        assembly=assembly,
        taxon_id=stage0_taxon_id,
        exempt_assemblies=exempt_assemblies,
        exempt_taxids=exempt_taxids,
    )
    assembly_targets = price_targets_by_assembly.get(assembly, {})

    arm_results: Dict[str, Any] = {}
    for arm_name in ("stage0", "stage1"):
        arm_data = arms.get(arm_name, {})
        arm_results[arm_name] = {
            "translation": {
                "called": evaluate_translation_arm(
                    arm_data, level="called", ec_entries=ec_entries, exemption=exemption
                ),
                "retained": evaluate_translation_arm(
                    arm_data, level="retained", ec_entries=ec_entries, exemption=exemption
                ),
            },
            "core_metabolism": {
                "called": evaluate_core_metabolism_arm(
                    arm_data, level="called", ec_entries=ec_entries, exemption=exemption
                ),
                "retained": evaluate_core_metabolism_arm(
                    arm_data, level="retained", ec_entries=ec_entries, exemption=exemption
                ),
            },
            "taxon_consistency": evaluate_taxon_consistency_arm(
                arm_data, taxon_constraints=taxon_constraints
            ),
            "truth_agreement": {
                "retained": evaluate_truth_agreement_arm(
                    arm_data,
                    manifest_targets_for_assembly=assembly_targets,
                    level="retained",
                    ec_entries=ec_entries,
                ),
                "called": evaluate_truth_agreement_arm(
                    arm_data,
                    manifest_targets_for_assembly=assembly_targets,
                    level="called",
                    ec_entries=ec_entries,
                ),
            },
        }

    return {
        "schema": "gspa.biology_cell_evaluation/1",
        "cell_index": cell_summary.get("cell_index"),
        "predictor": predictor,
        "assembly": assembly,
        "exemption": exemption,
        "stage0_equals_stage1_sha256": bool(
            cell_summary.get("stage0_equals_stage1_sha256", False)
        ),
        "arms": arm_results,
    }


def evaluate_collection_summaries(
    cell_summaries: Iterable[Dict[str, Any]],
    taxon_constraints: Optional[BacterialTaxonConstraintIndex] = None,
    ec_entries: Optional[Dict[str, Any]] = None,
    price_targets_by_assembly: Optional[Dict[str, Dict[str, Dict[str, Any]]]] = None,
    exempt_assemblies: Optional[Iterable[str]] = None,
    exempt_taxids: Optional[Iterable[str]] = None,
) -> Dict[str, Any]:
    """Evaluate a collection of cell summaries and aggregate truth agreement by predictor and stratum."""
    if ec_entries is None:
        ec_entries = load_ec_audit_entries()
    if taxon_constraints is None:
        taxon_constraints = load_bacterial_taxon_constraints(ec_entries=ec_entries)
    if price_targets_by_assembly is None:
        price_targets_by_assembly = load_price148_targets()

    evaluated_cells: List[Dict[str, Any]] = []
    for summary in cell_summaries:
        evaluated_cells.append(
            evaluate_cell_summary(
                summary,
                taxon_constraints=taxon_constraints,
                ec_entries=ec_entries,
                price_targets_by_assembly=price_targets_by_assembly,
                exempt_assemblies=exempt_assemblies,
                exempt_taxids=exempt_taxids,
            )
        )

    evaluated_cells.sort(key=lambda c: (c["predictor"], c["assembly"]))

    # Aggregate Price-148 truth agreement per predictor x arm x level across all evaluated genomes
    predictors = sorted({c["predictor"] for c in evaluated_cells})
    truth_by_predictor: Dict[str, Any] = {}
    for pred in predictors:
        pred_cells = [c for c in evaluated_cells if c["predictor"] == pred]
        truth_by_predictor[pred] = {}
        for arm_name in ("stage0", "stage1"):
            truth_by_predictor[pred][arm_name] = {}
            for level in ("retained", "called"):
                all_prot_records: List[Dict[str, Any]] = []
                for cell in pred_cells:
                    prots = (
                        cell["arms"][arm_name]["truth_agreement"][level]["proteins"]
                    )
                    for pid in sorted(prots.keys()):
                        all_prot_records.append(prots[pid])
                by_stratum = {}
                for stratum in STRATA_ORDER:
                    sub = [
                        r
                        for r in all_prot_records
                        if r["metabolic_stratum"] == stratum
                    ]
                    by_stratum[stratum] = _aggregate_truth_records(sub)
                truth_by_predictor[pred][arm_name][level] = {
                    "n_genomes_evaluated": len(pred_cells),
                    "overall": _aggregate_truth_records(all_prot_records),
                    "by_stratum": by_stratum,
                }

    return {
        "schema": "gspa.biology_collection_evaluation/1",
        "citations": {
            "exemption_policy": EXEMPTION_POLICY_CITATION,
            "translation": TRANSLATION_RULE_CITATION,
            "core_metabolism": CORE_METABOLISM_RULE_CITATION,
            "taxon_consistency": TAXON_CONSISTENCY_CITATION,
            "truth_agreement": TRUTH_AGREEMENT_CITATION,
            "mutual_exclusivity_note": MUTUAL_EXCLUSIVITY_OMISSION_NOTE,
        },
        "rules": {
            "translation_charging_routes": TRANSLATION_CHARGING_RULES,
            "core_bacterial_ecs": CORE_BACTERIAL_EC_RULES,
            "exempt_genera": EXEMPT_ENDOSYMBIONT_GENERA,
            "exempt_taxids": EXEMPT_ENDOSYMBIONT_TAXIDS,
        },
        "n_cells_evaluated": len(evaluated_cells),
        "cells": evaluated_cells,
        "truth_agreement_by_predictor": truth_by_predictor,
    }


def load_summaries_from_path(summaries_path: Path) -> List[Dict[str, Any]]:
    path = Path(summaries_path)
    if path.is_file():
        data = json.loads(path.read_text(encoding="utf-8"))
        if data.get("schema") == "gspa.biology_cell_summary/1":
            return [data]
        if isinstance(data.get("cells"), list):
            return list(data["cells"])
        raise ValueError(f"Unrecognized summary file schema in {path}")

    cells_dir = path / "cells" if (path / "cells").is_dir() else path
    summaries: List[Dict[str, Any]] = []
    for json_file in sorted(cells_dir.glob("*.json")):
        if json_file.name.endswith(".receipt.json") or json_file.name in (
            "receipts.json",
            "evaluation.json",
        ):
            continue
        data = json.loads(json_file.read_text(encoding="utf-8"))
        if data.get("schema") == "gspa.biology_cell_summary/1":
            summaries.append(data)
    return summaries


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--summaries",
        type=Path,
        required=True,
        help="Directory of extracted cell summary JSON files (or a single cell summary JSON)",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Path to write the evaluated biology JSON report",
    )
    parser.add_argument(
        "--exempt-assembly",
        action="append",
        default=[],
        help="Optional assembly accession to mark exempt (may be repeated)",
    )
    parser.add_argument(
        "--exempt-taxid",
        action="append",
        default=[],
        help="Optional NCBI Taxon ID to mark exempt (may be repeated)",
    )
    args = parser.parse_args(argv)

    summaries = load_summaries_from_path(args.summaries)
    evaluation = evaluate_collection_summaries(
        summaries,
        exempt_assemblies=args.exempt_assembly,
        exempt_taxids=args.exempt_taxid,
    )
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps(evaluation, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
