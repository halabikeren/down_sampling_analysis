from ete3 import Tree
from samplers.PDA import PDA
from pydantic import ValidationError

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    tree = Tree("(Pseudomonas_aeruginosa:0.385658,(Vibrio_cholerae:1.50709,((Haemophilus_influenzae:0.424726,Pasteurella_multocida:0.414948):1.66758,(((Proteus_mirabilis:0.693306,Arsenophonus_nasoniae:1.13448):0.399439,((Xenorhabdus_bovienii:0.227453,Xenorhabdus_nematophila:0.261374):0.459266,(Photorhabdus_asymbiotica:2.71318,Photorhabdus_luminescens:0.336152):0.0357045):0.0171979):0.137946,((Serratia_proteamaculans:0.351046,(Yersinia_pestis:0.435987,(Hamiltonella_defensa:1.51349,Regiella_insecticola:0.773942):0.387599):0.159113):0.15838,((((Pectobacterium_carotovorum:0.542096,(Pectobacterium_atrosepticum:0.974646,Pectobacterium_wasabiae:0.459581):0.0418491):1e-06,(Dickeya_dadantii_Ech703:0.666053,(Dickeya_dadantii_Ech586:0.42516,Dickeya_zeae:0.404284):1e-06):0.0134038):0.00553918,((Pantoea_ananatis:0.475677,(Erwinia_billingiae:0.283726,(Erwinia_tasmaniensis:0.18256,(Erwinia_pyrifoliae:0.0616878,(Erwinia_amylovora_ATCC_49946:1e-06,Erwinia_amylovora_CFBP1430:1e-06):0.113845):0.0932479):0.192241):0.128946):0.186528,((Cronobacter_turicensis:0.296508,Cronobacter_sakazakii:0.694746):0.00687224,((Escherichia_coli:0.382068,(Citrobacter_rodentium:0.212508,(Salmonella_enterica:0.401972,Citrobacter_koseri:0.81171):0.272659):0.013973):0.00668039,((Klebsiella_variicola:0.258098,Klebsiella_pneumoniae:0.221267):0.0916799,(Enterobacter_sp_638:0.215251,Enterobacter_cloacae:0.207577):0.244522):0.00821813):0.119389):0.0325466):0.0461099):0.0281421,((Edwardsiella_tarda:0.0257235,Edwardsiella_ictaluri:0.15111):0.353013,(Sodalis_glossinidius:0.365534,(Baumannia_cicadellinicola:0.990899,(Ishikawaella_capsulata:1.86325,((((Buchnera_aphidicola_str_Sg:0.730076,Buchnera_aphidicola_str_APS:0.573):0.781011,(Buchnera_aphidicola_str_Cc:3.72234,Buchnera_aphidicola_str_Bp:1.94707):0.601071):0.509047,(Riesia_pediculicola:4.02404,Wigglesworthia_glossinidia:2.69957):0.7056):0.460529,(Blochmannia_pennsylvanicus:1.02226,Blochmannia_floridanus:1.73116):1.10015):0.533078):0.942166):1.08871):0.551261):0.0820946):1e-06):1e-06):1e-06):1e-06):0.00202107,Xanthomonas_axonopodis:0.433488);")
    taxon_to_weight = {'Pseudomonas_aeruginosa': 1,
                        'Vibrio_cholerae': 1,
                        'Haemophilus_influenzae': 1,
                        'Pasteurella_multocida': 1,
                        'Proteus_mirabilis': 1,
                        'Arsenophonus_nasoniae': 1,
                        'Xenorhabdus_bovienii': 1,
                        'Xenorhabdus_nematophila': 1,
                        'Photorhabdus_asymbiotica': 1,
                        'Photorhabdus_luminescens': 1,
                        'Serratia_proteamaculans': 1,
                        'Yersinia_pestis': 1,
                        'Hamiltonella_defensa': 1,
                        'Regiella_insecticola': 1,
                        'Pectobacterium_carotovorum': 1,
                        'Pectobacterium_atrosepticum': 1,
                        'Pectobacterium_wasabiae': 1,
                        'Dickeya_dadantii_Ech703': 1,
                        'Dickeya_dadantii_Ech586': 1,
                        'Dickeya_zeae': 1,
                        'Pantoea_ananatis': 1,
                        'Erwinia_billingiae': 1,
                        'Erwinia_tasmaniensis': 1,
                        'Erwinia_pyrifoliae': 1,
                        'Erwinia_amylovora_ATCC_49946': 1,
                        'Erwinia_amylovora_CFBP1430': 1,
                        'Cronobacter_turicensis': 1,
                        'Cronobacter_sakazakii': 1,
                        'Escherichia_coli': 1,
                        'Citrobacter_rodentium': 1,
                        'Salmonella_enterica': 1,
                        'Citrobacter_koseri': 1,
                        'Klebsiella_variicola': 1,
                        'Klebsiella_pneumoniae': 1,
                        'Enterobacter_sp_638': 1,
                        'Enterobacter_cloacae': 1,
                        'Edwardsiella_tarda': 1,
                        'Edwardsiella_ictaluri': 1,
                        'Sodalis_glossinidius': 1,
                        'Baumannia_cicadellinicola': 1,
                        'Ishikawaella_capsulata': 1,
                        'Buchnera_aphidicola_str_Sg': 1,
                        'Buchnera_aphidicola_str_APS': 1,
                        'Buchnera_aphidicola_str_Cc': 1,
                        'Buchnera_aphidicola_str_Bp': 1,
                        'Riesia_pediculicola': 1,
                        'Wigglesworthia_glossinidia': 1,
                        'Blochmannia_pennsylvanicus': 1,
                        'Blochmannia_floridanus': 1,
                        'Xanthomonas_axonopodis': 1}
    try:
        pda = PDA(tree=tree, taxon_to_weight=taxon_to_weight)
        pda.compute_sample(1)
    except ValidationError as e:
        print(e.json())


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
