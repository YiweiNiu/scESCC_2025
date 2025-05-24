
#####################################################

# genes

#####################################################

# dissociation-related genes
dissocitaion_genes.hs = c("MIR22HG", "MAFF", "ITPKC", "ATF3", "ZFP36L2", "SLC38A2", "HSP90AA1", "TPPP3", "NR4A1", "MYC", "SOCS3", "EIF5", "DDX5",
                          "DNAJB4", "IRF8", "MYD88", "ACTG1", "ARID5A", "BAG3", "UBC", "OSGIN1", "PER1", "NCKAP5L", "SRF", "CEBPB", "OXNAD1", "RIPK1",
                          "CEBPG", "TIPARP", "TUBB4B", "CEBPD", "LITAF", "TUBB6", "SBNO2", "IFRD1", "DNAJB1", "DES", "JUND", "CCNL1", "ZC3H12A",
                          "ZYX", "PPP1R15A", "NFKBIZ", "SKIL", "RHOH", "ERF", "IER2", "BRD2", "NPPC", "KLF6", "RAP1B", "MYADM", "ID3", "PNP", "RHOB",
                          "USP2", "DDX3X", "IER3", "EIF1", "SQSTM1", "LMNA", "PDE4B", "PNRC1", "SRSF5", "GADD45A", "FOS", "NCOA7", "PCF11", "HSP90AB1",
                          "SLC41A1", "KLF9", "EGR1", "IDI1", "MIDN", "HIPK3", "SERPINE1", "DNAJA1", "ODC1", "PHLDA1", "SAT1", "STAT3", "TRA2B",
                          "KLF4", "CSRNP1", "TPM3", "GEM", "WAC", "HSPE1", "TRA2A", "FOSB", "ERRFI1", "ATF4", "IER5", "PPP1CC", "MCL1", "KCNE4",
                          "HSPH1", "KLF2", "JUNB", "ZFP36L1", "TAGLN2", "ZFP36", "DUSP8", "SRSF7", "IL6", "BTG2", "BHLHE40", "ZFAND5", "HSPB1",
                          "BTG1", "MAFK", "TRIB1", "DCN", "NOP58", "MT2A", "MT1A", "ANKRD1", "SDC4", "H3F3B", "JUN", "DUSP1", "NFKBIA", "IRF1",
                          "GADD45G", "SLC10A6", "TNFAIP6", "HSPA1B", "HSPA1B", "RASSF1", "CXCL1", "GCC1", "HSPA5", "TNFAIP3", "HSPA8", "FOSL2",
                          "PXDC1", "EGR2", "CCN1", "ERFE", "NOCT")
dissocitaion_genes.mm = c("Mir22hg", "Maff", "Itpkc", "Atf3", "Zfp36l2", "Slc38a2", "Hsp90aa1", "Tppp3", "Nr4a1", "Myc", "Socs3", "Eif5", "Ddx5",
                          "Dnajb4", "Irf8", "Myd88", "Actg1", "Arid5a", "Bag3", "Ubc", "Osgin1", "Per1", "Nckap5l", "Srf", "Cebpb", "Oxnad1", "Ripk1",
                          "Cebpg", "Tiparp", "Tubb4b", "Cebpd", "Litaf", "Tubb6", "Sbno2", "Ifrd1", "Dnajb1", "Des", "Jund", "Ccnl1", "Zc3h12a",
                          "Zyx", "Ppp1r15a", "Nfkbiz", "Skil", "Rhoh", "Erf", "Ier2", "Brd2", "Nppc", "Klf6", "Rap1b", "Myadm", "Id3", "Pnp", "Rhob",
                          "Usp2", "Ddx3x", "Ier3", "Eif1", "Sqstm1", "Lmna", "Pde4b", "Pnrc1", "Srsf5", "Gadd45a", "Fos", "Ncoa7", "Pcf11", "Hsp90ab1",
                          "Slc41a1", "Klf9", "Egr1", "Idi1", "Midn", "Hipk3", "Serpine1", "Dnaja1", "Odc1", "Phlda1", "Sat1", "Stat3", "Tra2b",
                          "Klf4", "Csrnp1", "Tpm3", "Gem", "Wac", "Hspe1", "Tra2a", "Fosb", "Errfi1", "Atf4", "Ier5", "Ppp1cc", "Mcl1", "Kcne4",
                          "Hsph1", "Klf2", "Junb", "Zfp36l1", "Tagln2", "Zfp36", "Dusp8", "Srsf7", "Il6", "Btg2", "Bhlhe40", "Zfand5", "Hspb1",
                          "Btg1", "Mafk", "Trib1", "Dcn", "Nop58", "Mt2", "Mt1", "Ankrd1", "Sdc4", "H3f3b", "Jun", "Dusp1", "Nfkbia", "Irf1", "Gadd45g",
                          "Slc10a6", "Tnfaip6", "Hspa1b", "Hspa1a", "Rassf1", "Cxcl1", "Gcc1", "Hspa5", "Tnfaip3", "Hspa8", "Fosl2", "Pxdc1", "Egr2",
                          "Cyr61", "Fam132b", "Ccrn4l")

# response to heat (GO:0034605)
heat_genes.hs = c('DNAJC2', 'NF1', 'HSPA7', 'HSPA13', 'CETN1', 'NUP160', 'ATR', 'CAMK2B', 'CAMK2G',
                  'ARPP21', 'HSF1', 'BAG3', 'SUMO1', 'EIF2S1', 'NDC1', 'FKBP4', 'AAAS', 'EP300', 'NUP214',
                  'DNAJC7', 'ATP2A2', 'NUP93', 'MAPK1', 'CCAR2', 'RPTOR', 'HDAC2', 'STAC', 'IL1A', 'HSP90AB1',
                  'ST8SIA1', 'HSPA5', 'HIKESHI', 'NUP88', 'NUP210', 'NUP205', 'ZFAND1', 'HSBP1L1', 'RAE1', 'NUP133',
                  'HSPH1', 'TP53', 'CD34', 'NUP54', 'TRPV4', 'HSP90AA4P', 'HSP90AB2P', 'HSP90AB3P', 'PTGS2', 'CDKN1A',
                  'IRAK1', 'PTGES3', 'BAG1', 'HSPA8', 'NUP85', 'HSF5', 'MYOF', 'NUP42', 'NUP50', 'HSPA1L', 'ATM',
                  'CXCL10', 'CRYAB', 'CAMK2A', 'ANO1', 'NUP188', 'NUP58', 'TFEC', 'HMOX1', 'MLST8', 'CAMK2D',
                  'SEC13', 'RPA3', 'SLU7', 'NUP107', 'HSP90AA2P', 'HSFX1', 'BAG5', 'NUP62', 'MTOR', 'SIRT1',
                  'HSPA1B', 'HSPA14', 'MAPKAPK2', 'BAG2', 'PDCD6', 'SEH1L', 'NUP35', 'NUP37', 'NUP43', 'MAPK3',
                  'RPA2', 'DAXX', 'DNAJB6', 'HSPA6', 'BAG4', 'NUP153', 'ATXN3L', 'HSF4', 'HSPA2', 'YWHAE',
                  'DNAJB1', 'VCP', 'HSP90AB4P', 'POM121C', 'HSF2', 'HSPA9', 'CREBBP', 'TRPV1', 'POLR2D',
                  'PRKACA', 'FGF1', 'SCARA5', 'HSFX3', 'HSFY1', 'TCIM', 'HSPB8', 'RBBP7', 'POM121', 'NUP155',
                  'MAPT', 'CHORDC1', 'TPR', 'HSPA1A', 'HSFX4', 'DHX36', 'IER5', 'PDCL3', 'ATXN3', 'HTRA2',
                  'THBS1', 'LYN', 'HSP90AA1', 'NUP98', 'AKT1S1', 'GSK3B', 'STUB1', 'SLC52A3', 'RANBP2', 'HSBP1', 'CLPB', 'RPA1')

# mRNA catabolic processes (GO:0061014)
mRNA_catabolic.hs = c('PRR5L', 'CNOT8', 'RIDA', 'ZC3H12D', 'TNRC6A', 'CPEB3', 'YTHDF2', 'PUM1', 'ZC3H12A',
                      'POLR2G', 'GIGYF2', 'ROCK1', 'RBM24', 'TNRC6C', 'METTL14', 'MOV10', 'METTL3', 'METTL16',
                      'PLEKHN1', 'PABPC1', 'ZFP36L1', 'NANOS3', 'NANOS2', 'CNOT1', 'TOB1', 'HNRNPD', 'GTPBP1',
                      'ZC3HAV1', 'TRIM71', 'NANOS1', 'KHSRP', 'RC3H1', 'UPF1', 'PNPT1', 'FTO', 'CNOT6L', 'DHX36',
                      'CNOT7', 'YTHDF3', 'TNRC6B', 'MEX3D', 'ZFP36L2', 'ZFP36', 'BTG2', 'AGO2', 'HNRNPR', 'ROCK2')

# ribonuclease activity (GO:0004540)
ribonuclease_activity.hs = c('NOB1', 'EXOSC4', 'PNLDC1', 'ANG', 'RPP14', 'DCPS', 'CNOT8', 'ERN1', 'RIDA', 'SAMHD1',
                             'ERVK-10', 'ISG20', 'EXOSC1', 'SLFN14', 'ANGEL2', 'RNASE10', 'RNASE12', 'RNASE13', 'EDC3',
                             'DIS3', 'RNASE2', 'ZC3H12A', 'RNASE7', 'PAN3', 'TSEN2', 'RNASEH2A', 'CPSF3', 'ELAC1', 'DIS3L',
                             'ERVK-7', 'ERVK-8', 'HERVK_11', 'ERVK-9', 'DIS3L2', 'DICER1', 'CNOT2', 'ABCE1', 'EXOSC7',
                             'ISG20L2', 'MRPL44', 'AGO3', 'OAS2', 'NUDT16L1', 'RNASE3', 'ERVK-11', 'RNASE11', 'DCP2',
                             'DBR1', 'CWF19L1', 'EXOSC10', 'ERN2', 'TSEN34', 'EXD2', 'PIWIL2', 'ERI1', 'EXO1', 'USB1',
                             'ELAC2', 'SLFN13', 'POP5', 'TSEN54', 'SND1', 'PDE12', 'RNH1', 'RPP25', 'ENDOV', 'EXOSC5', 'EXOSC3', 'CNOT1', 'PPP1R8', 'ERVK-19', 'SMG6', 'RNASE9', 'TSNAX-DI', 'DGCR8', 'TSNAX', 'RCL1', 'EXOSC2', 'POP1', 'ERVK-25', 'EXOSC8', 'ERI2', 'ERVK-18', 'EXOSC6', 'RNASE1', 'RNASET2', 'ENDOU', 'PNPT1', 'RNASEL', 'HELZ2', 'PIWIL4', 'RNASEH2B', 'OASL', 'TSN', 'LACTB2', 'FEN1', 'DROSHA', 'ENDOG', 'PARN', 'RPP21', 'CNOT6L', 'RNASE4', 'NOCT', 'EXOG', 'RNASE6', 'HSPA1A', 'CNOT6', 'PIWIL1', 'RPP40', 'POP7', 'ERI3', 'AGO1', 'EXOSC9', 'TOE1', 'CNOT7', 'E7EUB6', 'REXO2', 'H0YAE9', 'RPP30', 'RPP38', 'ERVK-6', 'ANGEL1', 'OAS1', 'XRN2', 'APEX1', 'PRORP', 'RNASEH1', 'AZGP1', 'RNASE8', 'PAN2', 'YBEY', 'XRN1', 'OAS3', 'NPM1', 'TMBIM6', 'TCEA1', 'RNASEK', 'AGO2', 'NUDT16', 'POP4')

# cell division (GO:0051301)
cell_division.hs = c("PDCD6IP", "SMC4", "ACTR3", "ACTR2", "KIF20A", "CDK2", "CDK3", "CDK4", "CDK5", "CDK6", "CDK7", "CDK2AP2", "STAGMAEA", "KATNB1", "CCNO", "TUBA1B", "SYCP2", "ANAPC10", "RACK1", "NDC80", "MAD2L2", "TACC3", "CIB1", "ZNRD2", "CENPA", "SMC2", "CENPC", "USP16", "SPAG5", "STAMBP", "CENPE", "TXNIP", "CENPF", "RGS14", "CETN1", "CETNCETN3", "USP39", "CFL1", "NUDC", "STAG2", "GIPC1", "PLK2", "ARPP19", "NEK6", "SEPT9", "ENTR1", "RAB10", "JTB", "TXNL4A", "MAPRE2", "KIF2C", "RAB35", "RCC1", "CNTRL", "UBE2C", "TRIOBP", "KATNA1", "CIT", "ZWINCHEK2", "PMF1", "SNX18", "DCTN3", "CDCA5", "OIP5", "EFHC1", "CCDC124", "HAUS1", "CKS1B", "CKS2", "CNTROB", "MRGPRX2", "PARD3B", "ANAPC16", "CLTA", "CLTC", "NEDD1", "SEPT12", "PLK3", "MISP", "PLK5", "KDF1", "CCSAARL8A", "CHMP4B", "MITD1", "DIS3L2", "MPLKIP", "E2F7", "CSNK1A1", "KIF18B", "SPC24", "PHF13", "SEPT10", "PPP1R1C", "SGO2", "SGO1", "ROPN1B", "SPICE1", "CDCA2", "DCT", "SDE2", "DCTN1", "CDC14C", "SEPT1", "DYNC1H1", "DRD2DRD3", "E4F1", "ECT2", "KLHDC8B", "CENPV", "FLCN", "CENPX", "ENSA", "TUBB", "EPB41", "ZNF449", "EPB41L2", "SENP5", "ERCC2", "EREG", "ESRRB", "ETV5", "EVI2B", "STOX1", "CCNY", "SKA1", "REEP3", "SKA3", "FGF1FGF2", "FGF3", "FGF4", "FGF5", "FGF6", "FGF7", "FGF8", "FGF9", "FGF13", "FGFR1", "FGFR2", "VEGFD", "CEP164", "MAPRE1", "MAPRE3", "SIRT2", "TPX2", "PDS5B", "WAPL", "SPART", "CLASP2", "POGZ", "SMC5", "ANKLE2", "SEPT6", "ABRAXAS2", "SEPT8", "FBXL7", "PDS5A", "NCAPD3", "TTC28", "CLASP1", "HAUS5", "MAU2", "SPECC1L", "NCAPH", "ITGB3BP", "ZFYVE26", "CDK20", "ORC6", "CD2AP", "NUP62", "KIF4A", "CDC26", "HEPACAM2", "SYCE2", "SNX33", "ANAPC13", "ANAPC15", "AHCTF1", "ASPM", "NSL1", "WWTR1", "CHMP2B", "OR1A2", "FBXO5", "GNL3", "OPN1MW", "LATS2", "AATF", "CKAP2", "INTU", "VPS4A", "CHMP2A", "UBE2S", "CECR2", "GNAI1", "GNAI2", "GNAIGOLGA2", "SFN", "BOD1L2", "DLL1", "KIF4B", "NSMCE2", "ANK3", "GIT1", "SETD2", "NR3C1", "CHMP4A", "BABAM1", "RACGAP1", "ANAPC2", "GPSM2", "SAC3D1", "ANAPC4", "HELLS", "PIK3R4", "ANXA11", "HNRNPU", "HOXB4", "APC", "BIRCHTR2B", "DCDC1", "SEPT14", "MACC1", "STRA8", "IGF2", "SKA2", "NUP43", "IL1A", "IL1B", "INCENP", "ING2", "ARF1", "TOPAZ1", "CENPS", "KIF2A", "KIT", "ARF6", "KIF11", "KIFC1", "RHOA", "CENPW", "C10orf99", "INSC", "KMT5A", "RHOB", "RHOC", "STMN1", "LIG1", "LIG3", "LIG4", "LLGL2", "ARL3", "MIR145", "MAD2L1", "MAP4", "MDK", "MYC", "MYH9", "MYH10", "NAP1L2", "SEPT2", "NEDD9", "NEK1", "NEK2", "NEK3", "NKX3-1", "NOTCH1", "YBX1", "NUMA1", "OSM", "PAFAH1B1", "SYCP3", "CUZD1", "PAX6", "TAS2R13", "PARD6A", "SH3GLB1", "DYNC1LI1", "LEF1", "CHMP1A", "NUSAP1", "RXFP3", "FZR1", "SNX9", "ANAPC5", "ANAPC7", "SPOUT1", "CHMP5", "ANAPC11", "ZC3HC1", "PDGFA", "PDGFB", "CINP", "CHMP3", "CDK14", "PGF", "PIK3C3", "PIK3CB", "PIN1", "PLK1", "PELO", "MIS18A", "SEPT5", "SEPT4", "ANLN", "PIMREG", "POU3F2", "POU3F3", "POU5F1", "INO8MAP10", "PPBP", "ALKBH4", "HAUS6", "NDE1", "ERCC6L", "NSUN2", "NCAPG2", "TTC19", "SPDL1", "HAUS4", "TIPIN", "PPP1CA", "PPP1CB", "PPP1CC", "ZWILCH", "FIGN", "HAUS2", "CDCA8", "CEP55", "ARL8B", "MIS18BP1", "TRIM36", "HAUS7", "PRPF40A", "INTS13", "CHFR", "SEPT11", "PRKCE", "CENPJ", "PPP2R2D", "PKN2", "BIN3", "RCC2", "KLHL9", "SEPT3", "PDGFC", "TEX14", "SUSD2", "GKN1", "PARD3", "FMN2", "SPIRE1", "PDXP", "KNL1", "CHMP1B", "VANGL2", "PTCH1", "SPC25", "BIRC6", "KLHL42", "MICAL3", "PTN", "USP37", "MARK4", "BBS4", "RAD21", "RALA", "RALB", "RAN", "RASA1", "RB1", "RBBP8", "RBL1", "RBL2", "CCND1", "OPN1LW", "BCL2L1", "RFPLRPS3", "RTKN", "KIF13A", "PRDM15", "BLM", "NCAPG", "SFRP2", "TENT4B", "CXCR5", "SOX17", "SYCE3", "ANAPC1", "SHH", "GAREM1", "SON", "SOX5", "SPAST", "SPTBN1", "BRCA2", "SSTR5", "ZFP36L2", "NEK4", "AURKA", "AURKSVIL", "SYCP1", "BTC", "TACC1", "TAL1", "BUB1", "DYNLT3", "DYNLT1", "TEAD3", "BUB1B", "TERF1", "TGFA", "TGFB1", "TGFB2", "TGFB3", "THBS4", "TIAL1", "TOP1", "TOP2A", "TPR", "TSG101", "OPN1MW2", "UBE2I", "UVRAVEGFA", "VEGFB", "VEGFC", "VRK1", "WEE1", "WNT7A", "WNT9B", "ZNF16", "ZBTB16", "ZNF207", "EVI5", "TUBA1A", "MIS12", "NUP37", "BRCC3", "FSD1", "NOX5", "HAUS3", "OR2A4", "CHMP6", "E2F8", "CSPP1", "BORA", "MAP9MCMBP", "DSN1", "ANKRD53", "CALM1", "CENPT", "MYO19", "CEP63", "PDGFD", "REEP4", "CALM2", "CALM3", "TAS1R2", "HMGA2", "LBH", "FAM83D", "CDT1", "CABLES2", "SEH1L", "NCOA3", "USP9X", "SMC1A", "CDC7", "FZD7", "CDCANUF2", "MAD1L1", "SETDB2", "HMCN1", "BRIP1", "USP44", "PROK1", "DYRK3", "RAB11FIP4", "LZTS2", "SPIRE2", "CUL3", "PARD6G", "PARD6B", "KIF2B", "CAT", "PSRC1", "TUBA1C", "RAE1", "KLHL22", "MASTL", "ZFYVE19", "PKP4CCNB3", "DOCK7", "LRRCC1", "DIXDC1", "CDC14B", "CDC14A", "RUVBL1", "TP63", "NUMB", "TNKS", "BECN1", "CDC23", "RAB11A", "CCNK", "CDC123", "CDC16", "CCNA2", "CCNA1", "CCNB1", "TIMELESS", "CCND2", "CCND3", "WASLWNT3A", "LMLN", "CCNE1", "CCNF", "SAPCD2", "CCNG1", "CCNG2", "TMEM250", "BRSK2", "RNF8", "KLHL13", "CCNT1", "KNSTRN", "CCNT2", "PRC1", "UNC119", "USP8", "LATS1", "SMC3", "BOD1", "CCNB2", "CCNE2", "CTDP1", "ZNF830", "NEK9", "CABLES1", "CHMP7", "ARHGEF2", "ZW10", "BUB3", "AURKB", "PTTG1", "CHMP4C", "NUMBL", "ITGB1BP1", "HAUS8", "SYCE1", "KIF3B", "ACTR8", "RECQL5", "AC027237.1", "VPS4B", "BCAR1", "BABAM2", "KIF20B", "ESPL1", "RAB11FIP3", "KNTC1", "CCP110", "CKAP5", "IST1", "CUL7", "CDK1", "SEPT7", "WASHC5", "CDC6", "KLHL21", "CDC20", "NCAPD2", "KIF14", "CDC25A", "DCLRE1A", "CDC25B", "CDC25C", "CDC42")#
# keratinization (GO:0031424)
keratinization.hs = c("CDH3", "KRTAP25-1", "KRTAP4-9", "KRTAP21-3", "KRTAP9-7", "KRTAP16-1", "KRTAP9-6", "KRTAP29-1", "CDSN", "SPINK5", "PKP3", "KRT71", "KRT74", "KRT40", "RPTN", "KRTAP13-1", "KRT72", "KRT80", "KRT25", "DSG4", "CSTA", "KRT28", "SPRR4", "DSC1", "DSC2", "DSC3", "DSG1", "DSG2", "DSG3DSP", "KRT24", "KRT78", "LCE4A", "EVPL", "FLG", "KAZN", "CASP14", "LCE5A", "KRTAP15-1", "KLK5", "KRT23", "KLK13", "ABCA12", "LCE2B", "SFN", "KRTAP13-4", "KRT6C", "KRT73", "KRTAP8-1", "KRTAP11-1", "KRTAP19-1", "KRTAP13-2", "KRTAP13-3", "KRTAP23-1", "KRTAP6-1", "KRTAP6-2", "KRTAP6-3", "KRTAP19-2", "KRTAP19-3", "KRTAP19-4", "KRTAP19-5", "KRTAP19-6", "KRTAP19-7", "KRTAP20-1", "KRTAP20-2", "KRTAP21-1", "KRTAP21-2", "KRTAP22-1", "KRT79", "LIPM", "KRT27", "LCE1A", "LCE1B", "LCE1C", "LCE1D", "LCE1E", "LCE1F", "LCE2A", "LCE2LCE2D", "LCE3A", "LCE3B", "LCE3C", "LCE3E", "KRT26", "KRTAP12-2", "KRTAP12-1", "KRTAP10-10", "IVL", "JUP", "KRT77", "KRTAP5-9", "KRT1", "KRT2", "KRT3", "KRT4", "KRT5", "KRT6A", "KRT6KRT7", "KRT8", "KRT9", "KRT10", "KRT12", "KRT13", "KRT14", "KRT15", "KRTAP10-4", "KRTAP10-6", "KRTAP10-7", "KRTAP10-9", "KRTAP10-1", "KRTAP10-11", "KRTAP10-2", "KRTAP10-5   KRTAP10-8", "KRTAP10-3", "KRTAP12-3", "KRTAP12-4", "KRTAP10-12", "KRT16", "KRT17", "KRTAP5-1", "KRTAP5-3", "KRTAP5-4", "KRTAP5-10", "KRT18", "KRT19", "KRT31", "KRT3KRT33A", "KRT33B", "KRT34", "KRT35", "HRNR", "KRT81", "KRT82", "KRTAP26-1", "KRT83", "KRT84", "KRT85", "KRT86", "KRT39", "LOR", "SPINK6", "KLK14", "KLK12", "KRTAP5-5", "KRTAP5-2", "KRTAP5-6", "KRTAP5-7", "KRTAP5-11", "LCE6A", "FURIN", "PCSK6", "KRT76", "PPHLN1", "PI3", "PKP1", "PKP2", "KRT20", "PPL", "PRSS8", "CYP26B1", "KRTAP5-8", "CELA2A", "PERP", "SPINK9", "LIPK", "LIPN", "KRTAP24-1", "KRTAP27-1", "KRTAP4-11", "SPRR1A", "SPRR1B", "SPRR2A", "SPRR2B", "SPRR2D", "SPRR2E", "SPRR2F", "SPRR2G", "SPRR3", "ST14", "TGM1", "TGM3", "TCHH", "KRTAP4-8", "KRTAP1-4", "KRTAP19-8", "KRTAP9-1", "KRTAP2-3", "KRTAP1-3", "KRTAP1-1", "SHARPIN", "KRTAP9-9", "KRTAP4-6", "KRTAP2-1", "CAPN1", "CAPNS1", "KRTAP1-5", "KRTAP3-1", "KRTAP3-2", "KRTAP9-2", "KRTAP9-3", "KRTAP9-8", "KRTAP17-1", "TMEM79", "CNFN", "KRTAP4-4", "LCE3D", "PKP4", "KRTAP9-4", "KRTAP4-1", "KRTAP4-5", "KRTAP4-3", "KRTAP4-2", "KRTAP3-3", "KRTAP2-4", "KRT38", "KRT37", "KRT36", "KRT75", "TGM5")
# EPITHELIAL_CELL_PROLIFERATION (GO:0050673)



# canonical markers
markers.immune = c("PTPRC")
markers.tcells = c(
  markers.immune,
  # T
  "CD3D", "CD3G", "CD3E", "CD247",
  # CD4, CD8
  "CD4", "CD8A", "CD8B",
  # Treg
  "FOXP3", "IL2RA",
  # MAIT
  "TRAV1-2", "KLRB1",
  # gamma-delta T
  "TRDV1", "TRDV2", "TRGV9", "TRGV2",
  # Tfh
  'CXCR5', 'PDCD1', 'ICOS', 'BCL6')

markers.NK = c(markers.immune,
               'NCAM1', # CD56
               'FCGR3A', # CD16
               "KLRF1", "NKG7", "NCR1", "XCL2")

markers.ILC = c('IL7R', 'AREG', 'TNFSF13B')

markers.bcells = c(
  markers.immune,
  # B (CD19, CD20)
  "CD19", "MS4A1",
  # B1/transitional (CD24 high, CD38 high)
  'CD5', 'CD24', 'CD38', 'EBF1', 'TCF3',
  # naive
  "IGHM", "IGHD",
  # Breg
  "IL10", 'CD1D',
  # follicular B cells
  # CD23, CD21, CD22
  "FCER2", "CR2", "CD22",
  # GC
  'MME', 'TNFRSF13C', 'TNFRSF17', 'BCL6',
  # memory
  "CD27", 'CD40', 'SPIB', 'PAX5', 'POU2AF1'
)
markers.plasma = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "JCHAIN")

markers.monocyte = c(
  # CD14, CD16
  "CD14", "FCGR3A",
  'VCAN', 'FCN1'
)

markers.macrophage = c('CD163', "C1QA", 'CD68', "APOE", "AIF1",
                       # Scavenger receptor
                       "MARCO", 'MSR1',
                       # Mineralocorticoid Receptor
                       "NR3C2", 'NR3C1')
markers.mps = unique(c(markers.immune, markers.monocyte, markers.macrophage))

markers.DC = c(markers.immune,
               # plasmacytoid DCs
               "IL3RA", # CD123
               "CLEC4C",
               # cDC1
               'CLEC9A', 'ITGAX', 'IRF8', 'BATF3',
               # cDC2
               "CD1C", "CLEC10A", 'IRF4',
               # tDC
               "LAMP3", "IDO1", "FSCN1", "CD200",
               # activated DCs
               'CCR7', 'CD83', 'IL7R', 'ID2',
               #Langerhans cells
               'CD207', 'CD1A', 'ZBTB46'
               )

markers.mast = c(markers.immune,
                 # tryptase
                 "TPSAB1", "TPSB2", 'TPSD1', 'TPSG1',
                 # chymase
                 'CMA1',
                 # carboxypeptidases
                 "CPA3", 'CPN1', 'CPB2', 'CPE',
                 # Leukotriene
                 'LTC4S', 'LTA4H',
                 # Prostaglandin
                 'PTGS2', 'PTGS1',
                 # activation
                 'TNF','IL4', 'IL13',
                 # other
                 "MS4A2", 'KIT')

marker.neutrophil = c(markers.immune,
                      "ELANE", "S100A8", "S100A9", "MMP9")

markers.myeloid = unique(c(markers.immune,
                           "LYZ", "HLA-DRA",
                           # mono/macrophage
                           markers.mps,
                           # DC
                           markers.DC,
                           # mast
                           markers.mast))

markers.platelets = c("PPBP", "PF4", "ITGA2B", "ITGB3")
markers.erythrocyte = c("HBB", "HBA2", "HBA1")
markers.epithelia = c("KRT6A", "KRT5", "KRT13", "KRT4", "ECM1", "SPRR3", "CNFN", "CRNN")
markers.fibroblasts = c("COL1A1", "COL1A2", "COL3A1", "LUM", "DCN", "C1R")
markers.endothelia = c("PECAM1", # also expressed in monocytes
                       "CLDN5", "CDH5", "VWF", "FLT1", "FLT4")
markers.pericyte = c('RGS5', 'PDGFRB', 'MCAM',
                     # NG2
                     'CSPG4')

markers.mesothelia = c('LRRN4', 'UPK3B', 'CDH1', 'MSLN')


# signatures used to fast classify
markers.to.plot = list(
  immune = c("PTPRC"),
  tcells = c("CD3D", "CD3G", "CD3E", "CD247", "CD4", "CD8A", "CD8B"),
  NK = c("NCAM1", "NKG7", "KLRF1", "NCR1", "XCL2"),
  bcells = c("CD19", "MS4A1", "IGHM", "IGHD"),
  plasma = c("IGHG1", "IGHG4", "IGHA1", "IGHA2"),
  MP = c("CD14", "FCGR3A", 'HLA-DRA', "LYZ", 'VCAN', 'FCN1', "APOE", "AIF1", "MARCO"),
  DC = c("IL3RA", "CLEC4C", "LILRA4"),
  mast = c("TPSAB1", "TPSB2", "CPA3", "MS4A2"),
  platelets = markers.platelets,
  erythrocytes = markers.erythrocyte,
  epithelia = markers.epithelia,
  fibroblasts = markers.fibroblasts,
  endothelia = markers.endothelia,
  pericyte = markers.pericyte,
  mesothelia = markers.mesothelia
)

#####################################################
# MHC
# from: https://www.genenames.org/data/genegroup/#!/group/588
#####################################################

MHC_genes = list(
  MHC1 = c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HLA-H", "HLA-J", "HLA-K", "HLA-L", "HLA-N", "HLA-P", "HLA-S", "HLA-T", "HLA-U", "HLA-V", "HLA-W", "HLA-X", "HLA-Y", "HLA-Z"),
  MHC2 = c("HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPA2", "HLA-DPA3", "HLA-DPB1", "HLA-DPB2", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DQB3", "HLA-DRA", "HLA-DRB1", "HLA-DRB2", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DRB6", "HLA-DRB7", "HLA-DRB8", "HLA-DRB9")
)

#####################################################
# Immune checkpoint ligands of antigen-presenting cell
# from: SnapShot: APC/T-Cell Immune Checkpoints
#####################################################

checkpoint_genes = list(
  ligand = c("PDCD1LG2", "CD274", "CD80", "CD86", "ICOSLG", "TNFRSF14", "CD40", "HHLA2", "PVR", "NECTIN2", "TIMD4", "CEACAM1", "LGALS9", "SLAMF1", "LY9", "CD48", "CD84", "SLAMF6", "SLAMF7", "CD70", "TNFSF4", "TNFSF8", "TNFSF9", "TNFSF14", "TNFSF18", "TNFSF15", "CALM1", "CLEC7A", "CLEC10A", "CD276", "C1QA", "C1QB", "C1QC", "VSIR", "BTN2A2", "BTN3A1", "BTNL2", "SIGLEC15", "VTCN1"),
  receptor = c("PDCD1", "CD80", "CD28", "CTLA4", "ICOS", "CD274", "BTLA", "TNFSF14", "CD40LG", "TMIGD2", "TIGIT", "CD226", "PVRIG", "HAVCR1", "HAVCR2", "SLAMF1", "LY9", "CD244", "CD84", "SLAMF6", "SLAMF7", "CD27", "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFRSF14", "TNFRSF18", "TNFRSF25", "ITGAL", "LGALS9", "LILRB1", "PTPRC", "TREML2", "LAG3", "LAIR1", "SELPLG")
)

#####################################################

# marker used to calculate score

#####################################################

epithelia_score_genes = c('EPCAM', 'ECM1', 'SPRR3', 'CNFN', 'CRNN', 'KRT1', 'KRT10', 'KRT13',
                          'KRT14', 'KRT15', 'KRT16', 'KRT17', 'KRT18', 'KRT19', 'KRT2', 'KRT222',
                          'KRT23', 'KRT24', 'KRT31', 'KRT4', 'KRT5', 'KRT6A', 'KRT6B', 'KRT6C',
                          'KRT7', 'KRT72', 'KRT73', 'KRT75', 'KRT78', 'KRT8', 'KRT80', 'KRT81',
                          'KRT86', 'KRTAP16-1', 'KRTCAP2', 'KRTCAP3', 'KRTDAP')

prolif_score_genes = c(
  # S
  "MCM5", "PCNA", "TYMS", "FEN1", "MCM7", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU", "HELLS", "RFC2", "POLR1B", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "MRPL36", "E2F8",
  # G2M
  "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "PIMREG", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "JPT1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"
)

#####################################################

# soupx genes

#####################################################

soupx_genes = list(
  ig = c('IGHA1', 'IGHA2', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'IGHD', 'IGHE',
         'IGHM', 'IGLC1', 'IGLC2', 'IGLC3', 'IGLC4', 'IGLC5', 'IGLC6', 'IGLC7', 'IGKC'),
  tcells = c("CD3D", "CD3G"),
  bcells = c('CD19', 'MS4A1'),
  mast = c("TPSAB1", "TPSB2"),
  erythrocyte = markers.erythrocyte,
  platelet = c("PPBP", "PF4"),
  epithelia = c('KRT6A', 'KRT13'),
  fibroblasts = markers.fibroblasts,
  endothelia = c("CLDN5", "FLT1", "CDH5"),
  mp = c('LYZ', 'CD14', 'AIF1')
)

#####################################################
# format transformation
#####################################################

seurat_2_cds <- function(srat_obj=NULL){
  # part one, gene annotations
  gene_annotation <- as.data.frame(rownames(srat_obj),
                                   row.names = rownames(srat_obj))
  colnames(gene_annotation) <- "gene_short_name"

  # part two, cell information
  cell_metadata <- srat_obj@meta.data

  # part three, counts sparse matrix
  New_matrix <- srat_obj@assays[["RNA"]]@counts
  expression_matrix <- New_matrix

  ### Construct the basic cds object
  cds <- new_cell_data_set(expression_matrix,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)
  ### Could be a space-holder, but essentially fills out louvain parameters
  cds@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

  ### Assign UMAP coordinate
  if ("umap" %in% names(srat_obj@reductions)) {
    cds@int_colData@listData[["reducedDims"]][["UMAP"]] <- srat_obj@reductions[["umap"]]@cell.embeddings
  } else {
    cds@int_colData@listData[["reducedDims"]][["UMAP"]] <- srat_obj@reductions[["UMAP"]]@cell.embeddings
  }

  ### Assign feature loading for downstream module analysis
  if ("pca" %in% names(srat_obj@reductions)) {
    cds@preprocess_aux$gene_loadings <- srat_obj@reductions[["pca"]]@feature.loadings
  } else {
    cds@preprocess_aux$gene_loadings <- srat_obj@reductions[["PCA"]]@feature.loadings
  }
  return(cds)
}



#####################################################

# filter_features

#####################################################

# filter out features that not in seurat_obj
filter_features <- function(seurat_obj=NULL, features=NULL) {
  tmp = unique(features)
  filt_f <- features[features %in% rownames(seurat_obj)]
  return(filt_f)
}


# percent above
# test = GetAssayData(seurat, assay = 'RNA', slot = 'data')[tmp,]
# percent_above(test)
percent_above <- function(data=NULL, threshold=0) {
  # data: gene x cell matrix
  res = apply(data, 1, FUN = Seurat:::PercentAbove, threshold=threshold)
  tibble(gene = names(res), res = res)
}



#####################################################

# SoupX

#####################################################

# read 10x counts, return a dataframe for soupx, which includes UMAP data and cluster info
.prepare_4_soupx <- function(count_dir = NULL) {
  seurat.data = Read10X(data.dir = count_dir)
  seurat <- CreateSeuratObject(counts = seurat.data)
  seurat = fast_cluster(seurat_obj = seurat, res = 1)
  data_2_fetch <- FetchData(seurat,
                            vars = c("seurat_clusters", "UMAP_1", "UMAP_2")) %>%
    dplyr::select(RD1 = UMAP_1, RD2 = UMAP_2, Cluster = seurat_clusters)
  return(data_2_fetch)
}

prepare_4_soupx <- function(data_dir = NULL, count_dir= NULL) {
  # extra data for soupx
  DR = .prepare_4_soupx(count_dir = count_dir)

  # construct SoupChannel
  sc = load10X(data_dir, keepDroplets = TRUE)
  # profiling the soup
  sc = estimateSoup(sc)
  # add dimenstion-reduction
  sc = setDR(sc, DR)
  # add cluster info
  sc = setClusters(sc, DR$Cluster)

  return (list(sc=sc, DR=DR))
}

#####################################################

# cell QC

#####################################################


# filter low-quality cells
lowQuality_filter <- function(seurat_obj, nmad = 3, gene_fixed_low=NULL, gene_fixed_high=NULL){
  # filter low quality
  reasons = scater::quickPerCellQC(seurat_obj@meta.data, lib_size = "nCount_RNA",
                                   n_features = "nFeature_RNA",
                                   percent_subsets=c("percent.mt", "percent.ribo",
                                                     "percent.dissociation", "percent.heat"),
                                   batch = seurat_obj$Source, nmad = nmad) %>%
    as.data.frame() %>%
    mutate(Source = seurat_obj$Source) %>%
    mutate(high_n_features = FALSE) %>% # add one column for gene_fixed_high filter
    dplyr::select(Source,
                  low_lib_size,
                  low_n_features, high_n_features,
                  high_percent.mt, high_percent.ribo,
                  high_percent.dissociation, high_percent.heat) # change column order

  # further filter by fixed threshold
  if(!is.null(gene_fixed_low)) {
    reasons$low_n_features = reasons$low_n_features | (seurat_obj$nFeature_RNA < gene_fixed_low)
  }
  if(!is.null(gene_fixed_high)) {
    reasons$high_n_features = seurat_obj$nFeature_RNA > gene_fixed_high
  }
  # discard
  reasons$discard = reasons$low_lib_size |
    reasons$low_n_features | reasons$high_n_features |
    reasons$high_percent.mt | reasons$high_percent.ribo |
    reasons$high_percent.dissociation | reasons$high_percent.heat

  # no filter, return reasons
  return(reasons)
}

# compute doublet score by scran::doubletCells
doublet_score <- function(i, seurat_obj) {
  # fetch one sample
  s = FetchData(object = seurat_obj, vars = 'Source')
  seuratObj <- seurat_obj[, which(x = s$Source == i)]
  ## fast_cluster using high resolution
  seuratObj = fast_cluster(seuratObj, res = 3)
  # doublets detection
  sce = as.SingleCellExperiment(seuratObj)
  #--- variance-modelling ---#
  set.seed(00010101)

  dec <- modelGeneVarByPoisson(sce)
  top <- getTopHVGs(dec, n = 1000)

  set.seed(100)
  dbl.dens <- doubletCells(sce, subset.row = top,
                           d = ncol(reducedDim(sce)))
  # summary(dbl.dens)
  sce$doubletScore <- log10(dbl.dens + 1)
  seuratObj = as.Seurat(sce)

  # doublet
  ds = seuratObj@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(median_doubletScore = median(doubletScore),
              num_cells = n())
  # median
  mu = median(ds$median_doubletScore)
  # mad, calculated only on values above the median to avoid the effects of zero-truncation
  mad = mad(ds$median_doubletScore[ds$median_doubletScore > mu])
  # This function is very similar to transform but it executes the transformations iteratively so that later transformations can use the columns created by earlier transformations. Like transform, unnamed components are silently dropped.
  a = ds %>%
    mutate(p = pnorm(median_doubletScore, mean = mu, sd = mad, lower.tail = F)) %>% # use pnorm to compute the 1-cdf
    mutate(p.adj = p.adjust(p, method = 'BH')) %>%
    filter(p.adj < 0.1)

  seuratObj$high_doubletScore = ifelse(seuratObj$seurat_clusters %in% a$seurat_clusters, TRUE, FALSE)
  return(seuratObj)
}

# compute doublet score by scds::cxds_bcds_hybrid
.doublet_scds <- function(seurat_obj) {
  ## fast_cluster using high resolution
  seuratObj = fast_cluster(seurat_obj, res = 3)
  # doublets detection
  sce = as.SingleCellExperiment(seuratObj)
  sce = cxds_bcds_hybrid(sce)
  seuratObj = as.Seurat(sce)

  # doublet
  ds = seuratObj@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(median_doubletScore = median(hybrid_score),
              num_cells = n())
  # median
  mu = median(ds$median_doubletScore)
  # mad, calculated only on values above the median to avoid the effects of zero-truncation
  mad = mad(ds$median_doubletScore[ds$median_doubletScore > mu])
  # This function is very similar to transform but it executes the transformations iteratively so that later transformations can use the columns created by earlier transformations. Like transform, unnamed components are silently dropped.
  a = ds %>%
    mutate(p = pnorm(median_doubletScore, mean = mu, sd = mad, lower.tail = F)) %>% # use pnorm to compute the 1-cdf
    mutate(p.adj = p.adjust(p, method = 'BH')) %>%
    filter(p.adj < 0.1)

  seuratObj$high_scds = ifelse(seuratObj$seurat_clusters %in% a$seurat_clusters, TRUE, FALSE)
  return(seuratObj)
}

# batch function for assembled-seurat obj
doublet_scds <- function(i, seurat_obj) {
  # fetch one sample
  s = FetchData(object = seurat_obj, vars = 'Source')
  seuratObj <- seurat_obj[, which(x = s$Source == i)]
  seuratObj <- .doublet_scds(seuratObj)
  return(seuratObj)
}

# doublet rate of 10x
doublet_rate <- function(cell_num) {
  x = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
  y = c(0.004, 0.008, 0.016, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076)
  fit = lm(y ~ x)
  predict(fit, data.frame(x = cell_num))
}


#####################################################

# PCA

#####################################################

# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# determine the number of PCs
#1. The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
#2. The point where the percent change in variation between the consecutive PCs is less than 0.1%.
#3. choose the smaller one
determine_pc_num <- function(seurat_obj) {
  # Determine percent of variation associated with each PC
  pct <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]

  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  return(c(co1, co2))
}

# get correlation of PC scores and confounders like percent.mt, percent.ribo etc.
confounder_check <- function(seurat_obj, features) {
  my_cor <- function(x, seurat_obj, features) {
    cor(x, seurat_obj@meta.data[, c(features)])
  }
  my_pheatmap <- function(mtx, zero_color=0) {
    paletteLength <- 100
    myColor <- colorRampPalette(c("blue","white","red"))(paletteLength)
    myBreaks <- c(seq(min(mtx), 0, length.out=ceiling(paletteLength/2) + 1),
                  seq(max(mtx)/paletteLength, max(mtx), length.out=floor(paletteLength/2)))
    pheatmap(mtx,
             cluster_rows = F, cluster_cols = F,
             show_rownames = F,
             color = myColor, breaks = myBreaks,
             fontsize_row = 8, fontsize_col = 12)
  }
  # get cell embeddings
  cell_embeddings = Embeddings(object = seurat_obj, reduction = "pca")
  # get correlations between PCs and confounders
  a = apply(cell_embeddings, 2, my_cor, seurat_obj = seurat_obj, features = features)
  rownames(a) = features
  my_pheatmap(t(a))
}


#####################################################

# cluster

#####################################################

# fast clustering using default parameters
fast_cluster <- function(seurat_obj, res = 0.8) {
  .seurat_obj = NormalizeData(seurat_obj, verbose = F)
  .seurat_obj = FindVariableFeatures(.seurat_obj, nfeatures = 2000, verbose = F)
  .seurat_obj = ScaleData(.seurat_obj, features = VariableFeatures(.seurat_obj), verbose = F)
  # remove existing
  if ('pca' %in% names(.seurat_obj@reductions)) {.seurat_obj[['pca']] <- NULL}
  if ('PCA' %in% names(.seurat_obj@reductions)) {.seurat_obj[['PCA']] <- NULL}
  .seurat_obj = RunPCA(.seurat_obj, features =  VariableFeatures(.seurat_obj), npcs = 50, verbose = F)
  pc_num = min(determine_pc_num(.seurat_obj))
  .seurat_obj <- FindNeighbors(.seurat_obj, dims = 1:pc_num, verbose = F)
  .seurat_obj <- FindClusters(.seurat_obj, res = res, verbose = F)
  # remove existing
  if ('umap' %in% names(.seurat_obj@reductions)) {.seurat_obj[['umap']] <- NULL}
  if ('UMAP' %in% names(.seurat_obj@reductions)) {.seurat_obj[['UMAP']] <- NULL}
  .seurat_obj <- RunUMAP(.seurat_obj, dims = 1:pc_num, verbose = F)
  return(.seurat_obj)
}

# sub-cluster cells
# criterion: n-1 groups with >=10 cells and >=10 differential expressed genes (lfc >= 1)
sub_cluster <- function(seurat_obj, res_init = 0.1, res_end = 2, res_step = 0.1, threshold_cells = 10, threshold_lfc = 1) {
  res = res_init
  while (res <= res_end) {
    print(paste0('resoulution: ', res))
    # find clusters
    seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = F)
    # cluster number
    cluster_num = length(unique(seurat_obj$seurat_clusters))
    # cluster assignment
    x = seurat_obj$seurat_clusters
    x = as.numeric(names(table(x)[table(x) > 10]))
    c1 = 0 # number of clusters under criterion 1
    for (i in x) {
      markers <- FindMarkers(seurat_obj, ident.1 = i, only.pos = TRUE, min.pct = 0.25, verbose = F)
      #print(i)
      #print(dim(markers))
      if (dim(markers[markers$avg_logFC >= threshold_lfc,])[1] >= 10) {
        c1 = c1 + 1
      }
    }
    if (c1 < cluster_num - 1){
      print(res)
      break
    } else {
      res = res + res_step
    }
  }
}

# clustering (sub-clustering)
get_cluster <- function(seurat_obj) {
  # get variable features
  get_vars <- function(seurat_obj) {
    vst_each = c()
    source = FetchData(seurat_obj, vars = c('Source'))
    for (i in unique(seurat_obj$Source)) {
      x = seurat_obj[,which(source == i)]
      if (dim(x)[2] < 100) {
        next
      } else {
        x = FindVariableFeatures(x, selection.method = 'vst', nfeatures = 500)
        vst_each = c(VariableFeatures(x), vst_each)
      }
    }
    VariableFeatures(seurat_obj) = unique(vst_each)
    return(seurat_obj)
  }
  # get pc
  x = get_vars(seurat_obj)
  vars.to.regress = c('nCount_RNA', 'percent.mt', 'percent.heat', 'percent.ribo', 'CC.Difference')
  x <- ScaleData(x, features = VariableFeatures(x), vars.to.regress = vars.to.regress, verbose = F)
  x = RunPCA(x, features = VariableFeatures(x), verbose = F)
  n_pc = min(determine_pc_num(x))
  # non-linear dimensional reduction (UMAP/tSNE)
  x = RunUMAP(x, dims = 1:n_pc, verbose = F)
  x = RunTSNE(x, dims = 1:n_pc, tsne.method = "FIt-SNE", nthreads = 4, max_iter = 2000, verbose = F)
  # clustering
  x = FindNeighbors(x, dims = 1:n_pc, verbose = F)
  resolutions = seq(0, 2, 0.2)
  x <- FindClusters(x, resolution = resolutions, verbose = F)
  return(list(n_pc = n_pc, resolutions = resolutions, seurat_obj = x))
}


#####################################################

# umap/tsne/diffusion-map

#####################################################

# run diffusion map on seurat obj
dm_on_seurat <- function(seurat_obj=NULL) {
  seuratObj = seurat_obj
  data = GetAssayData(seuratObj, assay = 'RNA', slot = "scale.data")
  n_pcs = min(determine_pc_num(seuratObj))
  dm <- destiny::DiffusionMap(t(data), n_pcs = n_pcs, n_eigs = 5)
  seuratObj@reductions$dm = CreateDimReducObject(embeddings = dm@eigenvectors, key = "DC_", assay = 'RNA')
  seuratObj@misc$dm = dm
  return(seuratObj)
}

get_projection <- function(seurat_obj=NULL, rm_mt_rb=FALSE, run_df = FALSE) {
  if (rm_mt_rb==TRUE) {
    tmp = rownames(seurat_obj)
    mt_rb = grep('^MT-|^RPL|^RPS', tmp, value = T)
    tmp = setdiff(tmp, mt_rb)
    x = seurat_obj[tmp, ]
  } else {
    x = seurat_obj
  }
  # get variable features
  get_vars <- function(seurat_obj) {
    vst_each = c()
    source = FetchData(seurat_obj, vars = c('Source'))
    for (i in unique(seurat_obj$Source)) {
      x = seurat_obj[,which(source == i)]
      if (dim(x)[2] < 100) {
        next
      } else {
        x = FindVariableFeatures(x, selection.method = 'vst', nfeatures = 500)
        vst_each = c(VariableFeatures(x), vst_each)
      }
    }
    VariableFeatures(seurat_obj) = unique(vst_each)
    return(seurat_obj)
  }
  # get pc
  x = get_vars(seurat_obj)
  message('Scale data')
  vars.to.regress = c('nCount_RNA', 'percent.mt', 'percent.heat', 'percent.ribo', 'CC.Difference')
  x <- ScaleData(x, features = VariableFeatures(x), vars.to.regress = vars.to.regress, verbose = F)
  message('Run PCA')
  x = RunPCA(x, features = VariableFeatures(x), verbose = F)
  n_pc = min(determine_pc_num(x))
  # non-linear dimensional reduction (UMAP/tSNE)
  message('Run UMAP')
  x = RunUMAP(x, dims = 1:n_pc, verbose = F)
  message('Run TSNE')
  x = RunTSNE(x, dims = 1:n_pc, verbose = F)
  if (run_df) {
    # diffusion-map
    message('diffusion-map')
    data = GetAssayData(x, assay = 'RNA', slot = "scale.data")
    dm <- destiny::DiffusionMap(t(data), n_pcs = n_pc, n_eigs = 5)
    x@reductions$dm = CreateDimReducObject(embeddings = dm@eigenvectors, key = "DC_", assay = 'RNA')
    x@misc$dm = dm
  }
  return(x)
}

# cell freq of Source
get_cell_freq <- function(srat = NULL, cell_prefix = NULL, groupby = "seurat_clusters") {
  require(dplyr)
  sam_info = srat@misc$sam_info %>%
    mutate(cell.filtered = as.numeric(cell.filtered))

  # get freq
  freq = srat@meta.data %>%
    group_by(!! sym(groupby), Source) %>%
    summarise(n = n()) %>%
    mutate(n = as.numeric(n)) %>%
    # change name of groupby
    dplyr::rename(group_v = !! sym(groupby))

  # add samples that do not exist
  tmp = lapply(unique(freq$group_v), FUN = function(x){
    tmp_sources = freq %>% filter(group_v == x) %>% pull(Source) %>% as.character()
    diff_source = as.character(setdiff(sam_info$Source, tmp_sources))
    #print(class(diff_source))
    if (length(diff_source) == 0) {
      return(data.frame())
    }
    return(data.frame(group_v = x, Source = diff_source, n = 0))
  })
  tmp = do.call(rbind, tmp)

  # bind
  freq = rbind(freq, tmp)
  # combine sample info
  freq = freq %>%
    left_join(sam_info, by = "Source") %>%
    mutate(prop = n/cell.filtered)
  # add prefix to groupby
  if (!is.null(cell_prefix)) {
    freq = freq %>%
      mutate(group_v = paste(cell_prefix, group_v, sep = "-"))
  }
  # wider
  freq = freq %>%
    pivot_wider(id_cols = "Source", names_from = group_v, values_from = "prop") %>%
    column_to_rownames("Source")
  return(freq)
}

# get cell freq from metadata table
get_cell_freq2 <- function(cell_info=NULL,
                           sam_info = NULL,
                           cell_prefix = NULL,
                           groupby = "seurat_clusters") {
  require(dplyr)
  sam_info = sam_info %>%
    mutate(cell.filtered = as.numeric(cell.filtered))

  # get freq
  freq = cell_info %>%
    group_by(!! sym(groupby), Source) %>%
    summarise(n = n()) %>%
    mutate(n = as.numeric(n)) %>%
    # change name of groupby
    dplyr::rename(group_v = !! sym(groupby))

  # add samples that do not exist
  tmp = lapply(unique(freq$group_v), FUN = function(x){
    tmp_sources = freq %>% filter(group_v == x) %>% pull(Source) %>% as.character()
    diff_source = as.character(setdiff(sam_info$Source, tmp_sources))
    #print(class(diff_source))
    if (length(diff_source) == 0) {
      return(data.frame())
    }
    return(data.frame(group_v = x, Source = diff_source, n = 0))
  })
  tmp = do.call(rbind, tmp)

  # bind
  freq = rbind(freq, tmp)
  # combine sample info
  freq = freq %>%
    left_join(sam_info, by = "Source") %>%
    mutate(prop = n/cell.filtered)
  # add prefix to groupby
  if (!is.null(cell_prefix)) {
    freq = freq %>%
      mutate(group_v = paste(cell_prefix, group_v, sep = "-"))
  }
  # wider
  freq = freq %>%
    pivot_wider(id_cols = "Source", names_from = group_v, values_from = "prop") %>%
    column_to_rownames("Source")
  return(freq)
}


# plot
plot_freq_cor <- function(cor_mat = NULL,
                          ann_col = NA,
                          col_ann = NA) {
  require(pheatmap)
  paletteLength <- 50
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(seq(min(cor_mat), 0, length.out=ceiling(paletteLength/2) + 1),
                seq(max(cor_mat)/paletteLength, max(cor_mat),
                    length.out=floor(paletteLength/2)))
  # Plot the heatmap
  x = pheatmap::pheatmap(cor_mat, clustering_method = "average",
                         color=myColor, breaks=myBreaks,
                         annotation_col = ann_col, annotation_colors = col_ann)
  return(x)
}


