import csv
 # TODO RA:  Is this needed. Also feels like this data could be stored in a better format. 

val1 = ['3-{[(4-methylphenyl)sulfonyl]amino}propyl pyridin-4-ylcarbamate', '4-(1h-Imidazol-4-Yl)-3-(5-Ethyl-2,4-Dihydroxy-Phenyl)-1h-Pyrazole',
        'ONO-2952', 'Asciminib', 'Lifirafenib', 'Metoserpate', 'Aminophenazone', 'Pirfenidone', 'Pelitinib',
        'PF-05089771', 'Moxaverine', 'Azosemide', '5-[3-(2-METHOXYPHENYL)-1H-PYRROLO[2,3-B]PYRIDIN-5-YL]-N,N-DIMETHYLPYRIDINE-3-CARBOXAMIDE',
        'Mefloquine', 'Sofinicline', '4-(2-(3-Propyl-(1,2,4)oxadiazol-5-yl)-vinyl)-benzene-1,2-diol', 'Capadenoson',
        'Azintamide', 'Carbimazole', 'Tenonitrozole', '5-(5-chloro-2,4-dihydroxyphenyl)-N-ethyl-4-[4-(morpholin-4-ylmethyl)phenyl]isoxazole-3-carboxamide',
        '3,5-Diaminophthalhydrazide', 'Pirlindole', '5,8-dimethoxy-1,4-dimethylquinolin-2(1H)-one',
        'N-[2-(6-AMINO-4-METHYLPYRIDIN-2-YL)ETHYL]-4-CYANOBENZAMIDE', '1-[2-(3-ACETYL-2-HYDROXY-6-METHOXY-PHENYL)-CYCLOPROPYL]-3-(5-CYANO-PYRIDIN-2-YL)-THIOUREA',
        'Tipifarnib', 'Balipodect', 'PF-00356231', 'N-{3-[(7ar,12as,12bs)-7-Oxo-1,3,4,6,7,7a,12a,12b-Octahydroindolo[2,3-a]Quinolizin-12(2h)-Yl]Propyl}Propane-2-Sulfonamide',
        '4-METHYL-PENTANOIC ACID {1-[4-GUANIDINO-1-(THIAZOLE-2-CARBONYL)-BUTYLCARBAMOYL]-2-METHYL-PROPYL}-AMIDE',
        'Cyclo(his-pro)', 'Toreforant', 'Pruvanserin', '[PHENYLALANINYL-PROLINYL]-[2-(PYRIDIN-4-YLAMINO)-ETHYL]-AMINE',
        'Disperse Blue 106', 'Ispronicline', '5-{4-[(3,5-DIFLUOROBENZYL)AMINO]PHENYL}-6-ETHYLPYRIMIDINE-2,4-DIAMINE',
        '4-{4-[4-(3-AMINOPROPOXY)PHENYL]-1H-PYRAZOL-5-YL}-6-CHLOROBENZENE-1,3-DIOL', 'Fipamezole',
        '5-Bromonicotinamide', '3-methyl-N-(pyridin-4-ylmethyl)imidazo[1,2-a]pyrazin-8-amine', '5-Methylpyrrole',
        'AZD-9496', 'Zuranolone', '5-[1-(4-methoxyphenyl)-1H-benzimidazol-6-yl]-1,3,4-oxadiazole-2(3H)-thione',
        'AZD-1236', 'Naratriptan', 'Harmaline']

test1 = ['Rezatomidine', 'Florasulam', '9-HYDROXY-6-(3-HYDROXYPROPYL)-4-(2-METHOXYPHENYL)PYRROLO[3,4-C]CARBAZOLE-1,3(2H,6H)-DIONE',
         'Dasolampanel etibutil', 'N-[2-(5-methyl-4H-1,2,4-triazol-3-yl)phenyl]-7H-pyrrolo[2,3-d]pyrimidin-4-amine',
         'Verucerfont', 'Bunazosin', '(2R)-4-[(8R)-8-METHYL-2-(TRIFLUOROMETHYL)-5,6-DIHYDRO[1,2,4]TRIAZOLO[1,5-A]PYRAZIN-7(8H)-YL]-4-OXO-1-(2,4,5-TRIFLUOROPHENYL)BUTAN-2-AMINE',
         'Etrasimod', 'GSK-3117391', 'Clortermine', 'Epibatidine', 'Molidustat', 'Oxaprozin', 'Vebreltinib',
         '2-Isobutyl-3-methoxypyrazine', 'CB-103', '7-Methoxy-8-[1-(Methylsulfonyl)-1h-Pyrazol-4-Yl]Naphthalene-2-Carboximidamide',
         'Fenpyroximate', 'Trilaciclib', '1-(5-BROMO-PYRIDIN-2-YL)-3-[2-(6-FLUORO-2-HYDROXY-3-PROPIONYL-PHENYL)-CYCLOPROPYL]-UREA',
         'Ifidancitinib', 'AMG-319', 'Gandotinib', 'Navoximod', 'Xylose-Derived Imidazole', 'Oglemilast', 'AZD-5069', 'Dapiprazole',
         'RAD-140', 'N-[(1R)-3-(4-HYDROXYPHENYL)-1-METHYLPROPYL]-2-(2-PHENYL-1H-INDOL-3-YL)ACETAMIDE', 'Cefatrizine',
         'N-((1R,2R)-2-(5-CHLORO-1H-INDOLE-2-CARBOXAMIDO)CYCLOHEXYL)-5-METHYL-4,5,6,7-TETRAHYDROTHIAZOLO[5,4-C]PYRIDINE-2-CARBOXAMIDE',
         'Sapanisertib', 'Aptazapine', 'PX-12', 'Vosaroxin', '4-Iodopyrazole', '2,4-Diamino-6-Phenyl-5,6,7,8,-Tetrahydropteridine',
         'N,N-DIMETHYL(5-(PYRIDIN-3-YL)FURAN-2-YL)METHANAMINE', "N-{(3R,4S)-4-[(6-amino-4-methylpyridin-2-yl)methyl]pyrrolidin-3-yl}-N'-(3-chlorobenzyl)ethane-1,2-diamine",
         'D-157495', 'N-{2-methyl-5-[(6-phenylpyrimidin-4-yl)amino]phenyl}methanesulfonamide', 'Ethionamide',
         'GSK-2881078', 'PF-06282999', 'Vorolanib', 'BMS-919373']

train1 = ['BTRX-246040', 'Bendamustine', 'Letrazuril', 'HYDROXY[3-(6-METHYLPYRIDIN-2-YL)PROPYL]FORMAMIDE',
          '4-(2-(1H-IMIDAZOL-4-YL)ETHYLAMINO)-2-(PHENYLAMINO)PYRAZOLO[1,5-A][1,3,5]TRIAZINE-8-CARBONITRILE',
          'Nelotanserin', 'Cicletanine', '5-CHLORO-THIOPHENE-2-CARBOXYLIC ACID ((3S,4S)-4-FLUORO- 1-{[2-FLUORO-4-(2-OXO-2H-PYRIDIN-1-YL)-PHENYLCARBAMOYL]-METHYL}-PYRROLIDIN-3-YL)-AMIDE',
          'Zilpaterol', 'Duvelisib', 'F-15599', 'Pindolol', 'Reldesemtiv', 'Methylisothiazolinone', 'MK-8245',
          'Brensocatib', '2-Amino-6-Chloropyrazine', '(2-AMINO-1,3-OXAZOL-5-YL)-(3-BROMOPHENYL)METHANONE', 'Foslinanib',
          'GSK-356278', 'Lanabecestat', '9-Deazahypoxanthine', 'INCB-057643', 'Surinabant', 'PXT 3003',
          '3-(6-HYDROXY-NAPHTHALEN-2-YL)-BENZO[D]ISOOXAZOL-6-OL', 'Oltipraz', 'ONO-8539', 'Chlorzoxazone', 'Deferiprone',
          'SB-705498', 'Nidufexor', '7-amino-2-tert-butyl-4-{[2-(1H-imidazol-4-yl)ethyl]amino}pyrido[2,3-d]pyrimidine-6-carboxamide',
          'Imiglitazar', 'Rimacalib', 'Tetrazolyl Histidine', 'JNJ-54175446', 'Foliglurax', 'CH-5132799',
          '4-{2-[(7-amino-2-furan-2-yl[1,2,4]triazolo[1,5-a][1,3,5]triazin-5-yl)amino]ethyl}phenol',
          'N-[(2R)-2-{[(2S)-2-(1,3-benzoxazol-2-yl)pyrrolidin-1-yl]carbonyl}hexyl]-N-hydroxyformamide',
          'Acid yellow 54 free acid', '3-phenyl-5-(1H-pyrazol-3-yl)isoxazole', 'MK-5108', 'GLPG-1205', 'Pioglitazone',
          '4-[1-allyl-7-(trifluoromethyl)-1H-indazol-3-yl]benzene-1,3-diol', 'MK-212', 'Bamifylline',
          '2-(Sec-Butyl)Thiazole', 'Rivanicline', '5-benzyl-1,3-thiazol-2-amine', 'Sitamaquine',
          '(2R,3R)-N^1^-[(1S)-2,2-DIMETHYL-1-(METHYLCARBAMOYL)PROPYL]-N^4^-HYDROXY-2-(2-METHYLPROPYL)-3-{[(1,3-THIAZOL-2-YLCARBONYL)AMINO]METHYL}BUTANEDIAMIDE',
          'Cinchocaine', 'Simurosertib', 'BTRX-335140', 'Afuresertib', 'Enitociclib', 'CCX-140',
          "N-{(3S,4S)-4-[(6-AMINO-4-METHYLPYRIDIN-2-YL)METHYL]PYRROLIDIN-3-YL}-N'-(4-CHLOROBENZYL)ETHANE-1,2-DIAMINE",
          '(7as,12ar,12bs)-1,2,3,4,7a,12,12a,12b-Octahydroindolo[2,3-a]Quinolizin-7(6h)-One', 'CAN-508', 'Cinalukast',
          '(1s,2s)-1-Amino-1-(1,3-Thiazol-2-Yl)Propan-2-Ol', 'Zolpidem', 'Galeterone', 'Mercaptopurine', 'GZ-389988',
          'Varespladib methyl', 'Farampator', 'PD-173952', 'Pirbuterol', 'Amsacrine', 'AZD-0328', 'Umifenovir',
          '(S)-wiskostatin', '4-(6-{[(1R)-1-(hydroxymethyl)propyl]amino}imidazo[1,2-b]pyridazin-3-yl)benzoic acid',
          'CC-115', 'Protionamide', 'Linrodostat', 'LGH-447', 'Henatinib', 'Alprazolam', 'MK-886', 'MK-6186',
          '6-[2-(1H-INDOL-6-YL)ETHYL]PYRIDIN-2-AMINE', '3-(4-{2-[2-(2-bromo-acetylamino)-ethyldisulfanyl]-ethylcarbamoyl}-cyclohexylcarbamoyl)-pyrazine-2-carboxylic acid',
          '6-MORPHOLIN-4-YL-9H-PURINE', 'N-hydroxy-5-[(3-phenyl-5,6-dihydroimidazo[1,2-a]pyrazin-7(8H)-yl)carbonyl]thiophene-2-carboxamide',
          '5-(Aminomethyl)-6-(2,4-Dichlorophenyl)-2-(3,5-Dimethoxyphenyl)Pyrimidin-4-Amine', 'Tulrampator', 'JAB-3068',
          'Ganaplacide', 'Prenoxdiazine', '1-[(2S)-4-(5-phenyl-1H-pyrazolo[3,4-b]pyridin-4-yl)morpholin-2-yl]methanamine',
          'Trapidil', '4-(2-Thienyl)-1-(4-Methylbenzyl)-1h-Imidazole', 'Elvitegravir', 'Glumetinib', 'Fenamole',
          'Lonidamine', 'CXD101', 'GSK-2982772', 'Olomoucine', 'R-82913', 'Etodolac', 'Zonisamide', 'Amlexanox',
          'PF-05212377', 'Arotinolol', 'Thiacloprid', 'N-(2-Aminoethyl)-5-Chloroisoquinoline-8-Sulfonamide', 'Thiohexam',
          'LY-2811376', 'BOS172722', 'N-[2-(METHYLAMINO)ETHYL]-5-ISOQUINOLINESULFONAMIDE', 'Samuraciclib', 'Amuvatinib',
          'Pecavaptan', 'Bradanicline', '2-[(2-methoxy-5-methylphenoxy)methyl]pyridine', 'AZD-1940', 'Belotecan',
          'Cenobamate', '2-(1H-pyrrol-1-ylcarbonyl)benzene-1,3,5-triol', 'Latrepirdine', 'Carboxymethylthio-3-(3-Chlorophenyl)-1,2,4-Oxadiazol',
          'Terevalefim', 'Chlormidazole', 'Rupatadine', 'Tolmetin', 'Lanifibranor',
          'N1-CYCLOPENTYL-N2-(THIAZOL-2-YL)OXALAMIDE', 'Benzimidazole', 'GSK-239512', 'PF-06260414', 'Etomidate',
          'Azatadine', 'RO-5028442', 'Abiraterone', 'arazepide', 'Upadacitinib', '2X-121', 'CRA_10972', 'ABX-464',
          'Tegoprazan', 'Tizanidine', '4-{[5-chloro-4-(1H-indol-3-yl)pyrimidin-2-yl]amino}-N-ethylpiperidine-1-carboxamide',
          'Anastrozole', 'Uracil', 'Nedisertib', 'Cilostazol', '1-(3,5-DICHLOROPHENYL)-5-METHYL-1H-1,2,4-TRIAZOLE-3-CARBOXYLIC ACID',
          'Fluconazole', 'Ipatasertib', 'N-acetylhistamine', 'Tinostamustine', 'Setipiprant', 'R-1487', 'GDC-0134',
          'Tolimidone', '(2E)-N-hydroxy-3-[1-methyl-4-(phenylacetyl)-1H-pyrrol-2-yl]prop-2-enamide',
          '3-[(2,2-DIMETHYLPROPANOYL)AMINO]-N-1,3-THIAZOL-2-YLPYRIDINE-2-CARBOXAMIDE', 'Pozanicline', 'Veliparib',
          'Ondelopran', 'IQP-0528', 'Ranirestat', '7-Nitroindazole', 'Phosphonoacetohydroxamic Acid', 'Perampanel',
          'N-phenyl-1H-pyrrolo[2,3-b]pyridin-3-amine', 'Tazobactam', 'N~4~-methyl-N~4~-(3-methyl-1H-indazol-6-yl)-N~2~-(3,4,5-trimethoxyphenyl)pyrimidine-2,4-diamine',
          'Mercaptocarboxylate Inhibitor', 'LY-3200882', 'Flosequinan', 'PF-04691502', 'GSK-2018682', 'Cefazedone',
          '{3-[(5-CHLORO-1,3-BENZOTHIAZOL-2-YL)METHYL]-2,4-DIOXO-3,4-DIHYDROPYRIMIDIN-1(2H)-YL}ACETIC ACID',
          'Selgantolimod', '2,2,2-TRIFLUORO-1-{5-[(3-PHENYL-5,6-DIHYDROIMIDAZO[1,2-A]PYRAZIN-7(8H)-YL)CARBONYL]THIOPHEN-2-YL}ETHANE-1,1-DIOL',
          'Isradipine', 'Gluco-Phenylimidazole', "TRW3-(2-AMINO-3-HYDROXY-PROPYL)-6-(N'-CYCLOHEXYL-HYDRAZINO)OCTAHYDRO-INDOL-7-OL",
          'Tucidinostat', 'Pemigatinib', '5-CHLORO-6-METHYL-N-(2-PHENYLETHYL)-2-PYRIDIN-2-YLPYRIMIDIN-4-AMINE',
          '[1-(4-Fluorobenzyl)Cyclobutyl]Methyl (1s)-1-[Oxo(1h-Pyrazol-5-Ylamino)Acetyl]Pentylcarbamate',
          '(1S,3R,6S)-4-oxo-6-{4-[(2-phenylquinolin-4-yl)methoxy]phenyl}-5-azaspiro[2.4]heptane-1-carboxylic acid',
          '1-[(2S)-4-(5-BROMO-1H-PYRAZOLO[3,4-B]PYRIDIN-4-YL)MORPHOLIN-2-YL]METHANAMINE',
          'ethyl 3-[(E)-2-amino-1-cyanoethenyl]-6,7-dichloro-1-methyl-1H-indole-2-carboxylate',
          '4-(5,11-DIOXO-5H-INDENO[1,2-C]ISOQUINOLIN-6(11H)-YL)BUTANOATE', '2-Methoxy-3-Isopropylpyrazine',
          '3-(1-NAPHTHYLMETHOXY)PYRIDIN-2-AMINE', 'Di-2-pyridylketone 4-cyclohexyl-4-methyl-3-thiosemicarbazone',
          'Methimazole', 'S-777469', 'ORE-1001', 'Oliceridine', 'N-6022', 'GDC-0425',
          "4-Bromo-3-(5'-Carboxy-4'-Chloro-2'-Fluorophenyl)-1-Methyl-5-Trifluoromethyl-Pyrazol", 'Tivirapine',
          'Rufinamide', '6-BENZYL-1-BENZYLOXYMETHYL-5-ISOPROPYL URACIL', 'Avadomide',
          '1-[(4S)-4-amino-5-(1,3-benzothiazol-2-yl)-5-oxopentyl]guanidine', 'Cipargamin',
          "1-(6-CYANO-3-PYRIDYLCARBONYL)-5',8'-DIFLUOROSPIRO[PIPERIDINE-4,2'(1'H)-QUINAZOLINE]-4'-AMINE", 'JNJ-39393406',
          'Piribedil', 'Imidazole', 'Pipequaline', 'Lavoltidine', 'N-{3-[5-(1H-1,2,4-triazol-3-yl)-1H-indazol-3-yl]phenyl}furan-2-carboxamide',
          'Tinengotinib', 'Uracil C-13', 'Minodronic acid', 'Tovorafenib', 'PF-03463275',
          '3-(4-CHLOROPHENYL)-5-(METHYLTHIO)-4H-1,2,4-TRIAZOLE', 'AZD-1981', 'Revaprazan', 'Rabeximod',
          '3-Methyladenine', '5-Aminoisoquinoline', '6-HYDROXY-1,3-BENZOTHIAZOLE-2-SULFONAMIDE', '1,2,4-Triazole',
          '1-(2-Chlorophenyl)-3,5-Dimethyl-1h-Pyrazole-4-Carboxylic Acid Ethyl Ester',
          '6-CHLORO-9-HYDROXY-1,3-DIMETHYL-1,9-DIHYDRO-4H-PYRAZOLO[3,4-B]QUINOLIN-4-ONE', 'Dilmapimod', 'CGS-27023',
          'N-[2-methyl-5-(methylcarbamoyl)phenyl]-2-{[(1R)-1-methylpropyl]amino}-1,3-thiazole-5-carboxamide',
          'Nitenpyram', 'AMG-900', '6-CHLORO-4-(CYCLOHEXYLOXY)-3-PROPYLQUINOLIN-2(1H)-ONE', 'Roflumilast', 'Posizolid',
          'Carboxyamidotriazole', '5-(2-chlorophenyl)-1,3,4-thiadiazole-2-sulfonamide',
          '4-CHLORO-6-(4-PIPERAZIN-1-YL-1H-PYRAZOL-5-YL)BENZENE-1,3-DIOL',
          '(2R)-1-[(5,6-DIPHENYL-7H-PYRROLO[2,3-D]PYRIMIDIN-4-YL)AMINO]PROPAN-2-OL', 'Birabresib', 'Cimicoxib',
          'Snubh-nm-333 F-18', 'Polaprezinc', 'Pradigastat', "(4R)-7-chloro-9-methyl-1-oxo-1,2,4,9-tetrahydrospiro[beta-carboline-3,4'-piperidine]-4-carbonitrile",
          'Paltusotine', 'Diclazuril', 'Butalamine', 'Riluzole', 'N-METHYL-1-[4-(9H-PURIN-6-YL)PHENYL]METHANAMINE',
          'Aleplasinin', 'Linzagolix', 'Gimeracil', 'BMS-488043', 'Clonazolam', 'Ciclopirox', 'Naphthoquine',
          '4-[(3-BROMO-4-O-SULFAMOYLBENZYL)(4-CYANOPHENYL)AMINO]-4H-[1,2,4]-TRIAZOLE',
          '5-[(Z)-(5-Chloro-2-oxo-1,2-dihydro-3H-indol-3-ylidene)methyl]-N,2,4-trimethyl-1H-pyrrole-3-carboxamide',
          'AG-24322', 'Pamiparib', 'PF-05241328', 'Dalfampridine', '2-Pyridinethiol', 'Mivavotinib',
          '2-(4-ETHYLPIPERAZIN-1-YL)-4-(PHENYLAMINO)PYRAZOLO[1,5-A][1,3,5]TRIAZINE-8-CARBONITRILE', 'Metralindole',
          'Lamotrigine', 'Cinoxacin', 'Anagliptin', 'Losartan', 'Allopurinol', 'Frovatriptan', 'Deucravacitinib',
          '2-(2-chloropyridin-4-yl)-4-methyl-1H-isoindole-1,3(2H)-dione', 'Irbesartan',
          '3-[2-bromo-4-(1H-pyrazolo[3,4-c]pyridazin-3-ylmethyl)phenoxy]-5-methylbenzonitrile', 'Onatasertib',
          "N'-(5-CHLORO-1,3-BENZODIOXOL-4-YL)-N-(3-MORPHOLIN-4-YLPHENYL)PYRIMIDINE-2,4-DIAMINE", 'Dactolisib',
          'Capivasertib', 'Cenerimod', 'PTI-428', 'V116517', 'AMD-070',
          '6-(2,4-DIAMINO-6-ETHYLPYRIMIDIN-5-YL)-4-(3-METHOXYPROPYL)-2,2-DIMETHYL-2H-1,4-BENZOXAZIN-3(4H)-ONE',
          'Vincamine', '(2Z)-1-(5-Chloro-1H-indol-3-yl)-3-hydroxy-3-(1H-tetrazol-5-yl)-2-propen-1-one', 'Sisunatovir',
          'Quinacrine mustard', 'Verdiperstat', 'Celecoxib', '3-{3-[(DIMETHYLAMINO)METHYL]-1H-INDOL-7-YL}PROPAN-1-OL',
          'Ethacridine', '(2S)-2-AMINO-4-(METHYLSULFANYL)-1-(1,3-THIAZOL-2-YL)BUTANE-1,1-DIOL',
          'N-(2,6-dimethylphenyl)-5-phenylimidazo[1,5-a]pyrazin-8-amine', 'Ezutromid', 'Nedocromil', 'Ziresovir',
          'Tavapadon', 'Pardoprunox', 'K-134', 'Pyrrole-2-Carboxylate', 'Sertaconazole', '6-Methylpurine',
          '2-amino-5-[3-(1-ethyl-1H-pyrazol-5-yl)-1H-pyrrolo[2,3-b]pyridin-5-yl]-N,N-dimethylbenzamide', 'Linopirdine',
          'Vecabrutinib', 'Venglustat', '4-{[1-Methyl-5-(2-Methyl-Benzoimidazol-1-Ylmethyl)-1h-Benzoimidazol-2-Ylmethyl]-Amino}-Benzamidine',
          'N,4-dimethyl-3-[(1-phenyl-1H-pyrazolo[3,4-d]pyrimidin-4-yl)amino]benzamide', 'NP-G2-044', 'Olutasidenib',
          'Simpinicline', 'PP-121', 'PHENYLAMINOIMIDAZO(1,2-ALPHA)PYRIDINE', 'GW810781', 'PSI-697',
          'N-(4-{[(3S)-3-(dimethylamino)pyrrolidin-1-yl]carbonyl}phenyl)-5-fluoro-4-[2-methyl-1-(1-methylethyl)-1H-imidazol-5-yl]pyrimidin-2-amine',
          'Hydronidone', 'Pyrazole', 'Bromazepam', 'AZD-8418', '4-({[4-(3-METHYLBENZOYL)PYRIDIN-2-YL]AMINO}METHYL)BENZENECARBOXIMIDAMIDE',
          'Fanapanel', 'Dianicline', '[2-(5-Mercapto-[1,3,4]thiadiazol-2-ylcarbamoyl)-1-phenyl-ethyl]-carbamic acid benzyl ester',
          "(4R)-7,8-dichloro-1',9-dimethyl-1-oxo-1,2,4,9-tetrahydrospiro[beta-carboline-3,4'-piperidine]-4-carbonitrile",
          'E-7820', '7-(aminomethyl)-6-(2-chlorophenyl)-1-methyl-1H-benzimidazole-5-carbonitrile', 'AZD-1386',
          'Imaradenant', 'Pranoprofen', 'Tisopurine', 'Decernotinib', 'Triazolopyridine',
          '9-ACETYL-2,3,4,9-TETRAHYDRO-1H-CARBAZOL-1-ONE', 'PF-06700841', 'Rogaratinib', 'ABT-288', 'Buparlisib',
          '4-chloro-6-{5-[(2-morpholin-4-ylethyl)amino]-1,2-benzisoxazol-3-yl}benzene-1,3-diol', 'Cyanocinnoline',
          'Zamicastat', 'Cp403700, (S)-1-{2-[(5-Chloro-1h-Indole-2-Carbonyl)-Amino]-3-Phenyl-Propionyl}-Azetidine-3-Carboxylate',
          'Isatoic anhydride', 'UK-500001', 'Droxicam', 'LY-3023414', 'Urapidil',
          '5-(7-(4-(4,5-dihydro-2-oxazolyl)phenoxy)heptyl)-3-methyl isoxazole', '2-CHLORO-N-[(3R)-2-OXO-1,2,3,4-TETRAHYDROQUINOLIN-3-YL]-6H-THIENO[2,3-B]PYRROLE-5-CARBOXAMIDE',
          'Triapine', 'PF-04991532', '1-[1-(3-aminophenyl)-3-tert-butyl-1H-pyrazol-5-yl]-3-phenylurea',
          'Benzoylformic Acid', 'N~2~-1,3-BENZOXAZOL-2-YL-3-CYCLOHEXYL-N-{2-[(4-METHOXYPHENYL)AMINO]ETHYL}-L-ALANINAMIDE',
          'PF-05175157', '1,N6-Ethenoadenine', '6-[4-(2-fluorophenyl)-1,3-oxazol-5-yl]-N-(1-methylethyl)-1,3-benzothiazol-2-amine',
          'Rilematovir', '(2S)-2-{[3-(3-aminophenyl)imidazo[1,2-b]pyridazin-6-yl]amino}-3-methylbutan-1-ol',
          'Pantoprazole', 'Piclozotan', 'LGD2941', 'Ocinaplon', 'Theodrenaline', 'Pexacerfont', 'Acumapimod',
          'Pumosetrag', '2-(methylsulfanyl)-5-(thiophen-2-ylmethyl)-1H-imidazol-4-ol', 'TAK-831',
          '4-{4-[(5-hydroxy-2-methylphenyl)amino]quinolin-7-yl}-1,3-thiazole-2-carbaldehyde', 'Cadralazine',
          'GSK-1059615', '(2R)-1-{4-[(4-Anilino-5-bromo-2-pyrimidinyl)amino]phenoxy}-3-(dimethylamino)-2-propanol',
          'Harmine', 'Proquazone', 'Viminol', 'GSK2798745', 'Eletriptan', 'Imiquimod', 'Losmapimod']

val2 = ['BTRX-246040', 'Bendamustine', 'Letrazuril', 'HYDROXY[3-(6-METHYLPYRIDIN-2-YL)PROPYL]FORMAMIDE',
        '4-(2-(1H-IMIDAZOL-4-YL)ETHYLAMINO)-2-(PHENYLAMINO)PYRAZOLO[1,5-A][1,3,5]TRIAZINE-8-CARBONITRILE',
        'Nelotanserin', 'Cicletanine',
        '5-CHLORO-THIOPHENE-2-CARBOXYLIC ACID ((3S,4S)-4-FLUORO- 1-{[2-FLUORO-4-(2-OXO-2H-PYRIDIN-1-YL)-PHENYLCARBAMOYL]-METHYL}-PYRROLIDIN-3-YL)-AMIDE',
        'Zilpaterol', 'Duvelisib', 'F-15599', 'Pindolol', 'Reldesemtiv', 'Methylisothiazolinone', 'MK-8245',
        'Brensocatib', '2-Amino-6-Chloropyrazine', '(2-AMINO-1,3-OXAZOL-5-YL)-(3-BROMOPHENYL)METHANONE', 'Foslinanib',
        'GSK-356278', 'Lanabecestat', '9-Deazahypoxanthine', 'INCB-057643', 'Surinabant', 'PXT 3003',
        '3-(6-HYDROXY-NAPHTHALEN-2-YL)-BENZO[D]ISOOXAZOL-6-OL', 'Oltipraz', 'ONO-8539', 'Chlorzoxazone', 'Deferiprone',
        'SB-705498', 'Nidufexor',
        '7-amino-2-tert-butyl-4-{[2-(1H-imidazol-4-yl)ethyl]amino}pyrido[2,3-d]pyrimidine-6-carboxamide', 'Imiglitazar',
        'Rimacalib', 'Tetrazolyl Histidine', 'JNJ-54175446', 'Foliglurax', 'CH-5132799',
        '4-{2-[(7-amino-2-furan-2-yl[1,2,4]triazolo[1,5-a][1,3,5]triazin-5-yl)amino]ethyl}phenol',
        'N-[(2R)-2-{[(2S)-2-(1,3-benzoxazol-2-yl)pyrrolidin-1-yl]carbonyl}hexyl]-N-hydroxyformamide',
        'Acid yellow 54 free acid', '3-phenyl-5-(1H-pyrazol-3-yl)isoxazole', 'MK-5108', 'GLPG-1205', 'Pioglitazone',
        '4-[1-allyl-7-(trifluoromethyl)-1H-indazol-3-yl]benzene-1,3-diol', 'MK-212', 'Bamifylline'
        ]

test2 = ['2-(Sec-Butyl)Thiazole', 'Rivanicline', '5-benzyl-1,3-thiazol-2-amine', 'Sitamaquine',
         '(2R,3R)-N^1^-[(1S)-2,2-DIMETHYL-1-(METHYLCARBAMOYL)PROPYL]-N^4^-HYDROXY-2-(2-METHYLPROPYL)-3-{[(1,3-THIAZOL-2-YLCARBONYL)AMINO]METHYL}BUTANEDIAMIDE',
         'Cinchocaine', 'Simurosertib', 'BTRX-335140', 'Afuresertib', 'Enitociclib', 'CCX-140',
         "N-{(3S,4S)-4-[(6-AMINO-4-METHYLPYRIDIN-2-YL)METHYL]PYRROLIDIN-3-YL}-N'-(4-CHLOROBENZYL)ETHANE-1,2-DIAMINE",
         '(7as,12ar,12bs)-1,2,3,4,7a,12,12a,12b-Octahydroindolo[2,3-a]Quinolizin-7(6h)-One', 'CAN-508', 'Cinalukast',
         '(1s,2s)-1-Amino-1-(1,3-Thiazol-2-Yl)Propan-2-Ol', 'Zolpidem', 'Galeterone', 'Mercaptopurine', 'GZ-389988',
         'Varespladib methyl', 'Farampator', 'PD-173952', 'Pirbuterol', 'Amsacrine', 'AZD-0328', 'Umifenovir',
         '(S)-wiskostatin', '4-(6-{[(1R)-1-(hydroxymethyl)propyl]amino}imidazo[1,2-b]pyridazin-3-yl)benzoic acid',
         'CC-115', 'Protionamide', 'Linrodostat', 'LGH-447', 'Henatinib', 'Alprazolam', 'MK-886', 'MK-6186',
         '6-[2-(1H-INDOL-6-YL)ETHYL]PYRIDIN-2-AMINE',
         '3-(4-{2-[2-(2-bromo-acetylamino)-ethyldisulfanyl]-ethylcarbamoyl}-cyclohexylcarbamoyl)-pyrazine-2-carboxylic acid',
         '6-MORPHOLIN-4-YL-9H-PURINE',
         'N-hydroxy-5-[(3-phenyl-5,6-dihydroimidazo[1,2-a]pyrazin-7(8H)-yl)carbonyl]thiophene-2-carboxamide',
         '5-(Aminomethyl)-6-(2,4-Dichlorophenyl)-2-(3,5-Dimethoxyphenyl)Pyrimidin-4-Amine', 'Tulrampator', 'JAB-3068',
         'Ganaplacide', 'Prenoxdiazine',
         '1-[(2S)-4-(5-phenyl-1H-pyrazolo[3,4-b]pyridin-4-yl)morpholin-2-yl]methanamine', 'Trapidil',
         '4-(2-Thienyl)-1-(4-Methylbenzyl)-1h-Imidazole'
         ]

train2 = ['Elvitegravir', 'Glumetinib', 'Fenamole', 'Lonidamine', 'CXD101', 'GSK-2982772', 'Olomoucine', 'R-82913',
          'Etodolac', 'Zonisamide', 'Amlexanox', 'PF-05212377', 'Arotinolol', 'Thiacloprid',
          'N-(2-Aminoethyl)-5-Chloroisoquinoline-8-Sulfonamide', 'Thiohexam', 'LY-2811376', 'BOS172722',
          'N-[2-(METHYLAMINO)ETHYL]-5-ISOQUINOLINESULFONAMIDE', 'Samuraciclib', 'Amuvatinib', 'Pecavaptan',
          'Bradanicline', '2-[(2-methoxy-5-methylphenoxy)methyl]pyridine', 'AZD-1940', 'Belotecan', 'Cenobamate',
          '2-(1H-pyrrol-1-ylcarbonyl)benzene-1,3,5-triol', 'Latrepirdine',
          'Carboxymethylthio-3-(3-Chlorophenyl)-1,2,4-Oxadiazol', 'Terevalefim', 'Chlormidazole', 'Rupatadine',
          'Tolmetin', 'Lanifibranor', 'N1-CYCLOPENTYL-N2-(THIAZOL-2-YL)OXALAMIDE', 'Benzimidazole', 'GSK-239512',
          'PF-06260414', 'Etomidate', 'Azatadine', 'RO-5028442', 'Abiraterone', 'arazepide', 'Upadacitinib', '2X-121',
          'CRA_10972', 'ABX-464', 'Tegoprazan', 'Tizanidine',
          '4-{[5-chloro-4-(1H-indol-3-yl)pyrimidin-2-yl]amino}-N-ethylpiperidine-1-carboxamide', 'Anastrozole',
          'Uracil', 'Nedisertib', 'Cilostazol', '1-(3,5-DICHLOROPHENYL)-5-METHYL-1H-1,2,4-TRIAZOLE-3-CARBOXYLIC ACID',
          'Fluconazole', 'Ipatasertib', 'N-acetylhistamine', 'Tinostamustine', 'Setipiprant', 'R-1487', 'GDC-0134',
          'Tolimidone', '(2E)-N-hydroxy-3-[1-methyl-4-(phenylacetyl)-1H-pyrrol-2-yl]prop-2-enamide',
          '3-[(2,2-DIMETHYLPROPANOYL)AMINO]-N-1,3-THIAZOL-2-YLPYRIDINE-2-CARBOXAMIDE', 'Pozanicline', 'Veliparib',
          'Ondelopran', 'IQP-0528', 'Ranirestat', '7-Nitroindazole', 'Phosphonoacetohydroxamic Acid', 'Perampanel',
          'N-phenyl-1H-pyrrolo[2,3-b]pyridin-3-amine', 'Tazobactam',
          'N~4~-methyl-N~4~-(3-methyl-1H-indazol-6-yl)-N~2~-(3,4,5-trimethoxyphenyl)pyrimidine-2,4-diamine',
          'Mercaptocarboxylate Inhibitor', 'LY-3200882', 'Flosequinan', 'PF-04691502', 'GSK-2018682', 'Cefazedone',
          '{3-[(5-CHLORO-1,3-BENZOTHIAZOL-2-YL)METHYL]-2,4-DIOXO-3,4-DIHYDROPYRIMIDIN-1(2H)-YL}ACETIC ACID',
          'Selgantolimod',
          '2,2,2-TRIFLUORO-1-{5-[(3-PHENYL-5,6-DIHYDROIMIDAZO[1,2-A]PYRAZIN-7(8H)-YL)CARBONYL]THIOPHEN-2-YL}ETHANE-1,1-DIOL',
          'Isradipine', 'Gluco-Phenylimidazole',
          "TRW3-(2-AMINO-3-HYDROXY-PROPYL)-6-(N'-CYCLOHEXYL-HYDRAZINO)OCTAHYDRO-INDOL-7-OL", 'Tucidinostat',
          'Pemigatinib', '5-CHLORO-6-METHYL-N-(2-PHENYLETHYL)-2-PYRIDIN-2-YLPYRIMIDIN-4-AMINE',
          '[1-(4-Fluorobenzyl)Cyclobutyl]Methyl (1s)-1-[Oxo(1h-Pyrazol-5-Ylamino)Acetyl]Pentylcarbamate',
          '(1S,3R,6S)-4-oxo-6-{4-[(2-phenylquinolin-4-yl)methoxy]phenyl}-5-azaspiro[2.4]heptane-1-carboxylic acid',
          '1-[(2S)-4-(5-BROMO-1H-PYRAZOLO[3,4-B]PYRIDIN-4-YL)MORPHOLIN-2-YL]METHANAMINE',
          'ethyl 3-[(E)-2-amino-1-cyanoethenyl]-6,7-dichloro-1-methyl-1H-indole-2-carboxylate',
          '4-(5,11-DIOXO-5H-INDENO[1,2-C]ISOQUINOLIN-6(11H)-YL)BUTANOATE', '2-Methoxy-3-Isopropylpyrazine',
          '3-(1-NAPHTHYLMETHOXY)PYRIDIN-2-AMINE', 'Di-2-pyridylketone 4-cyclohexyl-4-methyl-3-thiosemicarbazone',
          'Methimazole', 'S-777469', 'ORE-1001', 'Oliceridine', 'N-6022', 'GDC-0425',
          "4-Bromo-3-(5'-Carboxy-4'-Chloro-2'-Fluorophenyl)-1-Methyl-5-Trifluoromethyl-Pyrazol", 'Tivirapine',
          'Rufinamide', '6-BENZYL-1-BENZYLOXYMETHYL-5-ISOPROPYL URACIL', 'Avadomide',
          '1-[(4S)-4-amino-5-(1,3-benzothiazol-2-yl)-5-oxopentyl]guanidine', 'Cipargamin',
          "1-(6-CYANO-3-PYRIDYLCARBONYL)-5',8'-DIFLUOROSPIRO[PIPERIDINE-4,2'(1'H)-QUINAZOLINE]-4'-AMINE",
          'JNJ-39393406', 'Piribedil', 'Imidazole', 'Pipequaline', 'Lavoltidine',
          'N-{3-[5-(1H-1,2,4-triazol-3-yl)-1H-indazol-3-yl]phenyl}furan-2-carboxamide', 'Tinengotinib', 'Uracil C-13',
          'Minodronic acid', 'Tovorafenib', 'PF-03463275', '3-(4-CHLOROPHENYL)-5-(METHYLTHIO)-4H-1,2,4-TRIAZOLE',
          'AZD-1981', 'Revaprazan', 'Rabeximod', '3-Methyladenine', '5-Aminoisoquinoline',
          '6-HYDROXY-1,3-BENZOTHIAZOLE-2-SULFONAMIDE', '1,2,4-Triazole',
          '1-(2-Chlorophenyl)-3,5-Dimethyl-1h-Pyrazole-4-Carboxylic Acid Ethyl Ester',
          '6-CHLORO-9-HYDROXY-1,3-DIMETHYL-1,9-DIHYDRO-4H-PYRAZOLO[3,4-B]QUINOLIN-4-ONE', 'Dilmapimod', 'CGS-27023',
          'N-[2-methyl-5-(methylcarbamoyl)phenyl]-2-{[(1R)-1-methylpropyl]amino}-1,3-thiazole-5-carboxamide',
          'Nitenpyram', 'AMG-900', '6-CHLORO-4-(CYCLOHEXYLOXY)-3-PROPYLQUINOLIN-2(1H)-ONE', 'Roflumilast', 'Posizolid',
          'Carboxyamidotriazole', '5-(2-chlorophenyl)-1,3,4-thiadiazole-2-sulfonamide',
          '4-CHLORO-6-(4-PIPERAZIN-1-YL-1H-PYRAZOL-5-YL)BENZENE-1,3-DIOL',
          '(2R)-1-[(5,6-DIPHENYL-7H-PYRROLO[2,3-D]PYRIMIDIN-4-YL)AMINO]PROPAN-2-OL', 'Birabresib', 'Cimicoxib',
          'Snubh-nm-333 F-18', 'Polaprezinc', 'Pradigastat',
          "(4R)-7-chloro-9-methyl-1-oxo-1,2,4,9-tetrahydrospiro[beta-carboline-3,4'-piperidine]-4-carbonitrile",
          'Paltusotine', 'Diclazuril', 'Butalamine', 'Riluzole', 'N-METHYL-1-[4-(9H-PURIN-6-YL)PHENYL]METHANAMINE',
          'Aleplasinin', 'Linzagolix', 'Gimeracil', 'BMS-488043', 'Clonazolam', 'Ciclopirox', 'Naphthoquine',
          '4-[(3-BROMO-4-O-SULFAMOYLBENZYL)(4-CYANOPHENYL)AMINO]-4H-[1,2,4]-TRIAZOLE',
          '5-[(Z)-(5-Chloro-2-oxo-1,2-dihydro-3H-indol-3-ylidene)methyl]-N,2,4-trimethyl-1H-pyrrole-3-carboxamide',
          'AG-24322', 'Pamiparib', 'PF-05241328', 'Dalfampridine', '2-Pyridinethiol', 'Mivavotinib',
          '2-(4-ETHYLPIPERAZIN-1-YL)-4-(PHENYLAMINO)PYRAZOLO[1,5-A][1,3,5]TRIAZINE-8-CARBONITRILE', 'Metralindole',
          'Lamotrigine', 'Cinoxacin', 'Anagliptin', 'Losartan', 'Allopurinol', 'Frovatriptan', 'Deucravacitinib',
          '2-(2-chloropyridin-4-yl)-4-methyl-1H-isoindole-1,3(2H)-dione', 'Irbesartan',
          '3-[2-bromo-4-(1H-pyrazolo[3,4-c]pyridazin-3-ylmethyl)phenoxy]-5-methylbenzonitrile', 'Onatasertib',
          "N'-(5-CHLORO-1,3-BENZODIOXOL-4-YL)-N-(3-MORPHOLIN-4-YLPHENYL)PYRIMIDINE-2,4-DIAMINE", 'Dactolisib',
          'Capivasertib', 'Cenerimod', 'PTI-428', 'V116517', 'AMD-070',
          '6-(2,4-DIAMINO-6-ETHYLPYRIMIDIN-5-YL)-4-(3-METHOXYPROPYL)-2,2-DIMETHYL-2H-1,4-BENZOXAZIN-3(4H)-ONE',
          'Vincamine', '(2Z)-1-(5-Chloro-1H-indol-3-yl)-3-hydroxy-3-(1H-tetrazol-5-yl)-2-propen-1-one', 'Sisunatovir',
          'Quinacrine mustard', 'Verdiperstat', 'Celecoxib', '3-{3-[(DIMETHYLAMINO)METHYL]-1H-INDOL-7-YL}PROPAN-1-OL',
          'Ethacridine', '(2S)-2-AMINO-4-(METHYLSULFANYL)-1-(1,3-THIAZOL-2-YL)BUTANE-1,1-DIOL',
          'N-(2,6-dimethylphenyl)-5-phenylimidazo[1,5-a]pyrazin-8-amine', 'Ezutromid', 'Nedocromil', 'Ziresovir',
          'Tavapadon', 'Pardoprunox', 'K-134', 'Pyrrole-2-Carboxylate', 'Sertaconazole', '6-Methylpurine',
          '2-amino-5-[3-(1-ethyl-1H-pyrazol-5-yl)-1H-pyrrolo[2,3-b]pyridin-5-yl]-N,N-dimethylbenzamide', 'Linopirdine',
          'Vecabrutinib', 'Venglustat',
          '4-{[1-Methyl-5-(2-Methyl-Benzoimidazol-1-Ylmethyl)-1h-Benzoimidazol-2-Ylmethyl]-Amino}-Benzamidine',
          'N,4-dimethyl-3-[(1-phenyl-1H-pyrazolo[3,4-d]pyrimidin-4-yl)amino]benzamide', 'NP-G2-044', 'Olutasidenib',
          'Simpinicline', 'PP-121', 'PHENYLAMINOIMIDAZO(1,2-ALPHA)PYRIDINE', 'GW810781', 'PSI-697',
          'N-(4-{[(3S)-3-(dimethylamino)pyrrolidin-1-yl]carbonyl}phenyl)-5-fluoro-4-[2-methyl-1-(1-methylethyl)-1H-imidazol-5-yl]pyrimidin-2-amine',
          'Hydronidone', 'Pyrazole', 'Bromazepam', 'AZD-8418',
          '4-({[4-(3-METHYLBENZOYL)PYRIDIN-2-YL]AMINO}METHYL)BENZENECARBOXIMIDAMIDE', 'Fanapanel', 'Dianicline',
          '[2-(5-Mercapto-[1,3,4]thiadiazol-2-ylcarbamoyl)-1-phenyl-ethyl]-carbamic acid benzyl ester',
          "(4R)-7,8-dichloro-1',9-dimethyl-1-oxo-1,2,4,9-tetrahydrospiro[beta-carboline-3,4'-piperidine]-4-carbonitrile",
          'E-7820', '7-(aminomethyl)-6-(2-chlorophenyl)-1-methyl-1H-benzimidazole-5-carbonitrile', 'AZD-1386',
          'Imaradenant', 'Pranoprofen', 'Tisopurine', 'Decernotinib', 'Triazolopyridine',
          '9-ACETYL-2,3,4,9-TETRAHYDRO-1H-CARBAZOL-1-ONE', 'PF-06700841', 'Rogaratinib', 'ABT-288', 'Buparlisib',
          '4-chloro-6-{5-[(2-morpholin-4-ylethyl)amino]-1,2-benzisoxazol-3-yl}benzene-1,3-diol', 'Cyanocinnoline',
          'Zamicastat',
          'Cp403700, (S)-1-{2-[(5-Chloro-1h-Indole-2-Carbonyl)-Amino]-3-Phenyl-Propionyl}-Azetidine-3-Carboxylate',
          'Isatoic anhydride', 'UK-500001', 'Droxicam', 'LY-3023414', 'Urapidil',
          '5-(7-(4-(4,5-dihydro-2-oxazolyl)phenoxy)heptyl)-3-methyl isoxazole',
          '2-CHLORO-N-[(3R)-2-OXO-1,2,3,4-TETRAHYDROQUINOLIN-3-YL]-6H-THIENO[2,3-B]PYRROLE-5-CARBOXAMIDE', 'Triapine',
          'PF-04991532', '1-[1-(3-aminophenyl)-3-tert-butyl-1H-pyrazol-5-yl]-3-phenylurea', 'Benzoylformic Acid',
          'N~2~-1,3-BENZOXAZOL-2-YL-3-CYCLOHEXYL-N-{2-[(4-METHOXYPHENYL)AMINO]ETHYL}-L-ALANINAMIDE', 'PF-05175157',
          '1,N6-Ethenoadenine', '6-[4-(2-fluorophenyl)-1,3-oxazol-5-yl]-N-(1-methylethyl)-1,3-benzothiazol-2-amine',
          'Rilematovir', '(2S)-2-{[3-(3-aminophenyl)imidazo[1,2-b]pyridazin-6-yl]amino}-3-methylbutan-1-ol',
          'Pantoprazole', 'Piclozotan', 'LGD2941', 'Ocinaplon', 'Theodrenaline', 'Pexacerfont', 'Acumapimod',
          'Pumosetrag', '2-(methylsulfanyl)-5-(thiophen-2-ylmethyl)-1H-imidazol-4-ol', 'TAK-831',
          '4-{4-[(5-hydroxy-2-methylphenyl)amino]quinolin-7-yl}-1,3-thiazole-2-carbaldehyde', 'Cadralazine',
          'GSK-1059615', '(2R)-1-{4-[(4-Anilino-5-bromo-2-pyrimidinyl)amino]phenoxy}-3-(dimethylamino)-2-propanol',
          'Harmine', 'Proquazone', 'Viminol', 'GSK2798745', 'Eletriptan', 'Imiquimod', 'Losmapimod',
          '3-{[(4-methylphenyl)sulfonyl]amino}propyl pyridin-4-ylcarbamate',
          '4-(1h-Imidazol-4-Yl)-3-(5-Ethyl-2,4-Dihydroxy-Phenyl)-1h-Pyrazole', 'ONO-2952', 'Asciminib', 'Lifirafenib',
          'Metoserpate', 'Aminophenazone', 'Pirfenidone', 'Pelitinib', 'PF-05089771', 'Moxaverine', 'Azosemide',
          '5-[3-(2-METHOXYPHENYL)-1H-PYRROLO[2,3-B]PYRIDIN-5-YL]-N,N-DIMETHYLPYRIDINE-3-CARBOXAMIDE', 'Mefloquine',
          'Sofinicline', '4-(2-(3-Propyl-(1,2,4)oxadiazol-5-yl)-vinyl)-benzene-1,2-diol', 'Capadenoson', 'Azintamide',
          'Carbimazole', 'Tenonitrozole',
          '5-(5-chloro-2,4-dihydroxyphenyl)-N-ethyl-4-[4-(morpholin-4-ylmethyl)phenyl]isoxazole-3-carboxamide',
          '3,5-Diaminophthalhydrazide', 'Pirlindole', '5,8-dimethoxy-1,4-dimethylquinolin-2(1H)-one',
          'N-[2-(6-AMINO-4-METHYLPYRIDIN-2-YL)ETHYL]-4-CYANOBENZAMIDE',
          '1-[2-(3-ACETYL-2-HYDROXY-6-METHOXY-PHENYL)-CYCLOPROPYL]-3-(5-CYANO-PYRIDIN-2-YL)-THIOUREA', 'Tipifarnib',
          'Balipodect', 'PF-00356231',
          'N-{3-[(7ar,12as,12bs)-7-Oxo-1,3,4,6,7,7a,12a,12b-Octahydroindolo[2,3-a]Quinolizin-12(2h)-Yl]Propyl}Propane-2-Sulfonamide',
          '4-METHYL-PENTANOIC ACID {1-[4-GUANIDINO-1-(THIAZOLE-2-CARBONYL)-BUTYLCARBAMOYL]-2-METHYL-PROPYL}-AMIDE',
          'Cyclo(his-pro)', 'Toreforant', 'Pruvanserin',
          '[PHENYLALANINYL-PROLINYL]-[2-(PYRIDIN-4-YLAMINO)-ETHYL]-AMINE', 'Disperse Blue 106', 'Ispronicline',
          '5-{4-[(3,5-DIFLUOROBENZYL)AMINO]PHENYL}-6-ETHYLPYRIMIDINE-2,4-DIAMINE',
          '4-{4-[4-(3-AMINOPROPOXY)PHENYL]-1H-PYRAZOL-5-YL}-6-CHLOROBENZENE-1,3-DIOL', 'Fipamezole',
          '5-Bromonicotinamide', '3-methyl-N-(pyridin-4-ylmethyl)imidazo[1,2-a]pyrazin-8-amine', '5-Methylpyrrole',
          'AZD-9496', 'Zuranolone', '5-[1-(4-methoxyphenyl)-1H-benzimidazol-6-yl]-1,3,4-oxadiazole-2(3H)-thione',
          'AZD-1236', 'Naratriptan', 'Harmaline', 'Rezatomidine', 'Florasulam',
          '9-HYDROXY-6-(3-HYDROXYPROPYL)-4-(2-METHOXYPHENYL)PYRROLO[3,4-C]CARBAZOLE-1,3(2H,6H)-DIONE',
          'Dasolampanel etibutil', 'N-[2-(5-methyl-4H-1,2,4-triazol-3-yl)phenyl]-7H-pyrrolo[2,3-d]pyrimidin-4-amine',
          'Verucerfont', 'Bunazosin',
          '(2R)-4-[(8R)-8-METHYL-2-(TRIFLUOROMETHYL)-5,6-DIHYDRO[1,2,4]TRIAZOLO[1,5-A]PYRAZIN-7(8H)-YL]-4-OXO-1-(2,4,5-TRIFLUOROPHENYL)BUTAN-2-AMINE',
          'Etrasimod', 'GSK-3117391', 'Clortermine', 'Epibatidine', 'Molidustat', 'Oxaprozin', 'Vebreltinib',
          '2-Isobutyl-3-methoxypyrazine', 'CB-103',
          '7-Methoxy-8-[1-(Methylsulfonyl)-1h-Pyrazol-4-Yl]Naphthalene-2-Carboximidamide', 'Fenpyroximate',
          'Trilaciclib', '1-(5-BROMO-PYRIDIN-2-YL)-3-[2-(6-FLUORO-2-HYDROXY-3-PROPIONYL-PHENYL)-CYCLOPROPYL]-UREA',
          'Ifidancitinib', 'AMG-319', 'Gandotinib', 'Navoximod', 'Xylose-Derived Imidazole', 'Oglemilast', 'AZD-5069',
          'Dapiprazole', 'RAD-140', 'N-[(1R)-3-(4-HYDROXYPHENYL)-1-METHYLPROPYL]-2-(2-PHENYL-1H-INDOL-3-YL)ACETAMIDE',
          'Cefatrizine',
          'N-((1R,2R)-2-(5-CHLORO-1H-INDOLE-2-CARBOXAMIDO)CYCLOHEXYL)-5-METHYL-4,5,6,7-TETRAHYDROTHIAZOLO[5,4-C]PYRIDINE-2-CARBOXAMIDE',
          'Sapanisertib', 'Aptazapine', 'PX-12', 'Vosaroxin', '4-Iodopyrazole',
          '2,4-Diamino-6-Phenyl-5,6,7,8,-Tetrahydropteridine', 'N,N-DIMETHYL(5-(PYRIDIN-3-YL)FURAN-2-YL)METHANAMINE',
          "N-{(3R,4S)-4-[(6-amino-4-methylpyridin-2-yl)methyl]pyrrolidin-3-yl}-N'-(3-chlorobenzyl)ethane-1,2-diamine",
          'D-157495', 'N-{2-methyl-5-[(6-phenylpyrimidin-4-yl)amino]phenyl}methanesulfonamide', 'Ethionamide',
          'GSK-2881078', 'PF-06282999', 'Vorolanib', 'BMS-919373'
          ]

val3 = ['Elvitegravir', 'Glumetinib', 'Fenamole', 'Lonidamine', 'CXD101', 'GSK-2982772', 'Olomoucine', 'R-82913',
        'Etodolac', 'Zonisamide', 'Amlexanox', 'PF-05212377', 'Arotinolol', 'Thiacloprid',
        'N-(2-Aminoethyl)-5-Chloroisoquinoline-8-Sulfonamide', 'Thiohexam', 'LY-2811376', 'BOS172722',
        'N-[2-(METHYLAMINO)ETHYL]-5-ISOQUINOLINESULFONAMIDE', 'Samuraciclib', 'Amuvatinib', 'Pecavaptan',
        'Bradanicline', '2-[(2-methoxy-5-methylphenoxy)methyl]pyridine', 'AZD-1940', 'Belotecan', 'Cenobamate',
        '2-(1H-pyrrol-1-ylcarbonyl)benzene-1,3,5-triol', 'Latrepirdine',
        'Carboxymethylthio-3-(3-Chlorophenyl)-1,2,4-Oxadiazol', 'Terevalefim', 'Chlormidazole', 'Rupatadine',
        'Tolmetin', 'Lanifibranor', 'N1-CYCLOPENTYL-N2-(THIAZOL-2-YL)OXALAMIDE', 'Benzimidazole', 'GSK-239512',
        'PF-06260414', 'Etomidate', 'Azatadine', 'RO-5028442', 'Abiraterone', 'arazepide', 'Upadacitinib', '2X-121',
        'CRA_10972', 'ABX-464', 'Tegoprazan'
        ]

test3 = ['Tizanidine', '4-{[5-chloro-4-(1H-indol-3-yl)pyrimidin-2-yl]amino}-N-ethylpiperidine-1-carboxamide',
         'Anastrozole', 'Uracil', 'Nedisertib', 'Cilostazol',
         '1-(3,5-DICHLOROPHENYL)-5-METHYL-1H-1,2,4-TRIAZOLE-3-CARBOXYLIC ACID', 'Fluconazole', 'Ipatasertib',
         'N-acetylhistamine', 'Tinostamustine', 'Setipiprant', 'R-1487', 'GDC-0134', 'Tolimidone',
         '(2E)-N-hydroxy-3-[1-methyl-4-(phenylacetyl)-1H-pyrrol-2-yl]prop-2-enamide',
         '3-[(2,2-DIMETHYLPROPANOYL)AMINO]-N-1,3-THIAZOL-2-YLPYRIDINE-2-CARBOXAMIDE', 'Pozanicline', 'Veliparib',
         'Ondelopran', 'IQP-0528', 'Ranirestat', '7-Nitroindazole', 'Phosphonoacetohydroxamic Acid', 'Perampanel',
         'N-phenyl-1H-pyrrolo[2,3-b]pyridin-3-amine', 'Tazobactam',
         'N~4~-methyl-N~4~-(3-methyl-1H-indazol-6-yl)-N~2~-(3,4,5-trimethoxyphenyl)pyrimidine-2,4-diamine',
         'Mercaptocarboxylate Inhibitor', 'LY-3200882', 'Flosequinan', 'PF-04691502', 'GSK-2018682', 'Cefazedone',
         '{3-[(5-CHLORO-1,3-BENZOTHIAZOL-2-YL)METHYL]-2,4-DIOXO-3,4-DIHYDROPYRIMIDIN-1(2H)-YL}ACETIC ACID',
         'Selgantolimod',
         '2,2,2-TRIFLUORO-1-{5-[(3-PHENYL-5,6-DIHYDROIMIDAZO[1,2-A]PYRAZIN-7(8H)-YL)CARBONYL]THIOPHEN-2-YL}ETHANE-1,1-DIOL',
         'Isradipine', 'Gluco-Phenylimidazole',
         "TRW3-(2-AMINO-3-HYDROXY-PROPYL)-6-(N'-CYCLOHEXYL-HYDRAZINO)OCTAHYDRO-INDOL-7-OL", 'Tucidinostat',
         'Pemigatinib', '5-CHLORO-6-METHYL-N-(2-PHENYLETHYL)-2-PYRIDIN-2-YLPYRIMIDIN-4-AMINE',
         '[1-(4-Fluorobenzyl)Cyclobutyl]Methyl (1s)-1-[Oxo(1h-Pyrazol-5-Ylamino)Acetyl]Pentylcarbamate',
         '(1S,3R,6S)-4-oxo-6-{4-[(2-phenylquinolin-4-yl)methoxy]phenyl}-5-azaspiro[2.4]heptane-1-carboxylic acid',
         '1-[(2S)-4-(5-BROMO-1H-PYRAZOLO[3,4-B]PYRIDIN-4-YL)MORPHOLIN-2-YL]METHANAMINE',
         'ethyl 3-[(E)-2-amino-1-cyanoethenyl]-6,7-dichloro-1-methyl-1H-indole-2-carboxylate',
         '4-(5,11-DIOXO-5H-INDENO[1,2-C]ISOQUINOLIN-6(11H)-YL)BUTANOATE'
         ]

train3 = ['2-Methoxy-3-Isopropylpyrazine', '3-(1-NAPHTHYLMETHOXY)PYRIDIN-2-AMINE',
          'Di-2-pyridylketone 4-cyclohexyl-4-methyl-3-thiosemicarbazone', 'Methimazole', 'S-777469', 'ORE-1001',
          'Oliceridine', 'N-6022', 'GDC-0425',
          "4-Bromo-3-(5'-Carboxy-4'-Chloro-2'-Fluorophenyl)-1-Methyl-5-Trifluoromethyl-Pyrazol", 'Tivirapine',
          'Rufinamide', '6-BENZYL-1-BENZYLOXYMETHYL-5-ISOPROPYL URACIL', 'Avadomide',
          '1-[(4S)-4-amino-5-(1,3-benzothiazol-2-yl)-5-oxopentyl]guanidine', 'Cipargamin',
          "1-(6-CYANO-3-PYRIDYLCARBONYL)-5',8'-DIFLUOROSPIRO[PIPERIDINE-4,2'(1'H)-QUINAZOLINE]-4'-AMINE",
          'JNJ-39393406', 'Piribedil', 'Imidazole', 'Pipequaline', 'Lavoltidine',
          'N-{3-[5-(1H-1,2,4-triazol-3-yl)-1H-indazol-3-yl]phenyl}furan-2-carboxamide', 'Tinengotinib', 'Uracil C-13',
          'Minodronic acid', 'Tovorafenib', 'PF-03463275', '3-(4-CHLOROPHENYL)-5-(METHYLTHIO)-4H-1,2,4-TRIAZOLE',
          'AZD-1981', 'Revaprazan', 'Rabeximod', '3-Methyladenine', '5-Aminoisoquinoline',
          '6-HYDROXY-1,3-BENZOTHIAZOLE-2-SULFONAMIDE', '1,2,4-Triazole',
          '1-(2-Chlorophenyl)-3,5-Dimethyl-1h-Pyrazole-4-Carboxylic Acid Ethyl Ester',
          '6-CHLORO-9-HYDROXY-1,3-DIMETHYL-1,9-DIHYDRO-4H-PYRAZOLO[3,4-B]QUINOLIN-4-ONE', 'Dilmapimod', 'CGS-27023',
          'N-[2-methyl-5-(methylcarbamoyl)phenyl]-2-{[(1R)-1-methylpropyl]amino}-1,3-thiazole-5-carboxamide',
          'Nitenpyram', 'AMG-900', '6-CHLORO-4-(CYCLOHEXYLOXY)-3-PROPYLQUINOLIN-2(1H)-ONE', 'Roflumilast', 'Posizolid',
          'Carboxyamidotriazole', '5-(2-chlorophenyl)-1,3,4-thiadiazole-2-sulfonamide',
          '4-CHLORO-6-(4-PIPERAZIN-1-YL-1H-PYRAZOL-5-YL)BENZENE-1,3-DIOL',
          '(2R)-1-[(5,6-DIPHENYL-7H-PYRROLO[2,3-D]PYRIMIDIN-4-YL)AMINO]PROPAN-2-OL', 'Birabresib', 'Cimicoxib',
          'Snubh-nm-333 F-18', 'Polaprezinc', 'Pradigastat',
          "(4R)-7-chloro-9-methyl-1-oxo-1,2,4,9-tetrahydrospiro[beta-carboline-3,4'-piperidine]-4-carbonitrile",
          'Paltusotine', 'Diclazuril', 'Butalamine', 'Riluzole', 'N-METHYL-1-[4-(9H-PURIN-6-YL)PHENYL]METHANAMINE',
          'Aleplasinin', 'Linzagolix', 'Gimeracil', 'BMS-488043', 'Clonazolam', 'Ciclopirox', 'Naphthoquine',
          '4-[(3-BROMO-4-O-SULFAMOYLBENZYL)(4-CYANOPHENYL)AMINO]-4H-[1,2,4]-TRIAZOLE',
          '5-[(Z)-(5-Chloro-2-oxo-1,2-dihydro-3H-indol-3-ylidene)methyl]-N,2,4-trimethyl-1H-pyrrole-3-carboxamide',
          'AG-24322', 'Pamiparib', 'PF-05241328', 'Dalfampridine', '2-Pyridinethiol', 'Mivavotinib',
          '2-(4-ETHYLPIPERAZIN-1-YL)-4-(PHENYLAMINO)PYRAZOLO[1,5-A][1,3,5]TRIAZINE-8-CARBONITRILE', 'Metralindole',
          'Lamotrigine', 'Cinoxacin', 'Anagliptin', 'Losartan', 'Allopurinol', 'Frovatriptan', 'Deucravacitinib',
          '2-(2-chloropyridin-4-yl)-4-methyl-1H-isoindole-1,3(2H)-dione', 'Irbesartan',
          '3-[2-bromo-4-(1H-pyrazolo[3,4-c]pyridazin-3-ylmethyl)phenoxy]-5-methylbenzonitrile', 'Onatasertib',
          "N'-(5-CHLORO-1,3-BENZODIOXOL-4-YL)-N-(3-MORPHOLIN-4-YLPHENYL)PYRIMIDINE-2,4-DIAMINE", 'Dactolisib',
          'Capivasertib', 'Cenerimod', 'PTI-428', 'V116517', 'AMD-070',
          '6-(2,4-DIAMINO-6-ETHYLPYRIMIDIN-5-YL)-4-(3-METHOXYPROPYL)-2,2-DIMETHYL-2H-1,4-BENZOXAZIN-3(4H)-ONE',
          'Vincamine', '(2Z)-1-(5-Chloro-1H-indol-3-yl)-3-hydroxy-3-(1H-tetrazol-5-yl)-2-propen-1-one', 'Sisunatovir',
          'Quinacrine mustard', 'Verdiperstat', 'Celecoxib', '3-{3-[(DIMETHYLAMINO)METHYL]-1H-INDOL-7-YL}PROPAN-1-OL',
          'Ethacridine', '(2S)-2-AMINO-4-(METHYLSULFANYL)-1-(1,3-THIAZOL-2-YL)BUTANE-1,1-DIOL',
          'N-(2,6-dimethylphenyl)-5-phenylimidazo[1,5-a]pyrazin-8-amine', 'Ezutromid', 'Nedocromil', 'Ziresovir',
          'Tavapadon', 'Pardoprunox', 'K-134', 'Pyrrole-2-Carboxylate', 'Sertaconazole', '6-Methylpurine',
          '2-amino-5-[3-(1-ethyl-1H-pyrazol-5-yl)-1H-pyrrolo[2,3-b]pyridin-5-yl]-N,N-dimethylbenzamide', 'Linopirdine',
          'Vecabrutinib', 'Venglustat',
          '4-{[1-Methyl-5-(2-Methyl-Benzoimidazol-1-Ylmethyl)-1h-Benzoimidazol-2-Ylmethyl]-Amino}-Benzamidine',
          'N,4-dimethyl-3-[(1-phenyl-1H-pyrazolo[3,4-d]pyrimidin-4-yl)amino]benzamide', 'NP-G2-044', 'Olutasidenib',
          'Simpinicline', 'PP-121', 'PHENYLAMINOIMIDAZO(1,2-ALPHA)PYRIDINE', 'GW810781', 'PSI-697',
          'N-(4-{[(3S)-3-(dimethylamino)pyrrolidin-1-yl]carbonyl}phenyl)-5-fluoro-4-[2-methyl-1-(1-methylethyl)-1H-imidazol-5-yl]pyrimidin-2-amine',
          'Hydronidone', 'Pyrazole', 'Bromazepam', 'AZD-8418',
          '4-({[4-(3-METHYLBENZOYL)PYRIDIN-2-YL]AMINO}METHYL)BENZENECARBOXIMIDAMIDE', 'Fanapanel', 'Dianicline',
          '[2-(5-Mercapto-[1,3,4]thiadiazol-2-ylcarbamoyl)-1-phenyl-ethyl]-carbamic acid benzyl ester',
          "(4R)-7,8-dichloro-1',9-dimethyl-1-oxo-1,2,4,9-tetrahydrospiro[beta-carboline-3,4'-piperidine]-4-carbonitrile",
          'E-7820', '7-(aminomethyl)-6-(2-chlorophenyl)-1-methyl-1H-benzimidazole-5-carbonitrile', 'AZD-1386',
          'Imaradenant', 'Pranoprofen', 'Tisopurine', 'Decernotinib', 'Triazolopyridine',
          '9-ACETYL-2,3,4,9-TETRAHYDRO-1H-CARBAZOL-1-ONE', 'PF-06700841', 'Rogaratinib', 'ABT-288', 'Buparlisib',
          '4-chloro-6-{5-[(2-morpholin-4-ylethyl)amino]-1,2-benzisoxazol-3-yl}benzene-1,3-diol', 'Cyanocinnoline',
          'Zamicastat',
          'Cp403700, (S)-1-{2-[(5-Chloro-1h-Indole-2-Carbonyl)-Amino]-3-Phenyl-Propionyl}-Azetidine-3-Carboxylate',
          'Isatoic anhydride', 'UK-500001', 'Droxicam', 'LY-3023414', 'Urapidil',
          '5-(7-(4-(4,5-dihydro-2-oxazolyl)phenoxy)heptyl)-3-methyl isoxazole',
          '2-CHLORO-N-[(3R)-2-OXO-1,2,3,4-TETRAHYDROQUINOLIN-3-YL]-6H-THIENO[2,3-B]PYRROLE-5-CARBOXAMIDE', 'Triapine',
          'PF-04991532', '1-[1-(3-aminophenyl)-3-tert-butyl-1H-pyrazol-5-yl]-3-phenylurea', 'Benzoylformic Acid',
          'N~2~-1,3-BENZOXAZOL-2-YL-3-CYCLOHEXYL-N-{2-[(4-METHOXYPHENYL)AMINO]ETHYL}-L-ALANINAMIDE', 'PF-05175157',
          '1,N6-Ethenoadenine', '6-[4-(2-fluorophenyl)-1,3-oxazol-5-yl]-N-(1-methylethyl)-1,3-benzothiazol-2-amine',
          'Rilematovir', '(2S)-2-{[3-(3-aminophenyl)imidazo[1,2-b]pyridazin-6-yl]amino}-3-methylbutan-1-ol',
          'Pantoprazole', 'Piclozotan', 'LGD2941', 'Ocinaplon', 'Theodrenaline', 'Pexacerfont', 'Acumapimod',
          'Pumosetrag', '2-(methylsulfanyl)-5-(thiophen-2-ylmethyl)-1H-imidazol-4-ol', 'TAK-831',
          '4-{4-[(5-hydroxy-2-methylphenyl)amino]quinolin-7-yl}-1,3-thiazole-2-carbaldehyde', 'Cadralazine',
          'GSK-1059615', '(2R)-1-{4-[(4-Anilino-5-bromo-2-pyrimidinyl)amino]phenoxy}-3-(dimethylamino)-2-propanol',
          'Harmine', 'Proquazone', 'Viminol', 'GSK2798745', 'Eletriptan', 'Imiquimod', 'Losmapimod',
          '3-{[(4-methylphenyl)sulfonyl]amino}propyl pyridin-4-ylcarbamate',
          '4-(1h-Imidazol-4-Yl)-3-(5-Ethyl-2,4-Dihydroxy-Phenyl)-1h-Pyrazole', 'ONO-2952', 'Asciminib', 'Lifirafenib',
          'Metoserpate', 'Aminophenazone', 'Pirfenidone', 'Pelitinib', 'PF-05089771', 'Moxaverine', 'Azosemide',
          '5-[3-(2-METHOXYPHENYL)-1H-PYRROLO[2,3-B]PYRIDIN-5-YL]-N,N-DIMETHYLPYRIDINE-3-CARBOXAMIDE', 'Mefloquine',
          'Sofinicline', '4-(2-(3-Propyl-(1,2,4)oxadiazol-5-yl)-vinyl)-benzene-1,2-diol', 'Capadenoson', 'Azintamide',
          'Carbimazole', 'Tenonitrozole',
          '5-(5-chloro-2,4-dihydroxyphenyl)-N-ethyl-4-[4-(morpholin-4-ylmethyl)phenyl]isoxazole-3-carboxamide',
          '3,5-Diaminophthalhydrazide', 'Pirlindole', '5,8-dimethoxy-1,4-dimethylquinolin-2(1H)-one',
          'N-[2-(6-AMINO-4-METHYLPYRIDIN-2-YL)ETHYL]-4-CYANOBENZAMIDE',
          '1-[2-(3-ACETYL-2-HYDROXY-6-METHOXY-PHENYL)-CYCLOPROPYL]-3-(5-CYANO-PYRIDIN-2-YL)-THIOUREA', 'Tipifarnib',
          'Balipodect', 'PF-00356231',
          'N-{3-[(7ar,12as,12bs)-7-Oxo-1,3,4,6,7,7a,12a,12b-Octahydroindolo[2,3-a]Quinolizin-12(2h)-Yl]Propyl}Propane-2-Sulfonamide',
          '4-METHYL-PENTANOIC ACID {1-[4-GUANIDINO-1-(THIAZOLE-2-CARBONYL)-BUTYLCARBAMOYL]-2-METHYL-PROPYL}-AMIDE',
          'Cyclo(his-pro)', 'Toreforant', 'Pruvanserin',
          '[PHENYLALANINYL-PROLINYL]-[2-(PYRIDIN-4-YLAMINO)-ETHYL]-AMINE', 'Disperse Blue 106', 'Ispronicline',
          '5-{4-[(3,5-DIFLUOROBENZYL)AMINO]PHENYL}-6-ETHYLPYRIMIDINE-2,4-DIAMINE',
          '4-{4-[4-(3-AMINOPROPOXY)PHENYL]-1H-PYRAZOL-5-YL}-6-CHLOROBENZENE-1,3-DIOL', 'Fipamezole',
          '5-Bromonicotinamide', '3-methyl-N-(pyridin-4-ylmethyl)imidazo[1,2-a]pyrazin-8-amine', '5-Methylpyrrole',
          'AZD-9496', 'Zuranolone', '5-[1-(4-methoxyphenyl)-1H-benzimidazol-6-yl]-1,3,4-oxadiazole-2(3H)-thione',
          'AZD-1236', 'Naratriptan', 'Harmaline', 'Rezatomidine', 'Florasulam',
          '9-HYDROXY-6-(3-HYDROXYPROPYL)-4-(2-METHOXYPHENYL)PYRROLO[3,4-C]CARBAZOLE-1,3(2H,6H)-DIONE',
          'Dasolampanel etibutil', 'N-[2-(5-methyl-4H-1,2,4-triazol-3-yl)phenyl]-7H-pyrrolo[2,3-d]pyrimidin-4-amine',
          'Verucerfont', 'Bunazosin',
          '(2R)-4-[(8R)-8-METHYL-2-(TRIFLUOROMETHYL)-5,6-DIHYDRO[1,2,4]TRIAZOLO[1,5-A]PYRAZIN-7(8H)-YL]-4-OXO-1-(2,4,5-TRIFLUOROPHENYL)BUTAN-2-AMINE',
          'Etrasimod', 'GSK-3117391', 'Clortermine', 'Epibatidine', 'Molidustat', 'Oxaprozin', 'Vebreltinib',
          '2-Isobutyl-3-methoxypyrazine', 'CB-103',
          '7-Methoxy-8-[1-(Methylsulfonyl)-1h-Pyrazol-4-Yl]Naphthalene-2-Carboximidamide', 'Fenpyroximate',
          'Trilaciclib', '1-(5-BROMO-PYRIDIN-2-YL)-3-[2-(6-FLUORO-2-HYDROXY-3-PROPIONYL-PHENYL)-CYCLOPROPYL]-UREA',
          'Ifidancitinib', 'AMG-319', 'Gandotinib', 'Navoximod', 'Xylose-Derived Imidazole', 'Oglemilast', 'AZD-5069',
          'Dapiprazole', 'RAD-140', 'N-[(1R)-3-(4-HYDROXYPHENYL)-1-METHYLPROPYL]-2-(2-PHENYL-1H-INDOL-3-YL)ACETAMIDE',
          'Cefatrizine',
          'N-((1R,2R)-2-(5-CHLORO-1H-INDOLE-2-CARBOXAMIDO)CYCLOHEXYL)-5-METHYL-4,5,6,7-TETRAHYDROTHIAZOLO[5,4-C]PYRIDINE-2-CARBOXAMIDE',
          'Sapanisertib', 'Aptazapine', 'PX-12', 'Vosaroxin', '4-Iodopyrazole',
          '2,4-Diamino-6-Phenyl-5,6,7,8,-Tetrahydropteridine', 'N,N-DIMETHYL(5-(PYRIDIN-3-YL)FURAN-2-YL)METHANAMINE',
          "N-{(3R,4S)-4-[(6-amino-4-methylpyridin-2-yl)methyl]pyrrolidin-3-yl}-N'-(3-chlorobenzyl)ethane-1,2-diamine",
          'D-157495', 'N-{2-methyl-5-[(6-phenylpyrimidin-4-yl)amino]phenyl}methanesulfonamide', 'Ethionamide',
          'GSK-2881078', 'PF-06282999', 'Vorolanib', 'BMS-919373', 'BTRX-246040', 'Bendamustine', 'Letrazuril',
          'HYDROXY[3-(6-METHYLPYRIDIN-2-YL)PROPYL]FORMAMIDE',
          '4-(2-(1H-IMIDAZOL-4-YL)ETHYLAMINO)-2-(PHENYLAMINO)PYRAZOLO[1,5-A][1,3,5]TRIAZINE-8-CARBONITRILE',
          'Nelotanserin', 'Cicletanine',
          '5-CHLORO-THIOPHENE-2-CARBOXYLIC ACID ((3S,4S)-4-FLUORO- 1-{[2-FLUORO-4-(2-OXO-2H-PYRIDIN-1-YL)-PHENYLCARBAMOYL]-METHYL}-PYRROLIDIN-3-YL)-AMIDE',
          'Zilpaterol', 'Duvelisib', 'F-15599', 'Pindolol', 'Reldesemtiv', 'Methylisothiazolinone', 'MK-8245',
          'Brensocatib', '2-Amino-6-Chloropyrazine', '(2-AMINO-1,3-OXAZOL-5-YL)-(3-BROMOPHENYL)METHANONE', 'Foslinanib',
          'GSK-356278', 'Lanabecestat', '9-Deazahypoxanthine', 'INCB-057643', 'Surinabant', 'PXT 3003',
          '3-(6-HYDROXY-NAPHTHALEN-2-YL)-BENZO[D]ISOOXAZOL-6-OL', 'Oltipraz', 'ONO-8539', 'Chlorzoxazone',
          'Deferiprone', 'SB-705498', 'Nidufexor',
          '7-amino-2-tert-butyl-4-{[2-(1H-imidazol-4-yl)ethyl]amino}pyrido[2,3-d]pyrimidine-6-carboxamide',
          'Imiglitazar', 'Rimacalib', 'Tetrazolyl Histidine', 'JNJ-54175446', 'Foliglurax', 'CH-5132799',
          '4-{2-[(7-amino-2-furan-2-yl[1,2,4]triazolo[1,5-a][1,3,5]triazin-5-yl)amino]ethyl}phenol',
          'N-[(2R)-2-{[(2S)-2-(1,3-benzoxazol-2-yl)pyrrolidin-1-yl]carbonyl}hexyl]-N-hydroxyformamide',
          'Acid yellow 54 free acid', '3-phenyl-5-(1H-pyrazol-3-yl)isoxazole', 'MK-5108', 'GLPG-1205', 'Pioglitazone',
          '4-[1-allyl-7-(trifluoromethyl)-1H-indazol-3-yl]benzene-1,3-diol', 'MK-212', 'Bamifylline',
          '2-(Sec-Butyl)Thiazole', 'Rivanicline', '5-benzyl-1,3-thiazol-2-amine', 'Sitamaquine',
          '(2R,3R)-N^1^-[(1S)-2,2-DIMETHYL-1-(METHYLCARBAMOYL)PROPYL]-N^4^-HYDROXY-2-(2-METHYLPROPYL)-3-{[(1,3-THIAZOL-2-YLCARBONYL)AMINO]METHYL}BUTANEDIAMIDE',
          'Cinchocaine', 'Simurosertib', 'BTRX-335140', 'Afuresertib', 'Enitociclib', 'CCX-140',
          "N-{(3S,4S)-4-[(6-AMINO-4-METHYLPYRIDIN-2-YL)METHYL]PYRROLIDIN-3-YL}-N'-(4-CHLOROBENZYL)ETHANE-1,2-DIAMINE",
          '(7as,12ar,12bs)-1,2,3,4,7a,12,12a,12b-Octahydroindolo[2,3-a]Quinolizin-7(6h)-One', 'CAN-508', 'Cinalukast',
          '(1s,2s)-1-Amino-1-(1,3-Thiazol-2-Yl)Propan-2-Ol', 'Zolpidem', 'Galeterone', 'Mercaptopurine', 'GZ-389988',
          'Varespladib methyl', 'Farampator', 'PD-173952', 'Pirbuterol', 'Amsacrine', 'AZD-0328', 'Umifenovir',
          '(S)-wiskostatin', '4-(6-{[(1R)-1-(hydroxymethyl)propyl]amino}imidazo[1,2-b]pyridazin-3-yl)benzoic acid',
          'CC-115', 'Protionamide', 'Linrodostat', 'LGH-447', 'Henatinib', 'Alprazolam', 'MK-886', 'MK-6186',
          '6-[2-(1H-INDOL-6-YL)ETHYL]PYRIDIN-2-AMINE',
          '3-(4-{2-[2-(2-bromo-acetylamino)-ethyldisulfanyl]-ethylcarbamoyl}-cyclohexylcarbamoyl)-pyrazine-2-carboxylic acid',
          '6-MORPHOLIN-4-YL-9H-PURINE',
          'N-hydroxy-5-[(3-phenyl-5,6-dihydroimidazo[1,2-a]pyrazin-7(8H)-yl)carbonyl]thiophene-2-carboxamide',
          '5-(Aminomethyl)-6-(2,4-Dichlorophenyl)-2-(3,5-Dimethoxyphenyl)Pyrimidin-4-Amine', 'Tulrampator', 'JAB-3068',
          'Ganaplacide', 'Prenoxdiazine',
          '1-[(2S)-4-(5-phenyl-1H-pyrazolo[3,4-b]pyridin-4-yl)morpholin-2-yl]methanamine', 'Trapidil',
          '4-(2-Thienyl)-1-(4-Methylbenzyl)-1h-Imidazole'
          ]

val4 = ['2-Methoxy-3-Isopropylpyrazine', '3-(1-NAPHTHYLMETHOXY)PYRIDIN-2-AMINE',
        'Di-2-pyridylketone 4-cyclohexyl-4-methyl-3-thiosemicarbazone', 'Methimazole', 'S-777469', 'ORE-1001',
        'Oliceridine', 'N-6022', 'GDC-0425',
        "4-Bromo-3-(5'-Carboxy-4'-Chloro-2'-Fluorophenyl)-1-Methyl-5-Trifluoromethyl-Pyrazol", 'Tivirapine',
        'Rufinamide', '6-BENZYL-1-BENZYLOXYMETHYL-5-ISOPROPYL URACIL', 'Avadomide',
        '1-[(4S)-4-amino-5-(1,3-benzothiazol-2-yl)-5-oxopentyl]guanidine', 'Cipargamin',
        "1-(6-CYANO-3-PYRIDYLCARBONYL)-5',8'-DIFLUOROSPIRO[PIPERIDINE-4,2'(1'H)-QUINAZOLINE]-4'-AMINE", 'JNJ-39393406',
        'Piribedil', 'Imidazole', 'Pipequaline', 'Lavoltidine',
        'N-{3-[5-(1H-1,2,4-triazol-3-yl)-1H-indazol-3-yl]phenyl}furan-2-carboxamide', 'Tinengotinib', 'Uracil C-13',
        'Minodronic acid', 'Tovorafenib', 'PF-03463275', '3-(4-CHLOROPHENYL)-5-(METHYLTHIO)-4H-1,2,4-TRIAZOLE',
        'AZD-1981', 'Revaprazan', 'Rabeximod', '3-Methyladenine', '5-Aminoisoquinoline',
        '6-HYDROXY-1,3-BENZOTHIAZOLE-2-SULFONAMIDE', '1,2,4-Triazole',
        '1-(2-Chlorophenyl)-3,5-Dimethyl-1h-Pyrazole-4-Carboxylic Acid Ethyl Ester',
        '6-CHLORO-9-HYDROXY-1,3-DIMETHYL-1,9-DIHYDRO-4H-PYRAZOLO[3,4-B]QUINOLIN-4-ONE', 'Dilmapimod', 'CGS-27023',
        'N-[2-methyl-5-(methylcarbamoyl)phenyl]-2-{[(1R)-1-methylpropyl]amino}-1,3-thiazole-5-carboxamide',
        'Nitenpyram', 'AMG-900', '6-CHLORO-4-(CYCLOHEXYLOXY)-3-PROPYLQUINOLIN-2(1H)-ONE', 'Roflumilast', 'Posizolid',
        'Carboxyamidotriazole', '5-(2-chlorophenyl)-1,3,4-thiadiazole-2-sulfonamide',
        '4-CHLORO-6-(4-PIPERAZIN-1-YL-1H-PYRAZOL-5-YL)BENZENE-1,3-DIOL'
        ]

test4 = ['(2R)-1-[(5,6-DIPHENYL-7H-PYRROLO[2,3-D]PYRIMIDIN-4-YL)AMINO]PROPAN-2-OL', 'Birabresib', 'Cimicoxib',
         'Snubh-nm-333 F-18', 'Polaprezinc', 'Pradigastat',
         "(4R)-7-chloro-9-methyl-1-oxo-1,2,4,9-tetrahydrospiro[beta-carboline-3,4'-piperidine]-4-carbonitrile",
         'Paltusotine', 'Diclazuril', 'Butalamine', 'Riluzole', 'N-METHYL-1-[4-(9H-PURIN-6-YL)PHENYL]METHANAMINE',
         'Aleplasinin', 'Linzagolix', 'Gimeracil', 'BMS-488043', 'Clonazolam', 'Ciclopirox', 'Naphthoquine',
         '4-[(3-BROMO-4-O-SULFAMOYLBENZYL)(4-CYANOPHENYL)AMINO]-4H-[1,2,4]-TRIAZOLE',
         '5-[(Z)-(5-Chloro-2-oxo-1,2-dihydro-3H-indol-3-ylidene)methyl]-N,2,4-trimethyl-1H-pyrrole-3-carboxamide',
         'AG-24322', 'Pamiparib', 'PF-05241328', 'Dalfampridine', '2-Pyridinethiol', 'Mivavotinib',
         '2-(4-ETHYLPIPERAZIN-1-YL)-4-(PHENYLAMINO)PYRAZOLO[1,5-A][1,3,5]TRIAZINE-8-CARBONITRILE', 'Metralindole',
         'Lamotrigine', 'Cinoxacin', 'Anagliptin', 'Losartan', 'Allopurinol', 'Frovatriptan', 'Deucravacitinib',
         '2-(2-chloropyridin-4-yl)-4-methyl-1H-isoindole-1,3(2H)-dione', 'Irbesartan',
         '3-[2-bromo-4-(1H-pyrazolo[3,4-c]pyridazin-3-ylmethyl)phenoxy]-5-methylbenzonitrile', 'Onatasertib',
         "N'-(5-CHLORO-1,3-BENZODIOXOL-4-YL)-N-(3-MORPHOLIN-4-YLPHENYL)PYRIMIDINE-2,4-DIAMINE", 'Dactolisib',
         'Capivasertib', 'Cenerimod', 'PTI-428', 'V116517', 'AMD-070',
         '6-(2,4-DIAMINO-6-ETHYLPYRIMIDIN-5-YL)-4-(3-METHOXYPROPYL)-2,2-DIMETHYL-2H-1,4-BENZOXAZIN-3(4H)-ONE'
         ]

train4 = ['Vincamine', '(2Z)-1-(5-Chloro-1H-indol-3-yl)-3-hydroxy-3-(1H-tetrazol-5-yl)-2-propen-1-one', 'Sisunatovir',
          'Quinacrine mustard', 'Verdiperstat', 'Celecoxib', '3-{3-[(DIMETHYLAMINO)METHYL]-1H-INDOL-7-YL}PROPAN-1-OL',
          'Ethacridine', '(2S)-2-AMINO-4-(METHYLSULFANYL)-1-(1,3-THIAZOL-2-YL)BUTANE-1,1-DIOL',
          'N-(2,6-dimethylphenyl)-5-phenylimidazo[1,5-a]pyrazin-8-amine', 'Ezutromid', 'Nedocromil', 'Ziresovir',
          'Tavapadon', 'Pardoprunox', 'K-134', 'Pyrrole-2-Carboxylate', 'Sertaconazole', '6-Methylpurine',
          '2-amino-5-[3-(1-ethyl-1H-pyrazol-5-yl)-1H-pyrrolo[2,3-b]pyridin-5-yl]-N,N-dimethylbenzamide', 'Linopirdine',
          'Vecabrutinib', 'Venglustat',
          '4-{[1-Methyl-5-(2-Methyl-Benzoimidazol-1-Ylmethyl)-1h-Benzoimidazol-2-Ylmethyl]-Amino}-Benzamidine',
          'N,4-dimethyl-3-[(1-phenyl-1H-pyrazolo[3,4-d]pyrimidin-4-yl)amino]benzamide', 'NP-G2-044', 'Olutasidenib',
          'Simpinicline', 'PP-121', 'PHENYLAMINOIMIDAZO(1,2-ALPHA)PYRIDINE', 'GW810781', 'PSI-697',
          'N-(4-{[(3S)-3-(dimethylamino)pyrrolidin-1-yl]carbonyl}phenyl)-5-fluoro-4-[2-methyl-1-(1-methylethyl)-1H-imidazol-5-yl]pyrimidin-2-amine',
          'Hydronidone', 'Pyrazole', 'Bromazepam', 'AZD-8418',
          '4-({[4-(3-METHYLBENZOYL)PYRIDIN-2-YL]AMINO}METHYL)BENZENECARBOXIMIDAMIDE', 'Fanapanel', 'Dianicline',
          '[2-(5-Mercapto-[1,3,4]thiadiazol-2-ylcarbamoyl)-1-phenyl-ethyl]-carbamic acid benzyl ester',
          "(4R)-7,8-dichloro-1',9-dimethyl-1-oxo-1,2,4,9-tetrahydrospiro[beta-carboline-3,4'-piperidine]-4-carbonitrile",
          'E-7820', '7-(aminomethyl)-6-(2-chlorophenyl)-1-methyl-1H-benzimidazole-5-carbonitrile', 'AZD-1386',
          'Imaradenant', 'Pranoprofen', 'Tisopurine', 'Decernotinib', 'Triazolopyridine',
          '9-ACETYL-2,3,4,9-TETRAHYDRO-1H-CARBAZOL-1-ONE', 'PF-06700841', 'Rogaratinib', 'ABT-288', 'Buparlisib',
          '4-chloro-6-{5-[(2-morpholin-4-ylethyl)amino]-1,2-benzisoxazol-3-yl}benzene-1,3-diol', 'Cyanocinnoline',
          'Zamicastat',
          'Cp403700, (S)-1-{2-[(5-Chloro-1h-Indole-2-Carbonyl)-Amino]-3-Phenyl-Propionyl}-Azetidine-3-Carboxylate',
          'Isatoic anhydride', 'UK-500001', 'Droxicam', 'LY-3023414', 'Urapidil',
          '5-(7-(4-(4,5-dihydro-2-oxazolyl)phenoxy)heptyl)-3-methyl isoxazole',
          '2-CHLORO-N-[(3R)-2-OXO-1,2,3,4-TETRAHYDROQUINOLIN-3-YL]-6H-THIENO[2,3-B]PYRROLE-5-CARBOXAMIDE', 'Triapine',
          'PF-04991532', '1-[1-(3-aminophenyl)-3-tert-butyl-1H-pyrazol-5-yl]-3-phenylurea', 'Benzoylformic Acid',
          'N~2~-1,3-BENZOXAZOL-2-YL-3-CYCLOHEXYL-N-{2-[(4-METHOXYPHENYL)AMINO]ETHYL}-L-ALANINAMIDE', 'PF-05175157',
          '1,N6-Ethenoadenine', '6-[4-(2-fluorophenyl)-1,3-oxazol-5-yl]-N-(1-methylethyl)-1,3-benzothiazol-2-amine',
          'Rilematovir', '(2S)-2-{[3-(3-aminophenyl)imidazo[1,2-b]pyridazin-6-yl]amino}-3-methylbutan-1-ol',
          'Pantoprazole', 'Piclozotan', 'LGD2941', 'Ocinaplon', 'Theodrenaline', 'Pexacerfont', 'Acumapimod',
          'Pumosetrag', '2-(methylsulfanyl)-5-(thiophen-2-ylmethyl)-1H-imidazol-4-ol', 'TAK-831',
          '4-{4-[(5-hydroxy-2-methylphenyl)amino]quinolin-7-yl}-1,3-thiazole-2-carbaldehyde', 'Cadralazine',
          'GSK-1059615', '(2R)-1-{4-[(4-Anilino-5-bromo-2-pyrimidinyl)amino]phenoxy}-3-(dimethylamino)-2-propanol',
          'Harmine', 'Proquazone', 'Viminol', 'GSK2798745', 'Eletriptan', 'Imiquimod', 'Losmapimod',
          '3-{[(4-methylphenyl)sulfonyl]amino}propyl pyridin-4-ylcarbamate',
          '4-(1h-Imidazol-4-Yl)-3-(5-Ethyl-2,4-Dihydroxy-Phenyl)-1h-Pyrazole', 'ONO-2952', 'Asciminib', 'Lifirafenib',
          'Metoserpate', 'Aminophenazone', 'Pirfenidone', 'Pelitinib', 'PF-05089771', 'Moxaverine', 'Azosemide',
          '5-[3-(2-METHOXYPHENYL)-1H-PYRROLO[2,3-B]PYRIDIN-5-YL]-N,N-DIMETHYLPYRIDINE-3-CARBOXAMIDE', 'Mefloquine',
          'Sofinicline', '4-(2-(3-Propyl-(1,2,4)oxadiazol-5-yl)-vinyl)-benzene-1,2-diol', 'Capadenoson', 'Azintamide',
          'Carbimazole', 'Tenonitrozole',
          '5-(5-chloro-2,4-dihydroxyphenyl)-N-ethyl-4-[4-(morpholin-4-ylmethyl)phenyl]isoxazole-3-carboxamide',
          '3,5-Diaminophthalhydrazide', 'Pirlindole', '5,8-dimethoxy-1,4-dimethylquinolin-2(1H)-one',
          'N-[2-(6-AMINO-4-METHYLPYRIDIN-2-YL)ETHYL]-4-CYANOBENZAMIDE',
          '1-[2-(3-ACETYL-2-HYDROXY-6-METHOXY-PHENYL)-CYCLOPROPYL]-3-(5-CYANO-PYRIDIN-2-YL)-THIOUREA', 'Tipifarnib',
          'Balipodect', 'PF-00356231',
          'N-{3-[(7ar,12as,12bs)-7-Oxo-1,3,4,6,7,7a,12a,12b-Octahydroindolo[2,3-a]Quinolizin-12(2h)-Yl]Propyl}Propane-2-Sulfonamide',
          '4-METHYL-PENTANOIC ACID {1-[4-GUANIDINO-1-(THIAZOLE-2-CARBONYL)-BUTYLCARBAMOYL]-2-METHYL-PROPYL}-AMIDE',
          'Cyclo(his-pro)', 'Toreforant', 'Pruvanserin',
          '[PHENYLALANINYL-PROLINYL]-[2-(PYRIDIN-4-YLAMINO)-ETHYL]-AMINE', 'Disperse Blue 106', 'Ispronicline',
          '5-{4-[(3,5-DIFLUOROBENZYL)AMINO]PHENYL}-6-ETHYLPYRIMIDINE-2,4-DIAMINE',
          '4-{4-[4-(3-AMINOPROPOXY)PHENYL]-1H-PYRAZOL-5-YL}-6-CHLOROBENZENE-1,3-DIOL', 'Fipamezole',
          '5-Bromonicotinamide', '3-methyl-N-(pyridin-4-ylmethyl)imidazo[1,2-a]pyrazin-8-amine', '5-Methylpyrrole',
          'AZD-9496', 'Zuranolone', '5-[1-(4-methoxyphenyl)-1H-benzimidazol-6-yl]-1,3,4-oxadiazole-2(3H)-thione',
          'AZD-1236', 'Naratriptan', 'Harmaline', 'Rezatomidine', 'Florasulam',
          '9-HYDROXY-6-(3-HYDROXYPROPYL)-4-(2-METHOXYPHENYL)PYRROLO[3,4-C]CARBAZOLE-1,3(2H,6H)-DIONE',
          'Dasolampanel etibutil', 'N-[2-(5-methyl-4H-1,2,4-triazol-3-yl)phenyl]-7H-pyrrolo[2,3-d]pyrimidin-4-amine',
          'Verucerfont', 'Bunazosin',
          '(2R)-4-[(8R)-8-METHYL-2-(TRIFLUOROMETHYL)-5,6-DIHYDRO[1,2,4]TRIAZOLO[1,5-A]PYRAZIN-7(8H)-YL]-4-OXO-1-(2,4,5-TRIFLUOROPHENYL)BUTAN-2-AMINE',
          'Etrasimod', 'GSK-3117391', 'Clortermine', 'Epibatidine', 'Molidustat', 'Oxaprozin', 'Vebreltinib',
          '2-Isobutyl-3-methoxypyrazine', 'CB-103',
          '7-Methoxy-8-[1-(Methylsulfonyl)-1h-Pyrazol-4-Yl]Naphthalene-2-Carboximidamide', 'Fenpyroximate',
          'Trilaciclib', '1-(5-BROMO-PYRIDIN-2-YL)-3-[2-(6-FLUORO-2-HYDROXY-3-PROPIONYL-PHENYL)-CYCLOPROPYL]-UREA',
          'Ifidancitinib', 'AMG-319', 'Gandotinib', 'Navoximod', 'Xylose-Derived Imidazole', 'Oglemilast', 'AZD-5069',
          'Dapiprazole', 'RAD-140', 'N-[(1R)-3-(4-HYDROXYPHENYL)-1-METHYLPROPYL]-2-(2-PHENYL-1H-INDOL-3-YL)ACETAMIDE',
          'Cefatrizine',
          'N-((1R,2R)-2-(5-CHLORO-1H-INDOLE-2-CARBOXAMIDO)CYCLOHEXYL)-5-METHYL-4,5,6,7-TETRAHYDROTHIAZOLO[5,4-C]PYRIDINE-2-CARBOXAMIDE',
          'Sapanisertib', 'Aptazapine', 'PX-12', 'Vosaroxin', '4-Iodopyrazole',
          '2,4-Diamino-6-Phenyl-5,6,7,8,-Tetrahydropteridine', 'N,N-DIMETHYL(5-(PYRIDIN-3-YL)FURAN-2-YL)METHANAMINE',
          "N-{(3R,4S)-4-[(6-amino-4-methylpyridin-2-yl)methyl]pyrrolidin-3-yl}-N'-(3-chlorobenzyl)ethane-1,2-diamine",
          'D-157495', 'N-{2-methyl-5-[(6-phenylpyrimidin-4-yl)amino]phenyl}methanesulfonamide', 'Ethionamide',
          'GSK-2881078', 'PF-06282999', 'Vorolanib', 'BMS-919373', 'BTRX-246040', 'Bendamustine', 'Letrazuril',
          'HYDROXY[3-(6-METHYLPYRIDIN-2-YL)PROPYL]FORMAMIDE',
          '4-(2-(1H-IMIDAZOL-4-YL)ETHYLAMINO)-2-(PHENYLAMINO)PYRAZOLO[1,5-A][1,3,5]TRIAZINE-8-CARBONITRILE',
          'Nelotanserin', 'Cicletanine',
          '5-CHLORO-THIOPHENE-2-CARBOXYLIC ACID ((3S,4S)-4-FLUORO- 1-{[2-FLUORO-4-(2-OXO-2H-PYRIDIN-1-YL)-PHENYLCARBAMOYL]-METHYL}-PYRROLIDIN-3-YL)-AMIDE',
          'Zilpaterol', 'Duvelisib', 'F-15599', 'Pindolol', 'Reldesemtiv', 'Methylisothiazolinone', 'MK-8245',
          'Brensocatib', '2-Amino-6-Chloropyrazine', '(2-AMINO-1,3-OXAZOL-5-YL)-(3-BROMOPHENYL)METHANONE', 'Foslinanib',
          'GSK-356278', 'Lanabecestat', '9-Deazahypoxanthine', 'INCB-057643', 'Surinabant', 'PXT 3003',
          '3-(6-HYDROXY-NAPHTHALEN-2-YL)-BENZO[D]ISOOXAZOL-6-OL', 'Oltipraz', 'ONO-8539', 'Chlorzoxazone',
          'Deferiprone', 'SB-705498', 'Nidufexor',
          '7-amino-2-tert-butyl-4-{[2-(1H-imidazol-4-yl)ethyl]amino}pyrido[2,3-d]pyrimidine-6-carboxamide',
          'Imiglitazar', 'Rimacalib', 'Tetrazolyl Histidine', 'JNJ-54175446', 'Foliglurax', 'CH-5132799',
          '4-{2-[(7-amino-2-furan-2-yl[1,2,4]triazolo[1,5-a][1,3,5]triazin-5-yl)amino]ethyl}phenol',
          'N-[(2R)-2-{[(2S)-2-(1,3-benzoxazol-2-yl)pyrrolidin-1-yl]carbonyl}hexyl]-N-hydroxyformamide',
          'Acid yellow 54 free acid', '3-phenyl-5-(1H-pyrazol-3-yl)isoxazole', 'MK-5108', 'GLPG-1205', 'Pioglitazone',
          '4-[1-allyl-7-(trifluoromethyl)-1H-indazol-3-yl]benzene-1,3-diol', 'MK-212', 'Bamifylline',
          '2-(Sec-Butyl)Thiazole', 'Rivanicline', '5-benzyl-1,3-thiazol-2-amine', 'Sitamaquine',
          '(2R,3R)-N^1^-[(1S)-2,2-DIMETHYL-1-(METHYLCARBAMOYL)PROPYL]-N^4^-HYDROXY-2-(2-METHYLPROPYL)-3-{[(1,3-THIAZOL-2-YLCARBONYL)AMINO]METHYL}BUTANEDIAMIDE',
          'Cinchocaine', 'Simurosertib', 'BTRX-335140', 'Afuresertib', 'Enitociclib', 'CCX-140',
          "N-{(3S,4S)-4-[(6-AMINO-4-METHYLPYRIDIN-2-YL)METHYL]PYRROLIDIN-3-YL}-N'-(4-CHLOROBENZYL)ETHANE-1,2-DIAMINE",
          '(7as,12ar,12bs)-1,2,3,4,7a,12,12a,12b-Octahydroindolo[2,3-a]Quinolizin-7(6h)-One', 'CAN-508', 'Cinalukast',
          '(1s,2s)-1-Amino-1-(1,3-Thiazol-2-Yl)Propan-2-Ol', 'Zolpidem', 'Galeterone', 'Mercaptopurine', 'GZ-389988',
          'Varespladib methyl', 'Farampator', 'PD-173952', 'Pirbuterol', 'Amsacrine', 'AZD-0328', 'Umifenovir',
          '(S)-wiskostatin', '4-(6-{[(1R)-1-(hydroxymethyl)propyl]amino}imidazo[1,2-b]pyridazin-3-yl)benzoic acid',
          'CC-115', 'Protionamide', 'Linrodostat', 'LGH-447', 'Henatinib', 'Alprazolam', 'MK-886', 'MK-6186',
          '6-[2-(1H-INDOL-6-YL)ETHYL]PYRIDIN-2-AMINE',
          '3-(4-{2-[2-(2-bromo-acetylamino)-ethyldisulfanyl]-ethylcarbamoyl}-cyclohexylcarbamoyl)-pyrazine-2-carboxylic acid',
          '6-MORPHOLIN-4-YL-9H-PURINE',
          'N-hydroxy-5-[(3-phenyl-5,6-dihydroimidazo[1,2-a]pyrazin-7(8H)-yl)carbonyl]thiophene-2-carboxamide',
          '5-(Aminomethyl)-6-(2,4-Dichlorophenyl)-2-(3,5-Dimethoxyphenyl)Pyrimidin-4-Amine', 'Tulrampator', 'JAB-3068',
          'Ganaplacide', 'Prenoxdiazine',
          '1-[(2S)-4-(5-phenyl-1H-pyrazolo[3,4-b]pyridin-4-yl)morpholin-2-yl]methanamine', 'Trapidil',
          '4-(2-Thienyl)-1-(4-Methylbenzyl)-1h-Imidazole', 'Elvitegravir', 'Glumetinib', 'Fenamole', 'Lonidamine',
          'CXD101', 'GSK-2982772', 'Olomoucine', 'R-82913', 'Etodolac', 'Zonisamide', 'Amlexanox', 'PF-05212377',
          'Arotinolol', 'Thiacloprid', 'N-(2-Aminoethyl)-5-Chloroisoquinoline-8-Sulfonamide', 'Thiohexam', 'LY-2811376',
          'BOS172722', 'N-[2-(METHYLAMINO)ETHYL]-5-ISOQUINOLINESULFONAMIDE', 'Samuraciclib', 'Amuvatinib', 'Pecavaptan',
          'Bradanicline', '2-[(2-methoxy-5-methylphenoxy)methyl]pyridine', 'AZD-1940', 'Belotecan', 'Cenobamate',
          '2-(1H-pyrrol-1-ylcarbonyl)benzene-1,3,5-triol', 'Latrepirdine',
          'Carboxymethylthio-3-(3-Chlorophenyl)-1,2,4-Oxadiazol', 'Terevalefim', 'Chlormidazole', 'Rupatadine',
          'Tolmetin', 'Lanifibranor', 'N1-CYCLOPENTYL-N2-(THIAZOL-2-YL)OXALAMIDE', 'Benzimidazole', 'GSK-239512',
          'PF-06260414', 'Etomidate', 'Azatadine', 'RO-5028442', 'Abiraterone', 'arazepide', 'Upadacitinib', '2X-121',
          'CRA_10972', 'ABX-464', 'Tegoprazan', 'Tizanidine',
          '4-{[5-chloro-4-(1H-indol-3-yl)pyrimidin-2-yl]amino}-N-ethylpiperidine-1-carboxamide', 'Anastrozole',
          'Uracil', 'Nedisertib', 'Cilostazol', '1-(3,5-DICHLOROPHENYL)-5-METHYL-1H-1,2,4-TRIAZOLE-3-CARBOXYLIC ACID',
          'Fluconazole', 'Ipatasertib', 'N-acetylhistamine', 'Tinostamustine', 'Setipiprant', 'R-1487', 'GDC-0134',
          'Tolimidone', '(2E)-N-hydroxy-3-[1-methyl-4-(phenylacetyl)-1H-pyrrol-2-yl]prop-2-enamide',
          '3-[(2,2-DIMETHYLPROPANOYL)AMINO]-N-1,3-THIAZOL-2-YLPYRIDINE-2-CARBOXAMIDE', 'Pozanicline', 'Veliparib',
          'Ondelopran', 'IQP-0528', 'Ranirestat', '7-Nitroindazole', 'Phosphonoacetohydroxamic Acid', 'Perampanel',
          'N-phenyl-1H-pyrrolo[2,3-b]pyridin-3-amine', 'Tazobactam',
          'N~4~-methyl-N~4~-(3-methyl-1H-indazol-6-yl)-N~2~-(3,4,5-trimethoxyphenyl)pyrimidine-2,4-diamine',
          'Mercaptocarboxylate Inhibitor', 'LY-3200882', 'Flosequinan', 'PF-04691502', 'GSK-2018682', 'Cefazedone',
          '{3-[(5-CHLORO-1,3-BENZOTHIAZOL-2-YL)METHYL]-2,4-DIOXO-3,4-DIHYDROPYRIMIDIN-1(2H)-YL}ACETIC ACID',
          'Selgantolimod',
          '2,2,2-TRIFLUORO-1-{5-[(3-PHENYL-5,6-DIHYDROIMIDAZO[1,2-A]PYRAZIN-7(8H)-YL)CARBONYL]THIOPHEN-2-YL}ETHANE-1,1-DIOL',
          'Isradipine', 'Gluco-Phenylimidazole',
          "TRW3-(2-AMINO-3-HYDROXY-PROPYL)-6-(N'-CYCLOHEXYL-HYDRAZINO)OCTAHYDRO-INDOL-7-OL", 'Tucidinostat',
          'Pemigatinib', '5-CHLORO-6-METHYL-N-(2-PHENYLETHYL)-2-PYRIDIN-2-YLPYRIMIDIN-4-AMINE',
          '[1-(4-Fluorobenzyl)Cyclobutyl]Methyl (1s)-1-[Oxo(1h-Pyrazol-5-Ylamino)Acetyl]Pentylcarbamate',
          '(1S,3R,6S)-4-oxo-6-{4-[(2-phenylquinolin-4-yl)methoxy]phenyl}-5-azaspiro[2.4]heptane-1-carboxylic acid',
          '1-[(2S)-4-(5-BROMO-1H-PYRAZOLO[3,4-B]PYRIDIN-4-YL)MORPHOLIN-2-YL]METHANAMINE',
          'ethyl 3-[(E)-2-amino-1-cyanoethenyl]-6,7-dichloro-1-methyl-1H-indole-2-carboxylate',
          '4-(5,11-DIOXO-5H-INDENO[1,2-C]ISOQUINOLIN-6(11H)-YL)BUTANOATE'
          ]

val5 = ['(2Z)-1-(5-Chloro-1H-indol-3-yl)-3-hydroxy-3-(1H-tetrazol-5-yl)-2-propen-1-one', 'Sisunatovir',
        'Quinacrine mustard', 'Verdiperstat', 'Celecoxib', '3-{3-[(DIMETHYLAMINO)METHYL]-1H-INDOL-7-YL}PROPAN-1-OL',
        'Ethacridine', '(2S)-2-AMINO-4-(METHYLSULFANYL)-1-(1,3-THIAZOL-2-YL)BUTANE-1,1-DIOL',
        'N-(2,6-dimethylphenyl)-5-phenylimidazo[1,5-a]pyrazin-8-amine', 'Ezutromid', 'Nedocromil', 'Ziresovir',
        'Tavapadon', 'Pardoprunox', 'K-134', 'Pyrrole-2-Carboxylate', 'Sertaconazole', '6-Methylpurine',
        '2-amino-5-[3-(1-ethyl-1H-pyrazol-5-yl)-1H-pyrrolo[2,3-b]pyridin-5-yl]-N,N-dimethylbenzamide', 'Linopirdine',
        'Vecabrutinib', 'Venglustat',
        '4-{[1-Methyl-5-(2-Methyl-Benzoimidazol-1-Ylmethyl)-1h-Benzoimidazol-2-Ylmethyl]-Amino}-Benzamidine',
        'N,4-dimethyl-3-[(1-phenyl-1H-pyrazolo[3,4-d]pyrimidin-4-yl)amino]benzamide', 'NP-G2-044', 'Olutasidenib',
        'Simpinicline', 'PP-121', 'PHENYLAMINOIMIDAZO(1,2-ALPHA)PYRIDINE', 'GW810781', 'PSI-697',
        'N-(4-{[(3S)-3-(dimethylamino)pyrrolidin-1-yl]carbonyl}phenyl)-5-fluoro-4-[2-methyl-1-(1-methylethyl)-1H-imidazol-5-yl]pyrimidin-2-amine',
        'Hydronidone', 'Pyrazole', 'Bromazepam', 'AZD-8418',
        '4-({[4-(3-METHYLBENZOYL)PYRIDIN-2-YL]AMINO}METHYL)BENZENECARBOXIMIDAMIDE', 'Fanapanel', 'Dianicline',
        '[2-(5-Mercapto-[1,3,4]thiadiazol-2-ylcarbamoyl)-1-phenyl-ethyl]-carbamic acid benzyl ester',
        "(4R)-7,8-dichloro-1',9-dimethyl-1-oxo-1,2,4,9-tetrahydrospiro[beta-carboline-3,4'-piperidine]-4-carbonitrile",
        'E-7820', '7-(aminomethyl)-6-(2-chlorophenyl)-1-methyl-1H-benzimidazole-5-carbonitrile', 'AZD-1386',
        'Imaradenant', 'Pranoprofen', 'Tisopurine', 'Decernotinib'
        ]

test5 = ['Triazolopyridine', '9-ACETYL-2,3,4,9-TETRAHYDRO-1H-CARBAZOL-1-ONE', 'PF-06700841', 'Rogaratinib', 'ABT-288',
         'Buparlisib', '4-chloro-6-{5-[(2-morpholin-4-ylethyl)amino]-1,2-benzisoxazol-3-yl}benzene-1,3-diol',
         'Cyanocinnoline', 'Zamicastat',
         'Cp403700, (S)-1-{2-[(5-Chloro-1h-Indole-2-Carbonyl)-Amino]-3-Phenyl-Propionyl}-Azetidine-3-Carboxylate',
         'Isatoic anhydride', 'UK-500001', 'Droxicam', 'LY-3023414', 'Urapidil',
         '5-(7-(4-(4,5-dihydro-2-oxazolyl)phenoxy)heptyl)-3-methyl isoxazole',
         '2-CHLORO-N-[(3R)-2-OXO-1,2,3,4-TETRAHYDROQUINOLIN-3-YL]-6H-THIENO[2,3-B]PYRROLE-5-CARBOXAMIDE', 'Triapine',
         'PF-04991532', '1-[1-(3-aminophenyl)-3-tert-butyl-1H-pyrazol-5-yl]-3-phenylurea', 'Benzoylformic Acid',
         'N~2~-1,3-BENZOXAZOL-2-YL-3-CYCLOHEXYL-N-{2-[(4-METHOXYPHENYL)AMINO]ETHYL}-L-ALANINAMIDE', 'PF-05175157',
         '1,N6-Ethenoadenine', '6-[4-(2-fluorophenyl)-1,3-oxazol-5-yl]-N-(1-methylethyl)-1,3-benzothiazol-2-amine',
         'Rilematovir', '(2S)-2-{[3-(3-aminophenyl)imidazo[1,2-b]pyridazin-6-yl]amino}-3-methylbutan-1-ol',
         'Pantoprazole', 'Piclozotan', 'LGD2941', 'Ocinaplon', 'Theodrenaline', 'Pexacerfont', 'Acumapimod',
         'Pumosetrag', '2-(methylsulfanyl)-5-(thiophen-2-ylmethyl)-1H-imidazol-4-ol', 'TAK-831',
         '4-{4-[(5-hydroxy-2-methylphenyl)amino]quinolin-7-yl}-1,3-thiazole-2-carbaldehyde', 'Cadralazine',
         'GSK-1059615', '(2R)-1-{4-[(4-Anilino-5-bromo-2-pyrimidinyl)amino]phenoxy}-3-(dimethylamino)-2-propanol',
         'Harmine', 'Proquazone', 'Viminol', 'GSK2798745', 'Eletriptan', 'Imiquimod', 'Losmapimod'
         ]

train5 = ['3-{[(4-methylphenyl)sulfonyl]amino}propyl pyridin-4-ylcarbamate',
          '4-(1h-Imidazol-4-Yl)-3-(5-Ethyl-2,4-Dihydroxy-Phenyl)-1h-Pyrazole', 'ONO-2952', 'Asciminib', 'Lifirafenib',
          'Metoserpate', 'Aminophenazone', 'Pirfenidone', 'Pelitinib', 'PF-05089771', 'Moxaverine', 'Azosemide',
          '5-[3-(2-METHOXYPHENYL)-1H-PYRROLO[2,3-B]PYRIDIN-5-YL]-N,N-DIMETHYLPYRIDINE-3-CARBOXAMIDE', 'Mefloquine',
          'Sofinicline', '4-(2-(3-Propyl-(1,2,4)oxadiazol-5-yl)-vinyl)-benzene-1,2-diol', 'Capadenoson', 'Azintamide',
          'Carbimazole', 'Tenonitrozole',
          '5-(5-chloro-2,4-dihydroxyphenyl)-N-ethyl-4-[4-(morpholin-4-ylmethyl)phenyl]isoxazole-3-carboxamide',
          '3,5-Diaminophthalhydrazide', 'Pirlindole', '5,8-dimethoxy-1,4-dimethylquinolin-2(1H)-one',
          'N-[2-(6-AMINO-4-METHYLPYRIDIN-2-YL)ETHYL]-4-CYANOBENZAMIDE',
          '1-[2-(3-ACETYL-2-HYDROXY-6-METHOXY-PHENYL)-CYCLOPROPYL]-3-(5-CYANO-PYRIDIN-2-YL)-THIOUREA', 'Tipifarnib',
          'Balipodect', 'PF-00356231',
          'N-{3-[(7ar,12as,12bs)-7-Oxo-1,3,4,6,7,7a,12a,12b-Octahydroindolo[2,3-a]Quinolizin-12(2h)-Yl]Propyl}Propane-2-Sulfonamide',
          '4-METHYL-PENTANOIC ACID {1-[4-GUANIDINO-1-(THIAZOLE-2-CARBONYL)-BUTYLCARBAMOYL]-2-METHYL-PROPYL}-AMIDE',
          'Cyclo(his-pro)', 'Toreforant', 'Pruvanserin',
          '[PHENYLALANINYL-PROLINYL]-[2-(PYRIDIN-4-YLAMINO)-ETHYL]-AMINE', 'Disperse Blue 106', 'Ispronicline',
          '5-{4-[(3,5-DIFLUOROBENZYL)AMINO]PHENYL}-6-ETHYLPYRIMIDINE-2,4-DIAMINE',
          '4-{4-[4-(3-AMINOPROPOXY)PHENYL]-1H-PYRAZOL-5-YL}-6-CHLOROBENZENE-1,3-DIOL', 'Fipamezole',
          '5-Bromonicotinamide', '3-methyl-N-(pyridin-4-ylmethyl)imidazo[1,2-a]pyrazin-8-amine', '5-Methylpyrrole',
          'AZD-9496', 'Zuranolone', '5-[1-(4-methoxyphenyl)-1H-benzimidazol-6-yl]-1,3,4-oxadiazole-2(3H)-thione',
          'AZD-1236', 'Naratriptan', 'Harmaline', 'Rezatomidine', 'Florasulam',
          '9-HYDROXY-6-(3-HYDROXYPROPYL)-4-(2-METHOXYPHENYL)PYRROLO[3,4-C]CARBAZOLE-1,3(2H,6H)-DIONE',
          'Dasolampanel etibutil', 'N-[2-(5-methyl-4H-1,2,4-triazol-3-yl)phenyl]-7H-pyrrolo[2,3-d]pyrimidin-4-amine',
          'Verucerfont', 'Bunazosin',
          '(2R)-4-[(8R)-8-METHYL-2-(TRIFLUOROMETHYL)-5,6-DIHYDRO[1,2,4]TRIAZOLO[1,5-A]PYRAZIN-7(8H)-YL]-4-OXO-1-(2,4,5-TRIFLUOROPHENYL)BUTAN-2-AMINE',
          'Etrasimod', 'GSK-3117391', 'Clortermine', 'Epibatidine', 'Molidustat', 'Oxaprozin', 'Vebreltinib',
          '2-Isobutyl-3-methoxypyrazine', 'CB-103',
          '7-Methoxy-8-[1-(Methylsulfonyl)-1h-Pyrazol-4-Yl]Naphthalene-2-Carboximidamide', 'Fenpyroximate',
          'Trilaciclib', '1-(5-BROMO-PYRIDIN-2-YL)-3-[2-(6-FLUORO-2-HYDROXY-3-PROPIONYL-PHENYL)-CYCLOPROPYL]-UREA',
          'Ifidancitinib', 'AMG-319', 'Gandotinib', 'Navoximod', 'Xylose-Derived Imidazole', 'Oglemilast', 'AZD-5069',
          'Dapiprazole', 'RAD-140', 'N-[(1R)-3-(4-HYDROXYPHENYL)-1-METHYLPROPYL]-2-(2-PHENYL-1H-INDOL-3-YL)ACETAMIDE',
          'Cefatrizine',
          'N-((1R,2R)-2-(5-CHLORO-1H-INDOLE-2-CARBOXAMIDO)CYCLOHEXYL)-5-METHYL-4,5,6,7-TETRAHYDROTHIAZOLO[5,4-C]PYRIDINE-2-CARBOXAMIDE',
          'Sapanisertib', 'Aptazapine', 'PX-12', 'Vosaroxin', '4-Iodopyrazole',
          '2,4-Diamino-6-Phenyl-5,6,7,8,-Tetrahydropteridine', 'N,N-DIMETHYL(5-(PYRIDIN-3-YL)FURAN-2-YL)METHANAMINE',
          "N-{(3R,4S)-4-[(6-amino-4-methylpyridin-2-yl)methyl]pyrrolidin-3-yl}-N'-(3-chlorobenzyl)ethane-1,2-diamine",
          'D-157495', 'N-{2-methyl-5-[(6-phenylpyrimidin-4-yl)amino]phenyl}methanesulfonamide', 'Ethionamide',
          'GSK-2881078', 'PF-06282999', 'Vorolanib', 'BMS-919373', 'BTRX-246040', 'Bendamustine', 'Letrazuril',
          'HYDROXY[3-(6-METHYLPYRIDIN-2-YL)PROPYL]FORMAMIDE',
          '4-(2-(1H-IMIDAZOL-4-YL)ETHYLAMINO)-2-(PHENYLAMINO)PYRAZOLO[1,5-A][1,3,5]TRIAZINE-8-CARBONITRILE',
          'Nelotanserin', 'Cicletanine',
          '5-CHLORO-THIOPHENE-2-CARBOXYLIC ACID ((3S,4S)-4-FLUORO- 1-{[2-FLUORO-4-(2-OXO-2H-PYRIDIN-1-YL)-PHENYLCARBAMOYL]-METHYL}-PYRROLIDIN-3-YL)-AMIDE',
          'Zilpaterol', 'Duvelisib', 'F-15599', 'Pindolol', 'Reldesemtiv', 'Methylisothiazolinone', 'MK-8245',
          'Brensocatib', '2-Amino-6-Chloropyrazine', '(2-AMINO-1,3-OXAZOL-5-YL)-(3-BROMOPHENYL)METHANONE', 'Foslinanib',
          'GSK-356278', 'Lanabecestat', '9-Deazahypoxanthine', 'INCB-057643', 'Surinabant', 'PXT 3003',
          '3-(6-HYDROXY-NAPHTHALEN-2-YL)-BENZO[D]ISOOXAZOL-6-OL', 'Oltipraz', 'ONO-8539', 'Chlorzoxazone',
          'Deferiprone', 'SB-705498', 'Nidufexor',
          '7-amino-2-tert-butyl-4-{[2-(1H-imidazol-4-yl)ethyl]amino}pyrido[2,3-d]pyrimidine-6-carboxamide',
          'Imiglitazar', 'Rimacalib', 'Tetrazolyl Histidine', 'JNJ-54175446', 'Foliglurax', 'CH-5132799',
          '4-{2-[(7-amino-2-furan-2-yl[1,2,4]triazolo[1,5-a][1,3,5]triazin-5-yl)amino]ethyl}phenol',
          'N-[(2R)-2-{[(2S)-2-(1,3-benzoxazol-2-yl)pyrrolidin-1-yl]carbonyl}hexyl]-N-hydroxyformamide',
          'Acid yellow 54 free acid', '3-phenyl-5-(1H-pyrazol-3-yl)isoxazole', 'MK-5108', 'GLPG-1205', 'Pioglitazone',
          '4-[1-allyl-7-(trifluoromethyl)-1H-indazol-3-yl]benzene-1,3-diol', 'MK-212', 'Bamifylline',
          '2-(Sec-Butyl)Thiazole', 'Rivanicline', '5-benzyl-1,3-thiazol-2-amine', 'Sitamaquine',
          '(2R,3R)-N^1^-[(1S)-2,2-DIMETHYL-1-(METHYLCARBAMOYL)PROPYL]-N^4^-HYDROXY-2-(2-METHYLPROPYL)-3-{[(1,3-THIAZOL-2-YLCARBONYL)AMINO]METHYL}BUTANEDIAMIDE',
          'Cinchocaine', 'Simurosertib', 'BTRX-335140', 'Afuresertib', 'Enitociclib', 'CCX-140',
          "N-{(3S,4S)-4-[(6-AMINO-4-METHYLPYRIDIN-2-YL)METHYL]PYRROLIDIN-3-YL}-N'-(4-CHLOROBENZYL)ETHANE-1,2-DIAMINE",
          '(7as,12ar,12bs)-1,2,3,4,7a,12,12a,12b-Octahydroindolo[2,3-a]Quinolizin-7(6h)-One', 'CAN-508', 'Cinalukast',
          '(1s,2s)-1-Amino-1-(1,3-Thiazol-2-Yl)Propan-2-Ol', 'Zolpidem', 'Galeterone', 'Mercaptopurine', 'GZ-389988',
          'Varespladib methyl', 'Farampator', 'PD-173952', 'Pirbuterol', 'Amsacrine', 'AZD-0328', 'Umifenovir',
          '(S)-wiskostatin', '4-(6-{[(1R)-1-(hydroxymethyl)propyl]amino}imidazo[1,2-b]pyridazin-3-yl)benzoic acid',
          'CC-115', 'Protionamide', 'Linrodostat', 'LGH-447', 'Henatinib', 'Alprazolam', 'MK-886', 'MK-6186',
          '6-[2-(1H-INDOL-6-YL)ETHYL]PYRIDIN-2-AMINE',
          '3-(4-{2-[2-(2-bromo-acetylamino)-ethyldisulfanyl]-ethylcarbamoyl}-cyclohexylcarbamoyl)-pyrazine-2-carboxylic acid',
          '6-MORPHOLIN-4-YL-9H-PURINE',
          'N-hydroxy-5-[(3-phenyl-5,6-dihydroimidazo[1,2-a]pyrazin-7(8H)-yl)carbonyl]thiophene-2-carboxamide',
          '5-(Aminomethyl)-6-(2,4-Dichlorophenyl)-2-(3,5-Dimethoxyphenyl)Pyrimidin-4-Amine', 'Tulrampator', 'JAB-3068',
          'Ganaplacide', 'Prenoxdiazine',
          '1-[(2S)-4-(5-phenyl-1H-pyrazolo[3,4-b]pyridin-4-yl)morpholin-2-yl]methanamine', 'Trapidil',
          '4-(2-Thienyl)-1-(4-Methylbenzyl)-1h-Imidazole', 'Elvitegravir', 'Glumetinib', 'Fenamole', 'Lonidamine',
          'CXD101', 'GSK-2982772', 'Olomoucine', 'R-82913', 'Etodolac', 'Zonisamide', 'Amlexanox', 'PF-05212377',
          'Arotinolol', 'Thiacloprid', 'N-(2-Aminoethyl)-5-Chloroisoquinoline-8-Sulfonamide', 'Thiohexam', 'LY-2811376',
          'BOS172722', 'N-[2-(METHYLAMINO)ETHYL]-5-ISOQUINOLINESULFONAMIDE', 'Samuraciclib', 'Amuvatinib', 'Pecavaptan',
          'Bradanicline', '2-[(2-methoxy-5-methylphenoxy)methyl]pyridine', 'AZD-1940', 'Belotecan', 'Cenobamate',
          '2-(1H-pyrrol-1-ylcarbonyl)benzene-1,3,5-triol', 'Latrepirdine',
          'Carboxymethylthio-3-(3-Chlorophenyl)-1,2,4-Oxadiazol', 'Terevalefim', 'Chlormidazole', 'Rupatadine',
          'Tolmetin', 'Lanifibranor', 'N1-CYCLOPENTYL-N2-(THIAZOL-2-YL)OXALAMIDE', 'Benzimidazole', 'GSK-239512',
          'PF-06260414', 'Etomidate', 'Azatadine', 'RO-5028442', 'Abiraterone', 'arazepide', 'Upadacitinib', '2X-121',
          'CRA_10972', 'ABX-464', 'Tegoprazan', 'Tizanidine',
          '4-{[5-chloro-4-(1H-indol-3-yl)pyrimidin-2-yl]amino}-N-ethylpiperidine-1-carboxamide', 'Anastrozole',
          'Uracil', 'Nedisertib', 'Cilostazol', '1-(3,5-DICHLOROPHENYL)-5-METHYL-1H-1,2,4-TRIAZOLE-3-CARBOXYLIC ACID',
          'Fluconazole', 'Ipatasertib', 'N-acetylhistamine', 'Tinostamustine', 'Setipiprant', 'R-1487', 'GDC-0134',
          'Tolimidone', '(2E)-N-hydroxy-3-[1-methyl-4-(phenylacetyl)-1H-pyrrol-2-yl]prop-2-enamide',
          '3-[(2,2-DIMETHYLPROPANOYL)AMINO]-N-1,3-THIAZOL-2-YLPYRIDINE-2-CARBOXAMIDE', 'Pozanicline', 'Veliparib',
          'Ondelopran', 'IQP-0528', 'Ranirestat', '7-Nitroindazole', 'Phosphonoacetohydroxamic Acid', 'Perampanel',
          'N-phenyl-1H-pyrrolo[2,3-b]pyridin-3-amine', 'Tazobactam',
          'N~4~-methyl-N~4~-(3-methyl-1H-indazol-6-yl)-N~2~-(3,4,5-trimethoxyphenyl)pyrimidine-2,4-diamine',
          'Mercaptocarboxylate Inhibitor', 'LY-3200882', 'Flosequinan', 'PF-04691502', 'GSK-2018682', 'Cefazedone',
          '{3-[(5-CHLORO-1,3-BENZOTHIAZOL-2-YL)METHYL]-2,4-DIOXO-3,4-DIHYDROPYRIMIDIN-1(2H)-YL}ACETIC ACID',
          'Selgantolimod',
          '2,2,2-TRIFLUORO-1-{5-[(3-PHENYL-5,6-DIHYDROIMIDAZO[1,2-A]PYRAZIN-7(8H)-YL)CARBONYL]THIOPHEN-2-YL}ETHANE-1,1-DIOL',
          'Isradipine', 'Gluco-Phenylimidazole',
          "TRW3-(2-AMINO-3-HYDROXY-PROPYL)-6-(N'-CYCLOHEXYL-HYDRAZINO)OCTAHYDRO-INDOL-7-OL", 'Tucidinostat',
          'Pemigatinib', '5-CHLORO-6-METHYL-N-(2-PHENYLETHYL)-2-PYRIDIN-2-YLPYRIMIDIN-4-AMINE',
          '[1-(4-Fluorobenzyl)Cyclobutyl]Methyl (1s)-1-[Oxo(1h-Pyrazol-5-Ylamino)Acetyl]Pentylcarbamate',
          '(1S,3R,6S)-4-oxo-6-{4-[(2-phenylquinolin-4-yl)methoxy]phenyl}-5-azaspiro[2.4]heptane-1-carboxylic acid',
          '1-[(2S)-4-(5-BROMO-1H-PYRAZOLO[3,4-B]PYRIDIN-4-YL)MORPHOLIN-2-YL]METHANAMINE',
          'ethyl 3-[(E)-2-amino-1-cyanoethenyl]-6,7-dichloro-1-methyl-1H-indole-2-carboxylate',
          '4-(5,11-DIOXO-5H-INDENO[1,2-C]ISOQUINOLIN-6(11H)-YL)BUTANOATE', '2-Methoxy-3-Isopropylpyrazine',
          '3-(1-NAPHTHYLMETHOXY)PYRIDIN-2-AMINE', 'Di-2-pyridylketone 4-cyclohexyl-4-methyl-3-thiosemicarbazone',
          'Methimazole', 'S-777469', 'ORE-1001', 'Oliceridine', 'N-6022', 'GDC-0425',
          "4-Bromo-3-(5'-Carboxy-4'-Chloro-2'-Fluorophenyl)-1-Methyl-5-Trifluoromethyl-Pyrazol", 'Tivirapine',
          'Rufinamide', '6-BENZYL-1-BENZYLOXYMETHYL-5-ISOPROPYL URACIL', 'Avadomide',
          '1-[(4S)-4-amino-5-(1,3-benzothiazol-2-yl)-5-oxopentyl]guanidine', 'Cipargamin',
          "1-(6-CYANO-3-PYRIDYLCARBONYL)-5',8'-DIFLUOROSPIRO[PIPERIDINE-4,2'(1'H)-QUINAZOLINE]-4'-AMINE",
          'JNJ-39393406', 'Piribedil', 'Imidazole', 'Pipequaline', 'Lavoltidine',
          'N-{3-[5-(1H-1,2,4-triazol-3-yl)-1H-indazol-3-yl]phenyl}furan-2-carboxamide', 'Tinengotinib', 'Uracil C-13',
          'Minodronic acid', 'Tovorafenib', 'PF-03463275', '3-(4-CHLOROPHENYL)-5-(METHYLTHIO)-4H-1,2,4-TRIAZOLE',
          'AZD-1981', 'Revaprazan', 'Rabeximod', '3-Methyladenine', '5-Aminoisoquinoline',
          '6-HYDROXY-1,3-BENZOTHIAZOLE-2-SULFONAMIDE', '1,2,4-Triazole',
          '1-(2-Chlorophenyl)-3,5-Dimethyl-1h-Pyrazole-4-Carboxylic Acid Ethyl Ester',
          '6-CHLORO-9-HYDROXY-1,3-DIMETHYL-1,9-DIHYDRO-4H-PYRAZOLO[3,4-B]QUINOLIN-4-ONE', 'Dilmapimod', 'CGS-27023',
          'N-[2-methyl-5-(methylcarbamoyl)phenyl]-2-{[(1R)-1-methylpropyl]amino}-1,3-thiazole-5-carboxamide',
          'Nitenpyram', 'AMG-900', '6-CHLORO-4-(CYCLOHEXYLOXY)-3-PROPYLQUINOLIN-2(1H)-ONE', 'Roflumilast', 'Posizolid',
          'Carboxyamidotriazole', '5-(2-chlorophenyl)-1,3,4-thiadiazole-2-sulfonamide',
          '4-CHLORO-6-(4-PIPERAZIN-1-YL-1H-PYRAZOL-5-YL)BENZENE-1,3-DIOL',
          '(2R)-1-[(5,6-DIPHENYL-7H-PYRROLO[2,3-D]PYRIMIDIN-4-YL)AMINO]PROPAN-2-OL', 'Birabresib', 'Cimicoxib',
          'Snubh-nm-333 F-18', 'Polaprezinc', 'Pradigastat',
          "(4R)-7-chloro-9-methyl-1-oxo-1,2,4,9-tetrahydrospiro[beta-carboline-3,4'-piperidine]-4-carbonitrile",
          'Paltusotine', 'Diclazuril', 'Butalamine', 'Riluzole', 'N-METHYL-1-[4-(9H-PURIN-6-YL)PHENYL]METHANAMINE',
          'Aleplasinin', 'Linzagolix', 'Gimeracil', 'BMS-488043', 'Clonazolam', 'Ciclopirox', 'Naphthoquine',
          '4-[(3-BROMO-4-O-SULFAMOYLBENZYL)(4-CYANOPHENYL)AMINO]-4H-[1,2,4]-TRIAZOLE',
          '5-[(Z)-(5-Chloro-2-oxo-1,2-dihydro-3H-indol-3-ylidene)methyl]-N,2,4-trimethyl-1H-pyrrole-3-carboxamide',
          'AG-24322', 'Pamiparib', 'PF-05241328', 'Dalfampridine', '2-Pyridinethiol', 'Mivavotinib',
          '2-(4-ETHYLPIPERAZIN-1-YL)-4-(PHENYLAMINO)PYRAZOLO[1,5-A][1,3,5]TRIAZINE-8-CARBONITRILE', 'Metralindole',
          'Lamotrigine', 'Cinoxacin', 'Anagliptin', 'Losartan', 'Allopurinol', 'Frovatriptan', 'Deucravacitinib',
          '2-(2-chloropyridin-4-yl)-4-methyl-1H-isoindole-1,3(2H)-dione', 'Irbesartan',
          '3-[2-bromo-4-(1H-pyrazolo[3,4-c]pyridazin-3-ylmethyl)phenoxy]-5-methylbenzonitrile', 'Onatasertib',
          "N'-(5-CHLORO-1,3-BENZODIOXOL-4-YL)-N-(3-MORPHOLIN-4-YLPHENYL)PYRIMIDINE-2,4-DIAMINE", 'Dactolisib',
          'Capivasertib', 'Cenerimod', 'PTI-428', 'V116517', 'AMD-070',
          '6-(2,4-DIAMINO-6-ETHYLPYRIMIDIN-5-YL)-4-(3-METHOXYPROPYL)-2,2-DIMETHYL-2H-1,4-BENZOXAZIN-3(4H)-ONE',
          'Vincamine'
          ]

with open('Whole_molecule_corrected_reagents.csv', 'r') as csv_file:
    csv_reader = list(csv.reader(csv_file, delimiter=','))[1:]
    # print(csv_reader)
    validation_set1 = []
    test_set1 = []
    training_set1 = []
    validation_set2 = []
    test_set2 = []
    training_set2 = []
    validation_set3 = []
    test_set3 = []
    training_set3 = []
    validation_set4 = []
    test_set4 = []
    training_set4 = []
    validation_set5 = []
    test_set5 = []
    training_set5 = []

    for row in csv_reader:
        for name in val1:
            if name == row[0]:
                dat = [row[0], row[5], row[3], row[8], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20]]
                validation_set1.append(dat)
    for row in csv_reader:
        for name in test1:
            if name == row[0]:
                dat = [row[0], row[5], row[3], row[8], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20]]
                test_set1.append(dat)
    for row in csv_reader:
        for name in train1:
            if name == row[0]:
                dat = [row[0], row[5], row[3], row[8], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20]]
                training_set1.append(dat)
    for row in csv_reader:
        for name in val2:
            if name == row[0]:
                dat = [row[0], row[5], row[3], row[8], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20]]
                validation_set2.append(dat)
    for row in csv_reader:
        for name in test2:
            if name == row[0]:
                dat = [row[0], row[5], row[3], row[8], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20]]
                test_set2.append(dat)
    for row in csv_reader:
        for name in train2:
            if name == row[0]:
                dat = [row[0], row[5], row[3], row[8], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20]]
                training_set2.append(dat)
    for row in csv_reader:
        for name in val3:
            if name == row[0]:
                dat = [row[0], row[5], row[3], row[8], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20]]
                validation_set3.append(dat)
    for row in csv_reader:
        for name in test3:
            if name == row[0]:
                dat = [row[0], row[5], row[3], row[8], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20]]
                test_set3.append(dat)
    for row in csv_reader:
        for name in train3:
            if name == row[0]:
                dat = [row[0], row[5], row[3], row[8], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20]]
                training_set3.append(dat)
    for row in csv_reader:
        for name in val4:
            if name == row[0]:
                dat = [row[0], row[5], row[3], row[8], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20]]
                validation_set4.append(dat)
    for row in csv_reader:
        for name in test4:
            if name == row[0]:
                dat = [row[0], row[5], row[3], row[8], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20]]
                test_set4.append(dat)
    for row in csv_reader:
        for name in train4:
            if name == row[0]:
                dat = [row[0], row[5], row[3], row[8], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20]]
                training_set4.append(dat)
    for row in csv_reader:
        for name in val5:
            if name == row[0]:
                dat = [row[0], row[5], row[3], row[8], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20]]
                validation_set5.append(dat)
    for row in csv_reader:
        for name in test5:
            if name == row[0]:
                dat = [row[0], row[5], row[3], row[8], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20]]
                test_set5.append(dat)
    for row in csv_reader:
        for name in train5:
            if name == row[0]:
                dat = [row[0], row[5], row[3], row[8], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20]]
                training_set5.append(dat)

    with open('validation1_sites_extra.csv', 'w', newline='') as vali:
        writer = csv.writer(vali)
        header = ['Compound_name', 'Reaction_Smiles', 'Site', 'hf_activation_energy', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0', 'Electron_density', 'Electrostatic_potential', 'Diamagnetic_shielding', 'Mulliken_charge', 'Bond_index', 'Gross_population', 'HF_C-H_Bond_Length', 'HF_TS_C-C_Bond_Length']
        writer.writerow(header)
        for row in validation_set1:
            writer.writerow(row)
    with open('test1_sites_extra.csv', 'w', newline='') as tes:
        writer = csv.writer(tes)
        header = ['Compound_name', 'Reaction_Smiles', 'Site', 'hf_activation_energy', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0', 'Electron_density', 'Electrostatic_potential', 'Diamagnetic_shielding', 'Mulliken_charge', 'Bond_index', 'Gross_population', 'HF_C-H_Bond_Length', 'HF_TS_C-C_Bond_Length']
        writer.writerow(header)
        for row in test_set1:
            writer.writerow(row)
    with open('training1_sites_extra.csv', 'w', newline='') as tra:
        writer = csv.writer(tra)
        header = ['Compound_name', 'Reaction_Smiles', 'Site', 'hf_activation_energy', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0', 'Electron_density', 'Electrostatic_potential', 'Diamagnetic_shielding', 'Mulliken_charge', 'Bond_index', 'Gross_population', 'HF_C-H_Bond_Length', 'HF_TS_C-C_Bond_Length']
        writer.writerow(header)
        for row in training_set1:
            writer.writerow(row)
    with open('validation2_sites_extra.csv', 'w', newline='') as vali:
        writer = csv.writer(vali)
        header = ['Compound_name', 'Reaction_Smiles', 'Site', 'hf_activation_energy', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0', 'Electron_density', 'Electrostatic_potential', 'Diamagnetic_shielding', 'Mulliken_charge', 'Bond_index', 'Gross_population', 'HF_C-H_Bond_Length', 'HF_TS_C-C_Bond_Length']
        writer.writerow(header)
        for row in validation_set2:
            writer.writerow(row)

    with open('test2_sites_extra.csv', 'w', newline='') as tes:
        writer = csv.writer(tes)
        header = ['Compound_name', 'Reaction_Smiles', 'Site', 'hf_activation_energy', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0', 'Electron_density', 'Electrostatic_potential', 'Diamagnetic_shielding', 'Mulliken_charge', 'Bond_index', 'Gross_population', 'HF_C-H_Bond_Length', 'HF_TS_C-C_Bond_Length']
        writer.writerow(header)
        for row in test_set2:
            writer.writerow(row)

    with open('training2_sites_extra.csv', 'w', newline='') as tra:
        writer = csv.writer(tra)
        header = ['Compound_name', 'Reaction_Smiles', 'Site', 'hf_activation_energy', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0', 'Electron_density', 'Electrostatic_potential', 'Diamagnetic_shielding', 'Mulliken_charge', 'Bond_index', 'Gross_population', 'HF_C-H_Bond_Length', 'HF_TS_C-C_Bond_Length']
        writer.writerow(header)
        for row in training_set2:
            writer.writerow(row)
    with open('validation3_sites_extra.csv', 'w', newline='') as vali:
        writer = csv.writer(vali)
        header = ['Compound_name', 'Reaction_Smiles', 'Site', 'hf_activation_energy', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0', 'Electron_density', 'Electrostatic_potential', 'Diamagnetic_shielding', 'Mulliken_charge', 'Bond_index', 'Gross_population', 'HF_C-H_Bond_Length', 'HF_TS_C-C_Bond_Length']
        writer.writerow(header)
        for row in validation_set3:
            writer.writerow(row)

    with open('test3_sites_extra.csv', 'w', newline='') as tes:
        writer = csv.writer(tes)
        header = ['Compound_name', 'Reaction_Smiles', 'Site', 'hf_activation_energy', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0', 'Electron_density', 'Electrostatic_potential', 'Diamagnetic_shielding', 'Mulliken_charge', 'Bond_index', 'Gross_population', 'HF_C-H_Bond_Length', 'HF_TS_C-C_Bond_Length']
        writer.writerow(header)
        for row in test_set3:
            writer.writerow(row)

    with open('training3_sites_extra.csv', 'w', newline='') as tra:
        writer = csv.writer(tra)
        header = ['Compound_name', 'Reaction_Smiles', 'Site', 'hf_activation_energy', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0', 'Electron_density', 'Electrostatic_potential', 'Diamagnetic_shielding', 'Mulliken_charge', 'Bond_index', 'Gross_population', 'HF_C-H_Bond_Length', 'HF_TS_C-C_Bond_Length']
        writer.writerow(header)
        for row in training_set3:
            writer.writerow(row)
    with open('validation4_sites_extra.csv', 'w', newline='') as vali:
        writer = csv.writer(vali)
        header = ['Compound_name', 'Reaction_Smiles', 'Site', 'hf_activation_energy', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0', 'Electron_density', 'Electrostatic_potential', 'Diamagnetic_shielding', 'Mulliken_charge', 'Bond_index', 'Gross_population', 'HF_C-H_Bond_Length', 'HF_TS_C-C_Bond_Length']
        writer.writerow(header)
        for row in validation_set4:
            writer.writerow(row)

    with open('test4_sites_extra.csv', 'w', newline='') as tes:
        writer = csv.writer(tes)
        header = ['Compound_name', 'Reaction_Smiles', 'Site', 'hf_activation_energy', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0', 'Electron_density', 'Electrostatic_potential', 'Diamagnetic_shielding', 'Mulliken_charge', 'Bond_index', 'Gross_population', 'HF_C-H_Bond_Length', 'HF_TS_C-C_Bond_Length']
        writer.writerow(header)
        for row in test_set4:
            writer.writerow(row)

    with open('training4_sites_extra.csv', 'w', newline='') as tra:
        writer = csv.writer(tra)
        header = ['Compound_name', 'Reaction_Smiles', 'Site', 'hf_activation_energy', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0', 'Electron_density', 'Electrostatic_potential', 'Diamagnetic_shielding', 'Mulliken_charge', 'Bond_index', 'Gross_population', 'HF_C-H_Bond_Length', 'HF_TS_C-C_Bond_Length']
        writer.writerow(header)
        for row in training_set4:
            writer.writerow(row)
    with open('validation5_sites_extra.csv', 'w', newline='') as vali:
        writer = csv.writer(vali)
        header = ['Compound_name', 'Reaction_Smiles', 'Site', 'hf_activation_energy', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0', 'Electron_density', 'Electrostatic_potential', 'Diamagnetic_shielding', 'Mulliken_charge', 'Bond_index', 'Gross_population', 'HF_C-H_Bond_Length', 'HF_TS_C-C_Bond_Length']
        writer.writerow(header)
        for row in validation_set5:
            writer.writerow(row)

    with open('test5_sites_extra.csv', 'w', newline='') as tes:
        writer = csv.writer(tes)
        header = ['Compound_name', 'Reaction_Smiles', 'Site', 'hf_activation_energy', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0', 'Electron_density', 'Electrostatic_potential', 'Diamagnetic_shielding', 'Mulliken_charge', 'Bond_index', 'Gross_population', 'HF_C-H_Bond_Length', 'HF_TS_C-C_Bond_Length']
        writer.writerow(header)
        for row in test_set5:
            writer.writerow(row)

    with open('training5_sites_extra.csv', 'w', newline='') as tra:
        writer = csv.writer(tra)
        header = ['Compound_name', 'Reaction_Smiles', 'Site', 'hf_activation_energy', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0', 'Electron_density', 'Electrostatic_potential', 'Diamagnetic_shielding', 'Mulliken_charge', 'Bond_index', 'Gross_population', 'HF_C-H_Bond_Length', 'HF_TS_C-C_Bond_Length']
        writer.writerow(header)
        for row in training_set5:
            writer.writerow(row)
