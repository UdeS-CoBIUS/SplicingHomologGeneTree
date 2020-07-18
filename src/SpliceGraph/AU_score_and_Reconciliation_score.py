#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
from ete3 import Tree
from ete3 import PhyloTree
import numpy as np 
import glob
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
figure(num=None, figsize=(30, 6), dpi=80, facecolor='w', edgecolor='k')

def numberOfDuplicationGtoS(tree):
    ntrees, ndups, sptrees = tree.get_speciation_trees()
    return ndups

def numberOfLostGtoS(tree):
    compt = 0;
    for node in tree.traverse("preorder"):
        if(node.name == "L"):
            compt += 1
    return compt

def splided_gene_family(id):	
	sourceToTarget = "initialSource/" + id + "_initialsource2target.txt"
	file = open(sourceToTarget, "r")
	dict_ = {}	
	lines = file.readlines()
	
	flag = False
	for l in lines:
		l = l.replace("\n", "")
		p = l.split(" ")
		g_id = p[1]
		c_id = p[0]
		if g_id in dict_.keys():
			dict_[g_id] +=1
			flag = True
		else:
			dict_[g_id] =1

	return flag

def compute_ratio_consel_test():
	for seq_used in ["_Aln_ByOrthoGroup.phy._seq_OrthoGroup_consel", "_Aln_ByEnsemblGeneTree.phy._seq_Ensembl_consel", "_Aln_ByIsoSelGeneTree.phy._seq_IsoSel_consel"]:
		listes = [[0,1], [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]
		au_sht_list = []
		au_lt_list = []
		au_ltdb_list = []
		au_isosel_list = []

		for e in  listes:
			weightStructure = e[0]
			weightSequence = e[1]
			sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))	
			au_sht_pass =0
			au_lt_pass =0
			au_ltdb_pass =0
			au_isosel_pass =0				
			total = 0
			for x in glob.glob("ensemblGeneTree/sequences_" + sub_path + "/*" + seq_used): 
				id = x.split("/")[-1].split("_")[0]
				if not(splided_gene_family(id)):
					continue
				consel_result = open(x, "r")
				lines = consel_result.readlines()
				
				if len(lines) == 9:
					au_sht = float(lines[4].split()[4])
					au_lt = float(lines[5].split()[4])
					au_ltdb = float(lines[6].split()[4])
					au_isosel = float(lines[7].split()[4])


					if au_sht >= 0.05:
						au_sht_pass += 1.0
					
					if au_lt >= 0.05:
						au_lt_pass += 1.0
					
					if au_ltdb >= 0.05:						
						au_ltdb_pass += 1.0
					
					if au_isosel >= 0.05:
						au_isosel_pass += 1.0
					total += 1
			
			au_sht_list.append(int(round(100.0*au_sht_pass/total, 0)))
			au_lt_list.append(int(round(100.0*au_lt_pass/total, 0)))
			au_ltdb_list.append(int(round(100.0*au_ltdb_pass/total, 0)))
			au_isosel_list.append(int(round(100.0*au_isosel_pass/total, 0)))

		barWidth = 0.2
		r1 = np.arange(len(au_ltdb_list))
		r2 = [x + barWidth for x in r1]
		r3 = [x + 2*barWidth for x in r1]
		r4 = [x + 3*barWidth for x in r1]

		
		# Create blue bars
		rects1 = plt.bar(r1, au_ltdb_list, width = barWidth, color = '#4285F4', edgecolor = '#4285F4', capsize=30, label='LT_db')

		# Create cyan bars
		rects2 = plt.bar(r2, au_lt_list, width = barWidth, color = '#EA4335', edgecolor = '#EA4335', capsize=30, label='LT')

		# Create cyan bars
		rects3 = plt.bar(r3, au_sht_list, width = barWidth, color = '#FBBC03', edgecolor = '#FBBC03', capsize=30, label='SHT')

		rects4 = plt.bar(r4, au_isosel_list, width = barWidth, color = '#6b392b', edgecolor = '#6b392b', capsize=30, label='IsoSel')


		# general layout
		plt.xticks([r + 2*barWidth for r in range(len(au_ltdb_list))], [r'$\alpha=0.0$', r'$\alpha=0.2$', r'$\alpha=0.4$', r'$\alpha=0.5$', r'$\alpha=0.6$', r'$\alpha=0.8$', r'$\alpha=1$'])
		plt.xticks(fontsize=80)
		plt.yticks(fontsize=60)
		plt.legend(prop={'size': 35})
		plt.legend(loc="best", borderaxespad=0., prop={'size': 35}, ncol=4)
		#plt.xlabel('xlabel', fontsize=22)
		plt.ylabel('Percentage of trees that which pass the AU test', fontsize=45)
		#plt.savefig("RI.png")
		# Show graphic
		def autolabel(rects):
		    """
		    Attach a text label above each bar displaying its height
		    """
		    for rect in rects:
		        height = rect.get_height()        
		        plt.text(rect.get_x() + rect.get_width()/2., 0.99*height,
		                 str(height)+  "%",
		                ha='center', va='bottom',  fontsize=30)
		plt.ylim(0, 100)
		autolabel(rects1)                
		autolabel(rects2)                
		autolabel(rects3)                
		autolabel(rects4)
		plt.show()		


	
def computeReconciliation():
	#1: LT, 2: :LT_db, 3:Isosel	
	for method in [1, 2, 3]:
		listes = [[0,1], [0.2,0.8], [0.4,0.6], [0.5,0.5], [0.6,0.4], [0.8,0.2], [1,0]]
		equal_trees_list = []
		different_tree_same_rconciliation_list = []
		sht_better_list = []
		other_better_list = []

		for e in  listes:
			weightStructure = e[0]
			weightSequence = e[1]
			sub_path = str(int(weightStructure*100)) + "_" + str(int(weightSequence*100))	
			equal_trees = 0
			different_tree_same_rconciliation = 0
			sht_better = 0
			other_better = 0
			total = 0
			for x in glob.glob("ensemblGeneTree/trees_" + sub_path + "/*"): 
				try:
					geneFamilyId = x.split("/")[-1].split("_")[0].split(".")[0]
					if not(splided_gene_family(geneFamilyId)):
						continue
					trees = open(x, "r")
					lines = trees.readlines()
					if len(lines) == 4:
						sht_tree_nw = lines[0]
						sht_tree = PhyloTree(sht_tree_nw.replace("\n", ""))

						other_tree_nw = lines[method]
						other_tree = PhyloTree(other_tree_nw.replace("\n", ""))
						prot_gene_spe = open("ensemblGeneTree/transcript_protein_gene/" + geneFamilyId + "_update.txt", "r")
						sptree = PhyloTree("(((Drosophila melanogaster:0.155009,Caenorhabditis elegans strain N2:0.188311)Ecdysozoa:0.0001,((Ciona intestinalis:0.151327,Ciona savignyi:0.154463)Ciona:0.0329096,((Eptatretus burgeri:0.128455,Petromyzon marinus:0.169105)Cyclostomata:0.000769775,((((((((((((Astyanax mexicanus Astyanax_mexicanus-2.0:0.00393522,Astyanax mexicanus Astyanax_mexicanus-1.0.2:0.00769478)Astyanax mexicanus:0.0783682,Pygocentrus nattereri:0.0888718)Characoidei:0.0164521,Electrophorus electricus:0.107985)Characiphysae:0.0045476,Ictalurus punctatus:0.100737)Characiphysae:0.00945255,Danio rerio:0.110624)Otophysi:0.0001,(Clupea harengus:0.102439,Denticeps clupeoides:0.119621)Clupeiformes:0.00740506)Otomorpha:0.00259501,(((((((((((Fundulus heteroclitus:0.0802258,Cyprinodon variegatus:0.0839142)Cyprinodontoidei:0.0001,(((Xiphophorus maculatus:0.00748124,Xiphophorus couchianus:0.0110388)Xiphophorus:0.0241812,Gambusia affinis:0.0308088)Poeciliinae:0.00988622,(((Poecilia latipinna:0.00440636,Poecilia formosa:0.00521364)Poecilia:0.00197579,Poecilia mexicana:0.00794921)Poecilia:0.0164401,Poecilia reticulata:0.0285797)Poecilia:0.0101232)Poeciliinae:0.043166)Cyprinodontoidei:0.0137751,Kryptolebias marmoratus:0.0930248)Cyprinodontiformes:0.012096,(((Oryzias latipes ASM223467v1:0.0148451,Oryzias latipes ASM223471v1:0.0165149)Oryzias latipes:0.00101579,Oryzias latipes ASM223469v1:0.0171142)Oryzias latipes:0.0445773,Oryzias melastigma:0.0589407)Oryzias:0.0512964)Atherinomorphae:0.0101824,Gouania willdenowi:0.107449)Ovalentaria:0.0001,(((((Amphiprion percula:0.00551478,Amphiprion ocellaris:0.00676522)Amphiprion:0.0269473,Acanthochromis polyacanthus:0.0400077)Pomacentridae:0.0212392,Stegastes partitus:0.05178)Pomacentridae:0.0186711,Parambassis ranga:0.0783139)Ovalentaria incertae sedis:0.00775579,((((((Maylandia zebra:0.00223424,Astatotilapia calliptera:0.00236576)Haplochromini:0.000578973,Pundamilia nyererei:0.00881603)Haplochromini:0.000671,Haplochromis burtoni:0.00745101)Haplochromini:0.00816067,Neolamprologus brichardi:0.0156223)Pseudocrenilabrinae:0.0105045,Oreochromis niloticus:0.0230369)Pseudocrenilabrinae:0.0381004,Amphilophus citrinellus:0.0554442)Cichlidae:0.0208118)Ovalentaria:0.0185306)Ovalentaria:0.0001,((((Tetraodon nigroviridis:0.067489,Takifugu rubripes:0.074401)Tetraodontidae:0.0352102,Mola mola:0.0977148)Tetraodontiformes:0.0018189,(((Cottoperca gobio:0.0749839,Gasterosteus aculeatus:0.0930361)Perciformes:0.00120023,Larimichthys crocea:0.0742248)Eupercaria:0.00681344,Labrus bergylta:0.0840581)Eupercaria:0.00662427)Eupercaria:0.0001,(((Cynoglossus semilaevis:0.095666,Scophthalmus maximus:0.100244)Pleuronectoidei:0.0001,((Seriola dumerili:0.0127718,Seriola lalandi dorsalis:0.0127882)Seriola:0.0431147,Lates calcarifer:0.0602003)Carangaria:0.0254459)Carangaria:0.0001,((Mastacembelus armatus:0.0747997,Monopterus albus:0.0860603)Synbranchiformes:0.00518692,(Anabas testudineus:0.0721382,Betta splendens:0.0850318)Anabantoidei:0.00851808)Anabantaria:0.00727071)Percomorphaceae:0.00517341)Percomorphaceae:0.00920235)Percomorphaceae:0.0164438,Hippocampus comes:0.11414)Percomorphaceae:0.00159408,Periophthalmus magnuspinnatus:0.114372)Percomorphaceae:0.00367202,Gadus morhua:0.116188)Acanthomorphata:0.00378108,(Hucho hucho:0.0868644,Esox lucius:0.0989656)Protacanthopterygii:0.032928)Euteleosteomorpha:0.0066721)Clupeocephala:0.0001,(Scleropages formosus:0.116167,Paramormyrops kingsleyae:0.119273)Osteoglossiformes:0.0150268)Osteoglossocephalai:0.00643571,Lepisosteus oculatus:0.156606)Neopterygii:0.00646505,Erpetoichthys calabaricus:0.13214)Actinopterygii:0.0001,(((((((((Erinaceus europaeus:0.101026,Sorex araneus:0.107684)Eulipotyphla:0.000657858,((((Pteropus vampyrus:0.0800211,Myotis lucifugus:0.0826389)Chiroptera:0.00170224,(Equus caballus:0.00638598,Equus asinus asinus:0.00665402)Equus:0.0757903)Laurasiatheria:0.00292536,(((((Ursus maritimus:0.00323691,Ursus americanus:0.00365309)Ursus:0.0138412,Ailuropoda melanoleuca:0.0173938)Ursidae:0.0376385,(Mustela putorius furo:0.0160511,Neovison vison:0.0160589)Mustelinae:0.0387644)Caniformia:0.00639456,((Canis lupus dingo:0.000957111,Canis lupus familiaris:0.00135289)Canis lupus:0.0104553,Vulpes vulpes:0.0120047)Canidae:0.0473409)Caniformia:0.00610117,((Panthera pardus:0.00387187,Panthera tigris altaica:0.00451813)Panthera:0.00754061,Felis catus:0.0116044)Felidae:0.0534155)Carnivora:0.0143547)Laurasiatheria:0.00315901,(((((Ovis aries:0.0110756,Capra hircus:0.0112644)Caprinae:0.0170002,((((Bos indicus x Bos taurus UOA_Brahman_1:0.00247744,Bos indicus x Bos taurus UOA_Angus_1:0.00264256)Bos indicus x Bos taurus:0.0001,Bos taurus:0.00186489)Bos:0.00292802,Bos mutus:0.00540361)Bos:0.0001,Bison bison bison:0.00517296)Bovinae:0.0223216)Bovidae:0.0438906,Tursiops truncatus:0.0658466)Cetartiodactyla:0.0100522,Vicugna pacos:0.0818253)Cetartiodactyla:0.00120462,((((((Sus scrofa strain jinhua:0.00184196,Sus scrofa strain meishan:0.00198804)Sus scrofa:0.000258996,Sus scrofa strain rongchang:0.002011)Sus scrofa:0.000304281,Sus scrofa strain tibetan:0.00244605)Sus scrofa:0.0001,Sus scrofa strain bamei:0.0021991)Sus scrofa:0.000107332,Sus scrofa strain wuzhishan:0.00276613)Sus scrofa:0.000516651,((Sus scrofa strain reference:0.00167198,Sus scrofa strain usmarc:0.00196802)Sus scrofa:0.0001,((Sus scrofa strain largewhite:0.00133982,Sus scrofa strain pietrain:0.00150018)Sus scrofa:3.39977e-05,((Sus scrofa strain hampshire:0.00123984,Sus scrofa strain berkshire:0.00154016)Sus scrofa:0.000120658,Sus scrofa strain landrace:0.00145934)Sus scrofa:4.41161e-05)Sus scrofa:0.000169627)Sus scrofa:0.00096732)Sus scrofa:0.0733376)Cetartiodactyla:0.00883307)Laurasiatheria:0.0155472)Laurasiatheria:0.0001,(((((((((Cavia porcellus:0.00818682,Cavia aperea:0.0125132)Cavia:0.0661817,Chinchilla lanigera:0.0705433)Hystricomorpha:0.00414421,Octodon degus:0.0759569)Hystricomorpha:0.00138677,((Heterocephalus glaber HetGla_female_1.0:0.000503489,Heterocephalus glaber HetGla_1.0:0.000766511)Heterocephalus glaber:0.0562146,Fukomys damarensis:0.0577254)Bathyergidae:0.0193335)Hystricomorpha:0.0207161,(((Urocitellus parryii:0.0121268,Spermophilus dauricus:0.0130632)Marmotini:0.000652288,Ictidomys tridecemlineatus:0.0138277)Marmotini:0.00565301,Marmota marmota marmota:0.0183762)Marmotini:0.0711477)Rodentia:0.00260726,(Dipodomys ordii:0.0915058,Castor canadensis:0.0926542)Castorimorpha:0.00511557)Rodentia:0.0044061,(((((((((Mus spicilegus:0.0085474,Mus spretus strain SPRET/EiJ:0.0089626)Mus:0.00109964,(((((((Mus musculus strain reference CL57BL6:0.0001,Mus musculus strain C57BL/6NJ:0.000756133)Mus musculus:0.000884866,Mus musculus strain NZO/HlLtJ:0.00152513)Mus musculus:0.0001,((((Mus musculus strain A/J:0.000329467,Mus musculus strain BALB/cJ:0.000950533)Mus musculus:6.02928e-05,((Mus musculus strain C3H/HeJ:0.000397156,Mus musculus strain CBA/J:0.000522844)Mus musculus:0.000220714,Mus musculus strain DBA/2J:0.00115929)Mus musculus:0.000316136)Mus musculus:0.00033581,Mus musculus strain AKR/J:0.000688249)Mus musculus:0.0001,(Mus musculus strain FVB/NJ:0.000912933,Mus musculus strain NOD/ShiLtJ:0.00148707)Mus musculus:0.000239914)Mus musculus:0.000368911)Mus musculus:0.0001,(Mus musculus strain 129S1/SvImJ:0.000465156,Mus musculus strain LP/J:0.000604844)Mus musculus:0.000970521)Mus musculus:0.000412647,Mus musculus domesticus strain WSB/EiJ:0.00194964)Mus musculus:0.00280959,Mus musculus musculus strain PWK/PhJ:0.00538439)Mus musculus:0.000400067,Mus musculus castaneus strain CAST/EiJ:0.00513373)Mus musculus:0.00549987)Mus:0.0109037,Mus caroli strain CAROLI_EIJ:0.0204399)Mus:0.0168539,Mus pahari strain PAHARI_EIJ:0.0365829)Mus:0.0231643,Rattus norvegicus:0.0629608)Murinae:0.0162398,Meriones unguiculatus:0.0763181)Muridae:0.000636067,(((((Cricetulus griseus CHOK1GS_HDv1:0.000317,Cricetulus griseus CriGri_1.0:0.000693)Cricetulus griseus:0.000777232,Cricetulus griseus CriGri-PICR:0.000667768)Cricetulus griseus:0.0478828,Mesocricetus auratus:0.0525298)Cricetinae:0.0160348,Peromyscus maniculatus bairdii:0.0637651)Cricetidae:0.00356177,Microtus ochrogaster:0.0703313)Cricetidae:0.0104833)Muroidea:0.0120886,Nannospalax galili:0.0881937)Muroidea:0.00981086,Jaculus jaculus:0.0948615)Myomorpha:0.007332)Rodentia:0.00367649,(Oryctolagus cuniculus:0.0778675,Ochotona princeps:0.0908825)Lagomorpha:0.0251812)Glires:0.0001,((((((((((Pan troglodytes:0.002216,Pan paniscus:0.003324)Pan:0.00429938,Homo sapiens:0.00661063)Homininae:0.00183587,Gorilla gorilla:0.00859768)Homininae:0.00840489,Pongo abelii:0.0171815)Hominidae:0.0027716,Nomascus leucogenys:0.0196173)Hominoidea:0.011061,(((Piliocolobus tephrosceles:0.00970118,Colobus angolensis palliatus:0.0111088)Colobinae:0.00169751,(Rhinopithecus roxellana:0.00209587,Rhinopithecus bieti:0.00301413)Rhinopithecus:0.00966999)Colobinae:0.00484466,((((Cercocebus atys:0.00568422,Mandrillus leucophaeus:0.00676578)Cercopithecinae:0.000535902,(Papio anubis:0.00407742,Theropithecus gelada:0.00410258)Cercopithecinae:0.0024316)Cercopithecinae:0.00109526,((Macaca fascicularis:0.00205989,Macaca mulatta:0.00220011)Macaca:0.000877645,Macaca nemestrina:0.00415235)Macaca:0.00437927)Cercopithecinae:0.00407974,Chlorocebus sabaeus:0.0116629)Cercopithecinae:0.00499439)Cercopithecidae:0.0135309)Catarrhini:0.0177941,(((Cebus capucinus imitator:0.0260843,Saimiri boliviensis boliviensis:0.0262457)Cebidae:0.00269488,Callithrix jacchus:0.0300851)Cebidae:0.0001,Aotus nancymaae:0.0252368)Platyrrhini:0.0237389)Simiiformes:0.0275616,Carlito syrichta:0.0797473)Haplorrhini:0.000100626,(((Prolemur simus:0.0371796,Propithecus coquereli:0.0381704)Lemuriformes:0.00513981,Microcebus murinus:0.0414252)Lemuriformes:0.0275709,Otolemur garnettii:0.0747342)Strepsirrhini:0.00703666)Primates:0.0158798,Tupaia belangeri:0.0929595)Euarchontoglires:0.00575616)Euarchontoglires:0.0001)Boreoeutheria:0.00121202,(((Loxodonta africana:0.0623719,Procavia capensis:0.0761881)Afrotheria:0.0167982,Echinops telfairi:0.0966318)Afrotheria:0.019805,(Dasypus novemcinctus:0.0777143,Choloepus hoffmanni:0.0782357)Xenarthra:0.0226405)Eutheria:0.00463665)Eutheria:0.0132897,((((Vombatus ursinus:0.0353722,Phascolarctos cinereus:0.0375378)Diprotodontia:0.0214841,Notamacropus eugenii:0.0630759)Diprotodontia:0.0122477,Sarcophilus harrisii:0.0663409)Marsupialia:0.00209439,Monodelphis domestica:0.0747086)Marsupialia:0.0424224)Theria:0.00745259,Ornithorhynchus anatinus:0.127283)Mammalia:0.00148067,(((((Pogona vitticeps:0.104868,Anolis carolinensis:0.105202)Iguania:0.0134654,Notechis scutatus:0.111395)Toxicofera:0.0001,Salvator merianae:0.115591)Episquamata:0.0122703,Sphenodon punctatus:0.127363)Lepidosauria:0.000969327,((((((((Gallus gallus:0.0369694,Meleagris gallopavo:0.0430606)Phasianidae:0.00645001,Coturnix japonica:0.043595)Phasianidae:0.00126886,Numida meleagris:0.0489561)Galliformes:0.0297436,(Anas platyrhynchos platyrhynchos:0.0284988,Anser brachyrhynchus:0.0314312)Anatidae:0.0385057)Galloanserae:0.00751628,(((Calidris pugnax:0.016933,Calidris pygmaea:0.018207)Calidris:0.0469876,Melopsittacus undulatus:0.0749724)Neognathae:0.0086076,((((Cyanistes caeruleus:0.014851,Parus major:0.015359)Paridae:0.030815,Ficedula albicollis:0.04559)Passeriformes:0.000514931,(((Lonchura striata domestica:0.0200049,Taeniopygia guttata:0.0208851)Estrildinae:0.0209595,Serinus canaria:0.0412755)Passeroidea:0.0001,(Zonotrichia albicollis:0.0116358,Junco hyemalis:0.0129942)Passerellidae:0.0283156)Passeriformes:0.00860953)Passeriformes:0.0173337,(Lepidothrix coronata:0.0119296,Manacus vitellinus:0.0130504)Pipridae:0.0536023)Passeriformes:0.0105191)Neognathae:0.0112947)Neognathae:0.00150368,((((Apteryx owenii:0.00131082,Apteryx haastii:0.00213918)Apteryx:0.00323517,Apteryx rowi:0.00448483)Apteryx:0.0331671,Dromaius novaehollandiae:0.0391011)Palaeognathae:0.0308303,Nothoprocta perdicaria:0.0671553)Palaeognathae:0.0167124)Aves:0.0215438,Crocodylus porosus:0.112377)Archosauria:0.00203969,(((Chelonoidis abingdonii:0.0225391,Gopherus agassizii:0.0228709)Testudinidae:0.0125796,Chrysemys picta bellii:0.0340404)Testudinoidea:0.0413783,Pelodiscus sinensis:0.0782918)Cryptodira:0.0404071)Archelosauria:0.0172966)Sauria:0.00707042)Amniota:0.0125733,Xenopus tropicalis:0.145985)Tetrapoda:0.0001,Latimeria chalumnae:0.11922)Sarcopterygii:0.017797)Euteleostomi:0.0001,Callorhinchus milii:0.14598)Gnathostomata:0.00454606)Vertebrata:0.0418931)Chordata:0.00968241)Bilateria:0.120573,Saccharomyces cerevisiae strain S288C:0.120573)ENSEMBLTREE:1;", format=1)
						for n in sptree.traverse("postorder"):
							if n.is_leaf():
								n.name = n.name.replace(" ", "")
								n.name = n.name.replace("/", "")
								n.name = n.name.replace("_", "")        
							else:
								n.name = n.name.replace(" ", "")
								n.name = n.name.replace("/", "")
								n.name = n.name.replace("_", "")        
						dictSpe = {}
						dictSpeInv = {}
						compt = 10
						for node in sptree.traverse("preorder"):
							if node.is_leaf():

								var = "A" + str(compt)
								compt = compt + 1
								dictSpe[node.name] = var
								dictSpeInv[var] = node.name
								node.name = var

						lines = prot_gene_spe.readlines()
						mapping_gene_species = {}
						for line in lines:
							line = line.replace("\n", "")
							parts = line.split("\t")
							species_name = parts[2]
							species_name = species_name.replace(" ", "")
							species_name = species_name.replace("/", "")
							species_name = species_name.replace("_", "")        			
							mapping_gene_species[parts[1]] = species_name

						rf, max_rf, common_leaves, parts_t1, parts_t2, v5, v6 = sht_tree.robinson_foulds(other_tree)
						if rf == 0:
							equal_trees += 1
						else:									
							
							dictTreeBest = {}
							dictOthertree = {}

							compt = 10
							for node in sht_tree.traverse("preorder"):
								if node.is_leaf():
									var = "B" + str(compt)
									compt = compt + 1
									dictTreeBest[var] = node.name
									node.name = var
							
								
							compt = 10
							for node in other_tree.traverse("preorder"):
								if node.is_leaf():
									var = "D" + str(compt)
									compt = compt + 1
									dictOthertree[var] = node.name
									node.name = var

							
							for node in sht_tree.traverse("postorder"):
								if node.is_leaf():					
									node.name =   dictSpe[mapping_gene_species[dictTreeBest[node.name]]] + "_" + node.name
							
							recon_tree1, events1 = sht_tree.reconcile(sptree)
							nbDup1 = numberOfDuplicationGtoS(recon_tree1)
							nbLost1 = numberOfLostGtoS(recon_tree1)
							reconciliation_score1 = nbDup1 + nbLost1


							for node in other_tree.traverse("postorder"):
								if node.is_leaf():							
									node.name = dictSpe[mapping_gene_species[dictOthertree[node.name]]]  + "_" + node.name

							recon_tree2, events2 = other_tree.reconcile(sptree)
							nbDup2 = numberOfDuplicationGtoS(recon_tree2)
							nbLost2 = numberOfLostGtoS(recon_tree2)
							reconciliation_score2 = nbDup2 + nbLost2

							if reconciliation_score1 == reconciliation_score2:
								different_tree_same_rconciliation += 1

							if reconciliation_score1 < reconciliation_score2:
								sht_better += 1

							if reconciliation_score2 < reconciliation_score1:
								other_better += 1

						total += 1
				except:
					pass
			equal_trees_list.append(int(round(100.0*equal_trees/total, 0)))
			different_tree_same_rconciliation_list.append(int(round(100.0*different_tree_same_rconciliation/total, 0)))
			sht_better_list.append(int(round(100.0*sht_better/total, 0)))
			other_better_list.append(int(round(100.0*other_better/total, 0)))

		barWidth = 0.2
		# The x position of bars
		r1 = [x - 0.1 for x in np.arange(len(equal_trees_list))]
		#r1 = np.arange(len(bars1))
		r2 = [x + barWidth for x in r1]
		r3 = [x + 2*barWidth for x in r1]
		r4 = [x + 3*barWidth for x in r1]

		if method == 1:
			other_label = "LT lowest reconciliation cost"
			color = "#EA4335"
		elif method == 2:
			other_label = "LT_db lowest reconciliation cost"
			color = "#4285F4"
		elif method == 3:
			other_label = "IsoSel lowest reconciliation cost"
			color = "#6b392b"

		# Create blue bars
		rects1 = plt.bar(r1, equal_trees_list, width = barWidth, color = '#5E9252', edgecolor = '#5E9252', capsize=30, label='Same tree')

		# Create cyan bars
		rects2 = plt.bar(r2, different_tree_same_rconciliation_list, width = barWidth, color = '#999999', edgecolor = '#999999', capsize=30, label='Different tree, same reconciliation cost')

		# Create cyan bars
		rects3 = plt.bar(r3, other_better_list, width = barWidth, color = color, edgecolor = color, capsize=30, label=other_label)

		rects4 = plt.bar(r4, sht_better_list, width = barWidth, color = '#FBBC03', edgecolor = '#FBBC03', capsize=30, label='SHT lowest reconciliation cost')
		# general layout
		# general layout
		plt.xticks([r + barWidth for r in range(len(equal_trees_list))], [r'$\alpha=0.0$', r'$\alpha=0.2$', r'$\alpha=0.4$', r'$\alpha=0.5$', r'$\alpha=0.6$', r'$\alpha=0.8$', r'$\alpha=1$'])  # will be used to label x-ticks)
		plt.xticks(fontsize=70)
		plt.yticks(fontsize=70)
		plt.legend(prop={'size': 35})
		plt.legend(loc="best", borderaxespad=0., prop={'size': 35}, ncol=2)

		#  plt.xlabel('Ratio tree sizes LT/LT_BD', fontsize=40)
		plt.ylabel('Percentage of gene families', fontsize=70)
		#plt.savefig("RI.png")
		# Show graphic


		def autolabel(rects, x):
		    """
		    Attach a text label above each bar displaying its height
		    """
		    for rect in rects:
		        height = rect.get_height()        
		        plt.text(rect.get_x() + rect.get_width()/2. , height + x,
		                 str(height)+  "%",
		                ha='center', va='bottom',  fontsize=30)
		plt.ylim(0, 100)
		autolabel(rects1, 0)                
		autolabel(rects2, 0)                
		autolabel(rects3, -0.5)     
		autolabel(rects4, 2.5)                
		plt.show()





#compute_ratio_consel_test()
computeReconciliation()
