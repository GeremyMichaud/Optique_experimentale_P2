from classe_analyse import Etalonnage, Echantillons


c_1 = 41 * 0.37 / (0.37 + 1 + 1 + 0.634)
c_2 = 41 * 0.732 / (1 + 1 + 0.268 + 0.732)
c_3 = 41 * 1.462 / (1 + 0.462 + 1 + 0.537)
c_4 = 41 * 1.098 / (1 + 0.098 + 1 + 0.907)
c_5 = 41 * 2.195 / (1 + 1 + 0.195 + 0.805)
c_6 = 41
c_33 = 31.7 / 3

# Assurez-vous d'avoir les bons paramètres de domaine des pics et de calibration dans le fichier "classe_analyse.py".
eau = Etalonnage("fichier_trait\\Vodka0_100s.TXT", "Eau", coupe = [280, -1], conc=0, pente=True, zero=True)
ethanol_25 = Etalonnage("fichier_trait\\Ethanol_C25_100s.TXT", "Ethanol à 25 %", conc=25, domaine_pic = [880, 900], coupe = [283, -1], divismin = 2, max_c=400, max_p=600, pente=True, graph=False)
ethanol_50 = Etalonnage("fichier_trait\\Ethanol_C50_100s.TXT", "Ethanol à 50 %", conc=50, domaine_pic = [880, 900], coupe = [280, -1], divismin = 2, max_c=800, max_p=300, pente=True, graph=False)
ethanol_75 = Etalonnage("fichier_trait\\Ethanol_C75_100s.TXT", "Ethanol à 75 %", conc=75, domaine_pic = [880, 900], coupe = [280, -1], divismin = 2, max_c=800, max_p=600, pente=True, graph=False)
ethanol_100 = Etalonnage("fichier_trait\\Ethanol_C100_100s.TXT", "Ethanol à 100 %", conc=100, domaine_pic = [880, 900], coupe = [280, -1], divismin = 2.5, max_c=800, max_p=170, pente=True, graph=False)
Etalonnage.graphique_combiné_ethanolage()



vodka_05 = Echantillons("fichier_trait\\Vodka1_100s.TXT", "Vodka à 5%", d=5, domaine_pic = [880, 900], coupe = [280, -1], max_c=70, max_p=50, pente=True, graph=False)
vodka_10 = Echantillons("fichier_trait\\Vodka2_100s.TXT", "Vodka à 10%", d=5, domaine_pic = [880, 900], coupe = [280, -1], max_c=100, max_p=50, pente=True, graph=True)
vodka_20 = Echantillons("fichier_trait\\Vodka3_100s.TXT", "Vodka à 20%", d=5, domaine_pic = [880, 900], coupe = [280, -1], max_c=100, max_p=50, pente=True, graph=False)
Echantillons.analyse()


# # lim=1395, p=800 
# # pic à 1060 [1000, 1100]
# eau = Etalonnage("spectre_eau10min0.TXT", "eau pure", zero=True, coupe = [520, -150], Eau=True, rawdog = False)
# sucre_68 = Etalonnage("spectre_sucre_b6810min0.TXT","sucre 68%", conc=69, d=8, divismin = 1.5, zero=False, max_c=4500, max_p=4000, coupe = [511, -150], pente=True, Eau=True, graph=False, rawdog =False)
# sucre_61 = Etalonnage("spectre_sucre_b6110min0.TXT", "sucre61%", conc=62, d=8, divismin = 1.5, zero=False, max_c=5200, max_p=4200, coupe = [516, -150], pente=True, Eau=True, graph=False, rawdog = False)
# sucre_54 = Etalonnage("spectre_sucre_b5410min0.TXT","sucre54%", conc=55, d=8, divismin = 1.5, zero=False, max_c=1800, max_p=4000, coupe = [512, -150], pente=True, Eau=True, graph=False, rawdog = False)
# sucre_48 = Etalonnage("spectre_sucre_b4810min0.TXT", "sucre48%", conc=48, d=8, divismin = 1.5, zero=False, max_c=1795, max_p=3800, coupe = [512, -150], pente=True, Eau=True, graph=False, rawdog = False)
# sucre_37 = Etalonnage("spectre_sucre_b3710min0.TXT", "sucre37%", conc=35, d=8, divismin = 6, zero=False, max_c=1200, max_p=2200, coupe = [512, -150], pente=True, Eau=True, graph=False, rawdog = False)
# sucre_25 = Etalonnage("spectre_sucre_b2510min0.TXT", "sucre25%", conc=23, d=8, divismin = 3, zero=False, max_c=800, max_p=1200, coupe = [512, -150], pente=True, Eau=True, graph=False, rawdog = False)
# sucre_11 = Etalonnage("spectre_sucre_b1110min0.TXT", "sucre11%", conc=11, d=8, divismin = 3, zero=False, max_c=750, max_p=400, coupe = [511, -150], pente=True, Eau=True, graph=False, rawdog = False)
# sucre_8 = Etalonnage("spectre_sucre_b810min0.TXT", "sucre 8%", conc=7, d=8, divismin = 1.8, zero=False, max_c=500, max_p=300, coupe = [512, -150], pente=True, Eau=True, graph=True, rawdog = False)
# Echantillons.graphique_combiné_ethanolage()

# guru_r = Echantillons("spectre_guru_rouge10min0.TXT", "Guru rouge", d=3, domaine_pic = [1050, 1100], divismin = 1.6, coupe = [518, -150], max_c=4600, max_p=5000, pente=True, Eau=True, graph=True, rawdog = False)
# #print(guru_r.concentration())
# guru_j = Echantillons("spectre_guru_jaune_10min0.TXT", "Guru jaune", d=3, domaine_pic = [1050, 1100], divismin = 3, coupe = [520, -150], max_c=16000, max_p=80200, pente=True, Eau=True, graph=False, rawdog = False)
# #print(guru_j.concentration())
# guru_rose = Echantillons("spectre_guru_rose10min0.TXT", "Guru rose", d=3, domaine_pic = [1050, 1100], divismin = 2, coupe = [550, -150], max_c=2800, max_p=1500, pente=True, Eau=True, graph=False, rawdog = False)
# #print(guru_j.concentration())
# guru_b = Echantillons("spectre_guru_blanche10min0.TXT", "Guru blanche", d=2, domaine_pic = [1050, 1100], divismin = 1.2, coupe = [567, -150], max_c=90000, Eau=True, graph=False, rawdog = False)
# #print(guru_b.concentration())
# guru_bleu = Echantillons("spectre_guru_bleu_10min1.TXT", "Guru bleu", d=2, domaine_pic = [1050, 1100], divismin = 1, coupe = [568, -150], max_c=90000, Eau=True, graph=False, rawdog = False)
# #print(guru_bleu.concentration())
# Echantillons.analyse()


# guru_r = Echantillons("spectre_guru_rouge10min0.TXT", 8, d=3, nom="Guru rouge", domaine_pic = [800, 920], divismin = 1.6, coupe = [520, -150], max_c=4600, max_p=2800, Eau=True, graph=False, rawdog = False)
# #print(guru_r.concentration())
# guru_j = Echantillons("spectre_guru_jaune_10min0.TXT", 3, d=3, nom="Guru jaune", domaine_pic = [800, 920], divismin = 2.93, coupe = [520, -150], max_c=16000, max_p=80200, Eau=True, graph=False, rawdog = False)
# #print(guru_j.concentration())
# guru_rose = Echantillons("spectre_guru_rose10min0.TXT", 3, d=3, nom="Guru rose", domaine_pic = [800, 920], divismin = 2, coupe = [532, -150], max_c=2600, max_p=1500, pente=True, Eau=True, graph=False, rawdog = False)
# #print(guru_j.concentration())
# guru_b = Echantillons("spectre_guru_blanche10min0.TXT", 1.5, d=2, nom="Guru blanche", domaine_pic = [800, 920], divismin = 1.2, coupe = [570, -150], max_c=90000, max_p=1004000, Eau=True, graph=False, rawdog = False)
# #print(guru_b.concentration())
# guru_bleu = Echantillons("spectre_guru_bleu_10min1.TXT", 0.3, d=2, nom="Guru bleu", domaine_pic = [800, 920], divismin = 1, coupe = [570, -150], max_c=90000, max_p=10030700, Eau=True, graph=False, rawdog = False)
# #print(guru_bleu.concentration())
# Echantillons.analyse()
#guru_bleu.graphique_combiné_echantillon()


#, coupe = 3, lim=24000, p=10000
# guru_r = Echantillons("spectre_guru_rouge1h2.TXT", 8, d=4, nom="Guru rouge", domaine_pic = [1050, 1100], divismin = 5, coupe = [523, -150], max_c=600, max_p=1700, pente=True, Eau=True, graph=False, rawdog = False, balance = 5)
# #print(guru_r.concentration())
# guru_j = Echantillons("spectre_guru_jaune50mjn0.TXT", 3, d=4, nom="Guru jaune", domaine_pic = [1050, 1100], divismin = 2, coupe = [550, -150], max_c=9000, max_p=11200, pente=True, Eau=True, graph=False, rawdog = False, balance = 5)
# #print(guru_j.concentration())
# guru_rose = Echantillons("spectre_guru_rose50min0.TXT", 3, d=4, nom="Guru rose", domaine_pic = [1050, 1100], divismin = 4, coupe = [550, -150], max_c=902, max_p=1800, pente=True, Eau=True, graph=False, rawdog = False, balance = 5)
# #print(guru_j.concentration())
# #guru_r.graphique_combiné_echantillon2()
# Echantillons.analyse()
