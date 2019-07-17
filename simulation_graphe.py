#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 12:16:55 2019

@author: willy
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 11:21:06 2019

@author: willy
"""
import time;
import numpy as np;
import pandas as pd;
import itertools as it;

import algorithme_couverture as algo_couv;
import algorithme_correction as algo_corr;
from pathlib import Path    # http://stackoverflow.com/questions/6004073/how-can-i-create-directories-recursively

import genererMatA as geneMatA;
import generations_mesures as mesures;
import fonctions_auxiliaires as fct_aux;

###############################################################################
#                   fonctions annexes ---> debut
###############################################################################
def determiner_aretes_cliques(cliques):
    """
    former un ensemble contenant les aretes des cliques.
    """
    aretes = []
    for Cu in cliques:
        aretes.extend( it.combinations(Cu,2) )
    return aretes;

###############################################################################
#                   fonctions annexes ---> fin
###############################################################################

###############################################################################
#                   calculer DH ou DC ---> debut
###############################################################################
def calculer_distance_hamming(aretes_cliques_LG, aretes_matE_k_alpha):
    """
    calculer la distance de Hamming en fonction des aretes des cliques du LG trouve
    et les aretes de matE_k_alpha.
    """
    froz_aretes_cliques_LG = set(map(frozenset, aretes_cliques_LG))
    froz_aretes_matE_k_alpha = set(map(frozenset, aretes_matE_k_alpha))
    
    set_arcs_diffs = froz_aretes_cliques_LG.union(froz_aretes_matE_k_alpha) \
                        - froz_aretes_cliques_LG.intersection(
                            froz_aretes_matE_k_alpha)
                        
    return len(set_arcs_diffs), set_arcs_diffs;
###############################################################################
#                   calculer DH ou DC ---> fin
###############################################################################

###############################################################################
#                   creation de matE ---> debut
###############################################################################
#def matriceE(dimMat, 
#             nbre_lien, 
#             chemin_matrices, 
#             chemin_datasets, 
#             nbre_ts, 
#             epsilon, 
#             effet_joule, 
#             test):
#    """
#    cette fonction permet de :
#        * generer aleatoirement une matrice matA
#        * generer les mesures du graphe matA en faisant une propagation de flots 
#            descendantes( des sources vers les puits)
#        * generer la matrice matE (qui doit etre un linegraph)
#        
#    IMPORTANT :============> epsilon = 0.75 <===============
#    dimMat : dimension de la matrice matA
#    seuil : valeur par defaut definissant une adjacence entre 2 sommets
#    chemin_datasets = "data/datasets/"
#    
#    nbre_ts = nombre de time series 
#    """
#    #generer reseau de flots (A DETERMINER) avec mesures 
#    df_matA = None
#    if test == "EN TEST":
#        df_matA = pd.read_csv("datas/data_test/matrices/df_matA_generer.csv", \
#                              index_col = "nodes");
#    else:
#        df_matA = geneMatA.genererMatriceA(dimMat, nbre_lien)
#    df_matA.to_csv(chemin_matrices+"df_matA_generer.csv")
#    dico_sommet_arete = mesures.nommage_arcs( df_matA )
#    mesures.create_datasets(df_matA, dico_sommet_arete,
#                            chemin_datasets, nbre_ts, effet_joule ) 
#    
#    #matrice du linegraphe du reseau de flot a determiner
#    listeArcs = fct_aux.liste_arcs(df_matA)
#    matE = mesures.creation_matE(dico_sommet_arete, listeArcs)
#    matE.to_csv(chemin_matrices+"matE_LG.csv")
#    return matE, df_matA, dico_sommet_arete;
###############################################################################
#                   creation de matE ---> fin
###############################################################################

###############################################################################
#                   sauver la distribution de k_erreur ---> debut
###############################################################################
def sauvegarder_execution_k_alpha(path_distr, 
                                 num_graph,
                                 k, alpha, dc, dh, 
                                 nbre_aretes_matE, 
                                 correl_dc_dh, 
                                 start):
    """
    sauvegarder les caracteristiques de l'execution dans un fichier
    """    
    f = open(path_distr \
             +"distribution_moyDistLine_moyHamming_k_"+str(k)+".txt","a")
    f.write(num_graph+";"\
            +str(k)+";"\
            +str(alpha)+";"\
            +str(dc)+";"\
            +str(dh)+";"\
            +str(nbre_aretes_matE)+";"+\
            str( round(correl_dc_dh, 2) )+";"+\
            str( round(time.time() - start, 2) )+"\n")
    f.close();
###############################################################################
#                   sauver la distribution de k_erreur ---> fin
###############################################################################
    
###############################################################################
#                   sauver les parametres de
#                   l'execution de k_erreur,alpha ---> debut
###############################################################################
def sauvegarder_parametres_execution(nom_graphe, k_erreur, alpha, dc, dh, 
                                     runtime, nbre_aretes_matE,
                                     dico_deleted_add_edges, 
                                     dico_correction, 
                                     dico_voisins_etats_supp,
                                     dico_voisins_etats_ajout, 
                                     sommets_1_vois, 
                                     sommets_1_non_vois):
    """
    sauvegarde les parametres de l'execution de k_erreur,alpha
    """
    dico_return = dict();
    dico_return["nom_graphe"] = nom_graphe;
    dico_return['k_erreur'] = k_erreur;
    dico_return['alpha'] = alpha;
    dico_return['dc'] = dc;
    dico_return['dh'] = dh;
    dico_return['runtime'] = runtime;
    dico_return["nbre_aretes_matE"] = nbre_aretes_matE;
    dico_return["nbre_aretes_restantes"] = \
        len(dico_correction["aretes_restantes"]);
    dico_return["aretes_restantes"] = dico_correction["aretes_restantes"];
    dico_return["aretes_supprimees"] = dico_deleted_add_edges["supprimer"];
    dico_return["aretes_ajoutees"] = dico_deleted_add_edges["ajouter"];
    dico_return["nbre_cliques"] = len(dico_correction["C"]);
    dico_return["cliques"] = dico_correction["C"];
    dico_return["nbre_cliques_old"] = len(dico_correction["C_old"]);
    dico_return["cliques_old"] = dico_correction["C_old"];
    dico_return["ordre_noeuds_traites"] =  dico_correction["ordre_noeuds_traites"];
    dico_return["sommets_par_cliqs"] = dico_correction["sommets_par_cliqs"];
    dico_return["sommets_par_cliqs_avec_aretes"] = \
        dico_correction["sommets_par_cliqs_avec_aretes"];
    dico_return["sommets_a_corriger"] = dico_correction["sommets_a_corriger"];
    dico_return["etats_sommets"] = dico_correction["etats_sommets"];
    dico_return["voisins_etats_aretes_supp"] = dico_voisins_etats_supp;
    dico_return["voisins_etats_aretes_ajout"] = dico_voisins_etats_ajout;
    dico_return["sommets_1_voisins_aretes_modifs"] = sommets_1_vois;
    dico_return["sommets_1_non_voisins_aretes_modifs"] = sommets_1_non_vois;

    return dico_return;

###############################################################################
#                   sauver les parametres de
#                   l'execution de k_erreur,alpha ---> fin
###############################################################################
    
###############################################################################
#                   determiner l etat et le voisinage des 
#                       sommets des aretes supprimees/ajoutees ---> debut 
###############################################################################
def determiner_voisins_etats_aretes_modifiees(matE_k_alpha, 
                                              dico_deleted_add_edges,
                                              etats_sommets):
    """
    determiner l etat et le voisinage des sommets 
    des aretes supprimees ou ajoutees.
    """
    dico_vois_ajout = dict();
    dico_vois_supp = dict();
    etats = [0,1,2,3,-1]
    for cle_action, aretes_modifiees in dico_deleted_add_edges.items():
        for arete in aretes_modifiees:
            vois_inits = matE_k_alpha.index[ matE_k_alpha[arete[0]] == 1].tolist()
            vois_fins = matE_k_alpha.index[ matE_k_alpha[arete[1]] == 1].tolist()
            voisins = set(vois_inits).union(vois_fins);
            dico_vois = dict();
            for vois in voisins:
                etat = etats_sommets[vois] if vois in etats_sommets else None;
                dico_vois[vois] = etat;
            dico_vois[arete[0]] = etats_sommets[arete[0]] \
                                    if arete[0] in etats_sommets else None;
            dico_vois[arete[1]] = etats_sommets[arete[1]] \
                                    if arete[1] in etats_sommets else None;
            dico_etats = dict();
            for etat in etats:
                dico_etats["etat"+str(etat)] = len(list(filter(lambda x: x==etat, 
                                                       dico_vois.values())))
            if cle_action == "supprimer":
                dico_vois_supp[arete] = (len(voisins), dico_etats, dico_vois);
            elif cle_action == "ajouter":
                dico_vois_ajout[arete] = (len(voisins), dico_etats, dico_vois);
            else:
                print("cle_action = {} doesn't exist".format(cle_action))
    
    return dico_vois_supp, dico_vois_ajout;
    pass    
###############################################################################
#                   determiner l etat et le voisinage des 
#                       sommets des aretes supprimees/ajoutees ---> fin
###############################################################################

###############################################################################
#              determiner les sommets a -1 voisins aux sommets des aretes modifiees,
#              puis couvert par + de 2 cliques  ---> fin
###############################################################################
def determiner_sommets_1_avec_conditions(etats_sommets,
                                         matE_k_alpha,
                                         dico_sommets_par_cliques,
                                         dico_deleted_add_edges):
    """
    determiner les listes de sommets a -1 :
        - etant voisins ou adjacents aux sommets des aretes modifiees
        - n'etant pas voisins ou adjacents aux sommets des aretes modifiees
    et couvert par plus de 2 cliques.
    """
    etat_1 = -1;
    
    sommets_1_vois, sommets_1_non_vois = list(), list();
    for cle_action, aretes_modifiees in dico_deleted_add_edges.items():
        for arete in aretes_modifiees:
            vois_inits = matE_k_alpha.index[ matE_k_alpha[arete[0]] == 1].tolist()
            vois_fins = matE_k_alpha.index[ matE_k_alpha[arete[1]] == 1].tolist()
            voisins = set(vois_inits).union(vois_fins);
            sommets_1 = [sommet for sommet, val in etats_sommets.items() \
                                  if val == etat_1];
            sommets_1_vois_tmp = set(sommets_1).intersection(voisins);
            sommets_1_non_vois_tmp = set(sommets_1) - voisins;
            for sommet_1_vois in sommets_1_vois_tmp:
                if len(dico_sommets_par_cliques[sommet_1_vois]) > 2:
                    sommets_1_vois.append(
                            (sommet_1_vois, 
                             len(dico_sommets_par_cliques[sommet_1_vois]))
                            )
            for sommet_1_non_vois in sommets_1_non_vois_tmp:
                if len(dico_sommets_par_cliques[sommet_1_non_vois]) > 2:
                    sommets_1_non_vois.append(
                            (sommet_1_non_vois, 
                             len(dico_sommets_par_cliques[sommet_1_non_vois]))
                            )
            
    return sommets_1_vois, sommets_1_non_vois;
                
###############################################################################
#             simulation graphe avec des k=1 aretes supprimees ---> debut
###############################################################################
def simulation_p_correl_k(matE_LG, 
                              matA_GR, 
                              dico_sommet_arete,
                              chemin_matrices,
                              chemin_datasets,
                              mode, 
                              p_correl, 
                              k_erreur, 
                              num_graph,
                              rep_base,
                              dico_parametres_new):
    """
    executer les algos de couverture et de correction sur des graphes 
    dans lesquels on a supprime (p_correl = 1) k aretes alpha (<alpha_max) fois.
    """
    path_distr = rep_base+"/../"+"distribution/";
    path = Path(path_distr); path.mkdir(parents=True, exist_ok=True)
    print("path_distr = {}".format(path_distr))
    
    list_returns = list();
    for alpha in range(dico_parametres_new["alpha_max"]):
        
        num_graph_alpha = num_graph +"_"+str(alpha)
        print("num_graph={}, k = {}, alpha = {} ==> debut".format( 
                num_graph, k_erreur, alpha))
        
        start = time.time();
        
        try:
            matE_k_alpha = None; 
            dico_proba_cases = dict();
        
            matE_k_alpha, \
            dico_deleted_add_edges = \
                            fct_aux.modif_k_cases(
                                matE_LG.copy(), 
                                k_erreur,
                                dico_parametres_new["methode_delete_add_edges"],
                                p_correl)
            dico_proba_cases = fct_aux.ajouter_proba_matE(
                                matE_k_alpha, 
                                dico_deleted_add_edges,
                                dico_parametres_new["loi_stats"], 
                                p_correl,
                                dico_parametres_new["correl_seuil"])
                                
            matE_k_alpha.to_csv(chemin_matrices\
                                +"matE_"+str(k_erreur)+"_"+str(alpha)+".csv")
            
            # algorithme de couverture
            # ajouter les aretes restantes trop nombreuses ==> OK 
            # Verifier algo de couverture ==> OK
            print("1")
            dico_couverture = algo_couv.algo_decomposition_en_cliques(
                                matE_k_alpha, 
                                dico_sommet_arete, 
                                seuil_U=10, 
                                epsilon=0.75,
                                chemin_datasets=chemin_datasets, 
                                chemin_matrices=chemin_matrices,
                                ascendant_1=True, simulation=True, 
                                dico_proba_cases=dico_proba_cases,
                                dico_parametres_new=dico_parametres_new
                                )
            
            dico_couverture["k_erreur"] = k_erreur;
            dico_couverture["C_old"] = dico_couverture["C"].copy();
            dico_couverture["sommets_a_corriger"] = \
                [k for k, v in dico_couverture["etats_sommets"].items() 
                    if v == -1];
            dico_couverture["nbre_sommets_a_corriger"] = \
                                    len(dico_couverture["sommets_a_corriger"]);
            # algorithme correction pour k = 1
            print("2")
            dico_correction = dico_couverture;
            if -1 in dico_couverture['etats_sommets'].values():
                aretes_matE_k_alpha = fct_aux.liste_arcs(matE_k_alpha)
                dico_gamma_noeud = fct_aux.gamma_noeud(
                                        matE_k_alpha, 
                                        aretes_matE_k_alpha)
                dico_correction = algo_corr.correction_cliques(
                                    dico_correction,
                                    aretes_matE_k_alpha,
                                    dico_gamma_noeud,
                                    dico_parametres_new)
            elif -1 not in dico_couverture['etats_sommets'].values() and \
                len(dico_correction['aretes_restantes']) > 0:
                dico_correction["C"] = \
                            dico_correction["C"] \
                            + dico_correction['aretes_restantes'];
                dico_correction["sommets_par_cliqs_avec_aretes"] = \
                fct_aux.couverture_par_sommets(
                    sommets_matE = list(dico_correction["etats_sommets"].keys()),
                    C = dico_correction["C"]);
            elif -1 not in dico_couverture['etats_sommets'].values() and \
                len(dico_correction['aretes_restantes']) == 0:
                dico_correction["sommets_par_cliqs_avec_aretes"] = \
                fct_aux.couverture_par_sommets(
                    sommets_matE = list(dico_correction["etats_sommets"].keys()),
                    C = dico_correction["C"]);
                        
            # calcul DH et DC
            print("3")
            aretes_cliques = determiner_aretes_cliques(dico_correction["C"]);
            aretes_matE_k_alpha = fct_aux.liste_arcs(matE_k_alpha);
            aretes_matE_LG = fct_aux.liste_arcs(matE_LG);
            
            print("4")
            dc, set_dc = calculer_distance_hamming(
                                aretes_cliques, 
                                aretes_matE_k_alpha);
            dh, set_dh = calculer_distance_hamming(
                                aretes_cliques, 
                                aretes_matE_LG);
            X1 = abs(dc - k_erreur);
            correl_dc_dh = abs(dh - X1)/(k_erreur + dc) if k_erreur+dc != 0 else -1;
            print("k={}, ".format(k_erreur)\
                  +"moy_dc={}, ".format(dc) \
                  +"moy_dh={}, ".format(dh) \
                  +"moy_dc-k={}, ".format(X1) \
                  +"moy_dc+k={}, ".format( k_erreur+dc ) \
                  +"corr={} ".format( round(correl_dc_dh,2) )
                  )
            
            # sauvegarde dans un fichier distribution
            sauvegarder_execution_k_alpha(path_distr, 
                                 num_graph,
                                 k_erreur, alpha, dc, dh, 
                                 len(aretes_matE_LG), 
                                 correl_dc_dh, 
                                 start)
            print("5")
            
            # voisins des sommets des aretes supprimees/ajoutees et leurs etats
            dico_voisins_etats_supp,\
            dico_voisins_etats_ajout = \
                                determiner_voisins_etats_aretes_modifiees(
                                    matE_k_alpha, 
                                    dico_deleted_add_edges,
                                    dico_couverture["etats_sommets"]
                                )
            print("6")
            # sommets a -1 etant les voisins des sommets de l arete modifies
            sommets_1_vois, \
            sommets_1_non_vois = determiner_sommets_1_avec_conditions(
                                dico_correction["etats_sommets"],
                                matE_k_alpha,
                                dico_correction["sommets_par_cliqs_avec_aretes"],
                                dico_deleted_add_edges
                                );
            print("7")
            
            # sauvegarde les parametres de l'execution
            num_graph_alpha = num_graph +"_"+str(alpha)
            list_returns.append(sauvegarder_parametres_execution(
                                            num_graph_alpha, k_erreur, alpha, dc, dh, 
                                            time.time() - start, 
                                            len(aretes_matE_LG), 
                                            dico_deleted_add_edges, 
                                            dico_correction, 
                                            dico_voisins_etats_supp,
                                            dico_voisins_etats_ajout, 
                                            sommets_1_vois, 
                                            sommets_1_non_vois))
            print("8")
        
        except  Exception as e :
            print("####### EmptyDataError {}".format(dico_parametres_new["coef_fct_cout"][2]) \
                  +" {},".format(dico_parametres_new["mode_select_noeuds_1"]) \
                  +" seuil={}".format(dico_parametres_new["correl_seuil"]) \
                  +" p={}".format(p_correl) \
                  +" k={}".format(k_erreur) \
                  +" num_graph={}".format(num_graph) \
                  +" alpha={}".format(alpha) \
                  +" e={}".format(e)
                  )
            
    return list_returns;
    pass
###############################################################################
#             simulation graphe avec des k=1 aretes supprimees ---> fin
###############################################################################