#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 16:06:00 2019

@author: willy
"""
import random;

import numpy as np;
###############################################################################
#                => debut
###############################################################################

###############################################################################
#                => fin
###############################################################################

###############################################################################
#                => debut
###############################################################################

###############################################################################
#                => fin
###############################################################################

###############################################################################
#               calcul de la compression d un sommet => debut
###############################################################################
def compression_sommet(id_sommet_1,
                       nom_sommet_1,
                       noms_sommets_1,
#                       cliques_sommet_1,
                       cliques_par_nom_sommets,
                       cliques_couvertures,
                       aretes_LG_k_alpha_cor,
                       sommets_LG,
                       mode_correction,
                       critere_correction,
                       number_items_pi1_pi2,
                       DBG
                       
                       id_sommet_1,
                       sommet_1,
                       sommets_1,
                       cliques_par_sommets,
                       cliques_couv_cor = cliques_couvertures_cor,
                            aretes_LG_k_alpha_cor = aretes_LG_k_alpha_cor,
                            mode_select_noeuds_1 = mode_select_noeuds_1,
                            critere_select_pi1_pi2 = critere_select_pi1_pi2,
                            number_items_pi1_pi2 = number_items_pi1_pi2,
                            DBG = DBG):
    """ retourne la compression d'un sommet sommet_z. 
    
    la compression est le triplet (pi1, pi2, ps) dans lequel 
        * pi1, pi2 sont des cliques qui fusionnent 
            - des cliques augmentantes C1, C2 ou 
            - des cliques contractables C1, C2 ou 
            - un ensemble S1 tel que S1 n'est contenu par aucune clique C1 ou C2
        * pi1, pi2 sont des augmentations
        * ps est un ensemble de sommets u tel que (z,u) doit etre supprime de aretes_LG_k_alpha
        
    """
    aretes_cliques = fct_aux.edges_in_cliques(cliques_couvertures);
    s_z = S_sommet(sommet_z = nom_sommet_1,
                   gamma_z = sommets_LG[nom_sommet_1].voisins,
                   aretes_LG_k_alpha = aretes_LG_k_alpha_cor,
                   cliques_couvertures = cliques_couvertures,
                  aretes_cliques = aretes_cliques);
    
    # determination de C1 = (C_1,C_2) avec C_1, C_2 contratables
    dico_C1_C2_S1 = dict(); cpt = 0;
    
#    cliques_sommet_1 = cliques_par_nom_sommets[nom_sommet_1];
    cliques_sommet_1 = set(cliques_par_nom_sommets[nom_sommet_1]);
    cliques_voisines_sommet_1 = clique_voisine_sommet_z(
                                    sommet_z = nom_sommet_1,
                                    C = cliques_couvertures,
                                    cliques_sommet_z = cliques_sommet_1)
    
    cliques_contractables_s = cliques_contractables(
                                nom_sommet_z = nom_sommet_1, 
                                aretes_LG_k_alpha = aretes_LG_k_alpha_cor, 
                                aretes_cliques = aretes_cliques, 
                                cliques_sommet_z = cliques_sommet_1, 
                                cliques_voisines_z = cliques_voisines_sommet_1,
                                C = cliques_couvertures,
                                DBG = DBG)

    random.shuffle(cliques_contractables_s)
    
    logger.debug(
        "****** compres =>" \
        +" nom_sommet_z={}".format(nom_sommet_1) \
        +" cliques_contractables_s={}".format(len(cliques_contractables_s)) \
        +" cliques_voisines_sommet_1 = {}".format(len(cliques_voisines_sommet_1))
        ) if DBG else None;
    
    for clique_C1_C2_Cx in cliques_contractables_s[:MAX_NUMBER]:
        # construction de dico_C1_C2_S1
        #        dico_C1_C2_S1[(cpt, (C1,C2,...), (C3,C4,...))] = {
        #          "cliques_contratables_1":(C1, C2),
        #          "cliques_contratables_2":(C3, C4),
        #          "clique_possible_1": ,
        #          "clique_possible_2": ,
        #                   }
        dico_C1_C2_S1[(cpt, clique_C1_C2_Cx, frozenset())] = {
                             "cliques_contractables_1" : clique_C1_C2_Cx,
                             "cliques_contractables_2" : frozenset(),
                             "clique_possible_1" : \
                                 frozenset.union(
                                            *clique_C1_C2_Cx).union(
                                                    frozenset({nom_sommet_1})
                                                                ),
                             "clique_possible_2" : frozenset()
                            }
        cpt += 1;
        
    ## *chercher les paires de cliques contractables tel que 
    ## *  |contr1 \cap contr2 |= 1
    logger.debug("****** compres => Avant " \
                 +" nom_sommet_z={}".format(nom_sommet_1) \
                 +" dico_C1_C2_S1={}".format(len(dico_C1_C2_S1))
                ) if DBG else None;      
          
    for clique_p1_p2 in it.combinations(cliques_contractables_s[:MAX_NUMBER], 2):
        clique_p1 = frozenset.union(*clique_p1_p2[0]);
        clique_p2 = frozenset.union(*clique_p1_p2[1]);
        if cpt > MAX_NUMBER:
            break;
        if len(clique_p1.intersection(clique_p2)) == 1 and \
            clique_p1.intersection(clique_p2) == frozenset({nom_sommet_1}):
            cpt += 1;
            dico_C1_C2_S1[(cpt, clique_p1, clique_p2)] = {
                            "cliques_contractables_1" : clique_p1_p2[0],
                            "cliques_contractables_2" : clique_p1_p2[1],
                            "clique_possible_1" : frozenset.union(
                                                    clique_p1).union(
                                                    frozenset({nom_sommet_1})
                                                                ),
                            "clique_possible_2" : frozenset.union(
                                                    clique_p2).union(
                                                    frozenset({nom_sommet_1})
                                                                )
                            }
    print("****** compres => Avant " \
                 +" nom_sommet_z={}".format(nom_sommet_1) \
                 +" dico_C1_C2_S1={}".format(len(dico_C1_C2_S1))
                ) if DBG else None;  
    logger.debug("****** compres => Avant " \
                 +" nom_sommet_z={}".format(nom_sommet_1) \
                 +" dico_C1_C2_S1={}".format(len(dico_C1_C2_S1))
                ) if DBG else None;   
            
    
    
    # determination de pi1_pi2_ps
    nb_prod_cartesien = len(dico_C1_C2_S1);
    nbre_elts_pi1_pi2 = math.ceil( nb_prod_cartesien * number_items_pi1_pi2);
    cpt_prod_cartesien = 0;
    dico_p1_p2_ps = dict();
    
    if not dico_C1_C2_S1:
        ens_cliq_a_supprimer, aretes_ps = set(), set();
        dico_sommets_corriges, dico_sommets_non_corriges = dict(), dict();
        cliques_par_nom_sommets_new = dict();
        
        cliqs_couv_new, \
        aretes_LG_k_alpha_new, \
        dico_sommets_corriges, \
        dico_sommets_non_corriges, \
        cliques_par_nom_sommets_new, \
        sommets_LG_new = \
            mise_a_jour_aretes_cliques(
                    nom_sommet_z = nom_sommet_1,
                    cliques_couvertures_new = set(cliques_couvertures).copy(), 
                    aretes_LG_k_alpha_new = set(aretes_LG_k_alpha_cor).copy(), 
                    aretes_ps = aretes_ps,
                    noms_sommets_1 = noms_sommets_1.copy(),
                    sommets_LG = sommets_LG,
                    cliques_par_nom_sommets = cliques_par_nom_sommets.copy())
        
        dico_p1_p2_ps[cpt_prod_cartesien] = {
                    "id_sommet_1": id_sommet_1,
                    "nom_sommet_1": nom_sommet_1,
                    "p1": frozenset(),
                    "p2": frozenset(),
                    "ps": frozenset(),
                    "S_z": s_z,
                    "aretes_ajoutees_p1": frozenset(),
                    "nbre_aretes_ajoutees_p1": np.inf,
                    "aretes_ajoutees_p2": frozenset(),
                    "nbre_aretes_ajoutees_p2": np.inf,
                    "aretes_supprimees_ps": frozenset(),
                    "nbre_aretes_supprimees_ps": np.inf,
                    "aretes_LG_k_alpha_new": aretes_LG_k_alpha_new,
                    "cliques_couvertures_new": cliqs_couv_new,
                    "sommets_LG_new": sommets_LG_new,
                    "sommets_corriges": dico_sommets_corriges,
                    "sommets_non_corriges": dico_sommets_non_corriges,
                    "cliques_par_nom_sommets_new": cliques_par_nom_sommets_new,
                    "cliques_supprimees": ens_cliq_a_supprimer,
                    "cliques_contractables_1": frozenset(),
                    "cliques_contractables_2": frozenset()
                            }
    else:
        for k_c1_c2_s1, val_cpt_c1_c2_s1 in dico_C1_C2_S1.items():
            cpt_prod_cartesien += 1;
            p1 = None; p2 = None;
            if cpt_prod_cartesien > nbre_elts_pi1_pi2:
                break;
                
            if val_cpt_c1_c2_s1["cliques_contractables_1"] and \
                not val_cpt_c1_c2_s1["cliques_contractables_2"]:
                p1 = val_cpt_c1_c2_s1["clique_possible_1"];
                p2 = frozenset();
            elif val_cpt_c1_c2_s1["cliques_contractables_1"] and \
                val_cpt_c1_c2_s1["cliques_contractables_2"]:
                p1 = val_cpt_c1_c2_s1["clique_possible_1"];
                p2 = val_cpt_c1_c2_s1["clique_possible_2"];
            else :
                print("IMPOSSIBLE cliques_contr_1 ={}, cliques_contr_2={}".format(
                      len(val_cpt_c1_c2_s1["cliques_contractables_1"]),
                      len(val_cpt_c1_c2_s1["cliques_contractables_2"])))
            
            if p1 is not None and p2 is not None :
                gamma_1 = sommets_LG[nom_sommet_1].voisins
                ps = gamma_1 \
                     - val_cpt_c1_c2_s1["clique_possible_1"].intersection(gamma_1) \
                     - val_cpt_c1_c2_s1["clique_possible_2"].intersection(gamma_1);
                
                aretes_ps = set( frozenset((nom_sommet_1, sommet_ps)) 
                                    for sommet_ps in ps
                                )
                aretes_p1 = set( map(frozenset, it.combinations(p1,2)) )
                aretes_ajoutees_p1 = aretes_differente(
                                        aretes_LG_k_alpha_cor, 
                                        aretes_p1);
            
                aretes_p2 = set( map(frozenset, it.combinations(p2,2)) )
                aretes_ajoutees_p2 = aretes_differente(
                                        aretes_LG_k_alpha_cor, 
                                        aretes_p2);
                                                
                aretes_LG_k_alpha_new = set(aretes_LG_k_alpha_cor).union(
                                                aretes_ajoutees_p1.union(
                                                    aretes_ajoutees_p2
                                                    )
                                                );
                                                
                cliques_couvertures_new = set(cliques_couvertures.copy());
                ens_cliq_a_supprimer = set();                                       
                for cliq_a_supps in [val_cpt_c1_c2_s1["cliques_contractables_1"],
                                     val_cpt_c1_c2_s1["cliques_contractables_2"]]:
                    for cliq_a_supp in cliq_a_supps:
                        ens_cliq_a_supprimer.add(cliq_a_supp);
                                           
                for cliq_couv_new in cliques_couvertures_new :
                    if cliq_couv_new.issubset(val_cpt_c1_c2_s1["clique_possible_1"]) or \
                        cliq_couv_new.issubset(val_cpt_c1_c2_s1["clique_possible_2"]) :
                        ens_cliq_a_supprimer.add(cliq_couv_new);
               
                cliques_couvertures_new.difference_update(ens_cliq_a_supprimer);
            
                cliques_couvertures_new.add( 
                                        val_cpt_c1_c2_s1["clique_possible_1"] );
                cliques_couvertures_new.add( 
                                        val_cpt_c1_c2_s1["clique_possible_2"] ) \
                          if val_cpt_c1_c2_s1["clique_possible_2"] else None;
        
        
                dico_sommets_corriges, dico_sommets_non_corriges = dict(), dict();
                cliques_par_nom_sommets_new = dict();
                
                cliqs_couv_new, \
                aretes_LG_k_alpha_new, \
                dico_sommets_corriges, \
                dico_sommets_non_corriges, \
                cliques_par_nom_sommets_new, \
                sommets_LG_new = \
                    mise_a_jour_aretes_cliques(
                        nom_sommet_z = nom_sommet_1,
                        cliques_couvertures_new = set(cliques_couvertures_new).copy(), 
                        aretes_LG_k_alpha_new = set(aretes_LG_k_alpha_new).copy(), 
                        aretes_ps = aretes_ps,
                        noms_sommets_1 = noms_sommets_1.copy(),
                        sommets_LG = sommets_LG,
                        cliques_par_nom_sommets = cliques_par_nom_sommets.copy())
                
                dico_p1_p2_ps[cpt_prod_cartesien] = {
                        "id_sommet_1": id_sommet_1,
                        "nom_sommet_1": nom_sommet_1,
                        "p1": val_cpt_c1_c2_s1["clique_possible_1"],
                        "p2": val_cpt_c1_c2_s1["clique_possible_2"],
                        "ps": ps,
                        "S_z": s_z,
                        "aretes_ajoutees_p1": aretes_ajoutees_p1,
                        "nbre_aretes_ajoutees_p1": len(aretes_ajoutees_p1),
                        "aretes_ajoutees_p2": aretes_ajoutees_p2,
                        "nbre_aretes_ajoutees_p2": len(aretes_ajoutees_p2),
                        "aretes_supprimees_ps": aretes_ps,
                        "nbre_aretes_supprimees_ps": len(aretes_ps),
                        "aretes_LG_k_alpha_new": aretes_LG_k_alpha_new.copy(),
                        "cliques_couvertures_new": cliqs_couv_new,
                        "sommets_LG_new": sommets_LG_new,
                        "sommets_corriges": dico_sommets_corriges,
                        "sommets_non_corriges": dico_sommets_non_corriges,
                        "cliques_par_nom_sommets_new": cliques_par_nom_sommets_new,
                        "cliques_supprimees" : ens_cliq_a_supprimer,
                        "cliques_contractables_1": set(val_cpt_c1_c2_s1["cliques_contractables_1"]),
                        "cliques_contractables_2": set(val_cpt_c1_c2_s1["cliques_contractables_2"])
                        } 
                
            pass # end for
        pass # end else
        
    logger.debug("****** compres ===> Fin compression " \
                 +" sommet_z : {}, ".format(nom_sommet_1) \
                 +" nbre_elts_pi1_pi2:{}, ".format(nbre_elts_pi1_pi2) \
                 +" dico_C1_C2_S1:{}, ".format(len(dico_C1_C2_S1)) \
                 +" dico_p1_p2_ps:{}".format(len(dico_p1_p2_ps))
          )  
    return dico_p1_p2_ps;
    pass
###############################################################################
#               calcul de la compression d un sommet => fin
###############################################################################

###############################################################################
#      critere selection compression: critere local et global => debut
###############################################################################
def critere_C2_C1_local(dico_compression, 
                        mode_correction,
                        critere_correction, 
                        DBG):
    """ 
    selectionner le dico selon C2 puis C1 parmi les compressions possibles 
        sommet a corriger (sommet a -1)
    
    C2 : le maximum de sommets corriges
        * choisir le sommet a -1 qui corrige le max de sommets a -1 possibles
    C1 : le minimum d'aretes corriges
    
    dico_compression : dictionnaire contenant les compressions (p1,p2,ps) du 
                        sommet sommet_z
    """
    max_c2 = 0;
    min_c1 = np.inf;
    dico_c1_c2 = dict();
    critere = "";
    
    if not dico_compression :
        #print("@@CritereLocal: dico_compression={}".format( len(dico_compression) ))
        return min_c1, max_c2, [];
    
    # definition de C2
    if critere_correction == "voisins_corriges":                               # C2
        critere = "C2";
        for cpt_prod_cartesien, dico_p1_p2_ps in dico_compression.items():
            if len(dico_p1_p2_ps["sommets_corriges"]) >= max_c2:
                max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
                nbre_aretes_corriges = \
                            dico_p1_p2_ps["nbre_aretes_ajoutees_p1"] + \
                            dico_p1_p2_ps["nbre_aretes_ajoutees_p2"] + \
                            dico_p1_p2_ps["nbre_aretes_supprimees_ps"];
                min_c1 = nbre_aretes_corriges if min_c1 >= nbre_aretes_corriges \
                                                else min_c1;
                if (min_c1,max_c2) not in dico_c1_c2:
                    dico_c1_c2[(min_c1,max_c2)] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2[(min_c1,max_c2)].append(dico_p1_p2_ps);
    
    # definition de C1
    elif critere_correction == "nombre_aretes_corrigees":                      # C1
        critere = "C1";
        for cpt_prod_cartesien, dico_p1_p2_ps in dico_compression.items():
            nbre_aretes_corriges = \
                        dico_p1_p2_ps["nbre_aretes_ajoutees_p1"] + \
                        dico_p1_p2_ps["nbre_aretes_ajoutees_p2"] + \
                        dico_p1_p2_ps["nbre_aretes_supprimees_ps"];
            min_c1 = nbre_aretes_corriges if min_c1 >= nbre_aretes_corriges \
                                          else min_c1;
            max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
            if (min_c1,max_c2) not in dico_c1_c2:
                dico_c1_c2[(min_c1,max_c2)] = [dico_p1_p2_ps];
            else:
                dico_c1_c2[(min_c1,max_c2)].append(dico_p1_p2_ps);
    
    # definition de C2 puis de C1
    elif critere_correction == "voisins_nombre_aretes_corrigees":              # C2_C1
        critere = "C2_C1"
        for cpt_prod_cartesien, dico_p1_p2_ps in dico_compression.items() :
            if len(dico_p1_p2_ps["sommets_corriges"]) >= max_c2 :
                max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
                nbre_aretes_corriges = \
                            dico_p1_p2_ps["nbre_aretes_ajoutees_p1"] + \
                            dico_p1_p2_ps["nbre_aretes_ajoutees_p2"] + \
                            dico_p1_p2_ps["nbre_aretes_supprimees_ps"];
                min_c1 = nbre_aretes_corriges if min_c1 >= nbre_aretes_corriges \
                                                else min_c1;
                if (min_c1,max_c2) not in dico_c1_c2:
                    dico_c1_c2[(min_c1,max_c2)] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2[(min_c1,max_c2)].append(dico_p1_p2_ps);
                    
    logger.debug("@@CritereLocal: critere={},".format(critere) \
                 +" min_c1={},".format(min_c1) \
                 +" max_c2={},".format(max_c2) \
                 +" cles_dico_c1_c2={},".format(set(dico_c1_c2.keys())) \
                 +" dico_c1_c2={},".format(len(dico_c1_c2[(min_c1,max_c2)])) \
                 +" dico_compression={}".format(len(dico_compression))) \
                 if DBG else None
                          
    if not dico_c1_c2:
        return min_c1, max_c2, [];
    else:
        return min_c1, max_c2, dico_c1_c2[(min_c1,max_c2)];
    pass

def rechercher_min_max(liste_tuples, critere):
    """ retourne la tuple (min, max)
    """
    min_c1 = np.inf;
    if len(liste_tuples) == 0:
        return np.inf
    
    #max_c2 = 0;
    if critere == "C1":
        return min(liste_tuples)
    elif critere == "C2":
        return max(liste_tuples)
    elif critere == "C2_C1":
        liste_intermediaires = [];
        min_c1 = min(liste_tuples)[0]
        for tuple_ in liste_tuples:
           if tuple_[0] == min_c1:
               liste_intermediaires.append(tuple_)
        return max(liste_intermediaires)

def critere_C2_C1_global(dico_compression,
                         mode_correction,
                         critere_correction,
                         DBG):
    """ recherche la compression optimale parmi tous les sommets a corriger 
        selon les criteres C1 et C2.
        
    C2 : le maximum de sommets corriges
        * choisir le sommet a -1 qui corrige le max de sommets a -1 possibles
    C1 : le minimum d'aretes corriges 
    
    methode de selection C2
        je cherche le min local de c1 pour tous les sommets a corriger
        parmi les min locaux, je cherche le max global de c2
        une fois la liste des (min_global,max_global), je prends le 1er element.
    """
    
    max_c2_global = 0;
    min_c1_global = np.inf;
    dico_c1_c2_global = dict();
    cle_min_max_c2 = None;
    
    critere = ""
    
    if critere_correction == "voisins_corriges":                               # C2
        critere = "C2"
        for id_sommet_z, dicos_p1_p2_ps in dico_compression.items():
            # selection de dico selon C1
            min_c1_local = dicos_p1_p2_ps[0];
            max_c2_local = dicos_p1_p2_ps[1];
            
            for dico_p1_p2_ps in dicos_p1_p2_ps[2]:
                nbre_aretes_corriges = \
                                dico_p1_p2_ps["nbre_aretes_ajoutees_p1"] \
                                + dico_p1_p2_ps["nbre_aretes_ajoutees_p2"] \
                                + dico_p1_p2_ps["nbre_aretes_supprimees_ps"];
                min_c1_local = nbre_aretes_corriges \
                                if min_c1_local >= nbre_aretes_corriges \
                                else min_c1_local;
                if (min_c1_local,max_c2_local) not in dico_c1_c2_global :
                    dico_c1_c2_global[(min_c1_local, 
                                       max_c2_local)] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2_global[(min_c1_local, 
                                       max_c2_local)].append(dico_p1_p2_ps);
                           
        # selection selon C2
        cle_min_max_c2 = rechercher_min_max(dico_c1_c2_global.keys(), "C2");
    
    elif critere_correction == "nombre_aretes_corrigees":                      # C1
        critere = "C1"
        for id_sommet_z, dicos_p1_p2_ps in dico_compression.items() :
            # selection de dico selon C2
            max_c2_local = dicos_p1_p2_ps[1];
            min_c1_local = dicos_p1_p2_ps[0];
            
            for dico_p1_p2_ps in dicos_p1_p2_ps[2] :
                nbre_sommets_corriges = len(dico_p1_p2_ps["sommets_corriges"]);
                max_c2_local = nbre_sommets_corriges \
                                if nbre_sommets_corriges > max_c2_local \
                                else max_c2_local;
                if (min_c1_local,max_c2_local) not in dico_c1_c2_global:
                    dico_c1_c2_global[(min_c1_local, 
                                       max_c2_local)] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2_global[(min_c1_local, 
                                       max_c2_local)].append(dico_p1_p2_ps);

        # selection selon C1
        cle_min_max_c2 = rechercher_min_max(dico_c1_c2_global.keys(), "C1");
        
    elif critere_correction == "voisins_nombre_aretes_corrigees":              # C2_C1
        critere = "C2_C1"
        for id_sommet_z, dicos_p1_p2_ps in dico_compression.items():
            min_c1_local = dicos_p1_p2_ps[0]; #np.inf
            max_c2_local = dicos_p1_p2_ps[1]; #0
            
            for dico_p1_p2_ps in dicos_p1_p2_ps[2] :
                nbre_sommets_corriges = len(dico_p1_p2_ps["sommets_corriges"]);
                max_c2_local = nbre_sommets_corriges \
                                if nbre_sommets_corriges > max_c2_local \
                                else max_c2_local;
                nbre_aretes_corriges = \
                                dico_p1_p2_ps["nbre_aretes_ajoutees_p1"] \
                                + dico_p1_p2_ps["nbre_aretes_ajoutees_p2"] \
                                + dico_p1_p2_ps["nbre_aretes_supprimees_ps"];
                min_c1_local = nbre_aretes_corriges \
                                if min_c1_local >= nbre_aretes_corriges \
                                else min_c1_local;
                if (min_c1_local,max_c2_local) not in dico_c1_c2_global:
                    dico_c1_c2_global[(min_c1_local, 
                                       max_c2_local)] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2_global[(min_c1_local, 
                                       max_c2_local)].append(dico_p1_p2_ps);
                
        # selection selon C2_C1
        cle_min_max_c2 = rechercher_min_max(dico_c1_c2_global.keys(), "C2_C1");
    
    numero_sol_c1_c2 = np.random.randint(
                        low=0, 
                        high=len(dico_c1_c2_global[cle_min_max_c2])
                        )
    
    min_c1_global = cle_min_max_c2[0];
    max_c2_global = cle_min_max_c2[1];
    
    logger.debug("@@CritereGlobal critere={}, ".format(critere) \
                 +" dico_compression={}, ".format(len(dico_compression)) \
                 +" cle_globale={}, ".format(set(dico_c1_c2_global.keys())) \
                 +" cle_min_max_c2={}, ".format(cle_min_max_c2) \
                 +" nbre_sol_c1_c2={}".format(
                         len(dico_c1_c2_global[cle_min_max_c2]))
                 ) if DBG else None;
                 
    return min_c1_global, \
            max_c2_global, \
            dico_c1_c2_global[cle_min_max_c2][numero_sol_c1_c2];
###############################################################################
#      critere selection compression: critere local et global => fin
###############################################################################

###############################################################################
#                application de la correction => debut
###############################################################################
def appliquer_correction(dico_sol_C2_C1,
                         sommets_1,
                         DBG):
    """ appliquer la compression choisie dans le graphe.
    """
    """
    {
    "id_sommet_1": id_sommet_1,
    "sommet_1": nom_sommet_1,
    "p1": val_cpt_c1_c2_s1["clique_possible_1"],
    "p2": val_cpt_c1_c2_s1["clique_possible_2"],
    "ps": ps,
    "S_z": s_z,
    "aretes_ajoutees_p1": aretes_ajoutees_p1,
    "nbre_aretes_ajoutees_p1": len(aretes_ajoutees_p1),
    "aretes_ajoutees_p2": aretes_ajoutees_p2,
    "nbre_aretes_ajoutees_p2": len(aretes_ajoutees_p2),
    "aretes_supprimees_ps": aretes_ps,
    "nbre_aretes_supprimees_ps": len(aretes_ps),
    "aretes_LG_k_alpha_new": aretes_LG_k_alpha_new.copy(),
    "cliques_couvertures_new": cliqs_couv_new,
    "sommets_LG_new": sommets_LG_new,
    "sommets_corriges": dico_sommets_corriges,
    "sommets_non_corriges": dico_sommets_non_corriges,
    "cliques_par_nom_sommets_new": cliques_par_nom_sommets_new,
    "cliques_supprimees" : ens_cliq_a_supprimer,
    "cliques_contractables_1" : set(val_cpt_c1_c2_s1["cliques_contractables_1"]),
    "cliques_contractables_2" : set(val_cpt_c1_c2_s1["cliques_contractables_2"])
    } 
    """
    
    cliques_couvertures = set();
    cliques_couvertures = dico_sol_C2_C1['cliques_couvertures_new'];
    aretes_LG_k_alpha = dico_sol_C2_C1['aretes_LG_k_alpha_new'];
    sommets_LG = dico_sol_C2_C1['sommets_LG_new'];
    cliques_par_nom_sommets = dico_sol_C2_C1['cliques_par_nom_sommets_new'];
    
    id_sommets_1 = set(dico_sol_C2_C1["sommets_corriges"].keys());
    id_sommets_1.add(dico_sol_C2_C1["id_sommet_1"]);
    sommets_corriges = dico_sol_C2_C1["sommets_corriges"].values();
    
    logger.debug(
            "*** Avant correction : id_sommets_1:{}, ".format(id_sommets_1) \
            +" sommets_corriges={}, ".format(sommets_corriges) \
            +" sommet_1={}".format(dico_sol_C2_C1["nom_sommet_1"])) \
        if DBG else None;
                
    sommets_1 = np.delete(sommets_1, list(id_sommets_1)).tolist();
    logger.debug("*** Apres correction : "
                  +"sommets_1 restants = {}".format(sommets_1)) \
        if DBG else None;
                     
    if set(sommets_1).intersection(set(sommets_corriges)) :
        print("---ERROR : sommets {} suppression : NOK -----".
              format(sommets_corriges))
                     
    """
    cliques_couv_new,\
            aretes_LG_k_alpha_cor_new,\
            sommets_LG,\
            noms_sommets_1
    """
    return cliques_couvertures,\
            aretes_LG_k_alpha,\
            sommets_LG,\
            cliques_par_nom_sommets,\
            sommets_1;
###############################################################################
#               application de la correction => fin
###############################################################################

###############################################################################
#                correction aleatoire des sommets_1 => debut
###############################################################################
def correction_aleatoire(sommets_1, 
                         aretes_LG_k_alpha_cor,
                         dico_correction,
                         dico_parametres_new):
    """
    correction des sommets a corriger par une selection aleatoire de ceux-ci.
    """
    cliques_couvertures_cor = dico_correction["C"].copy();
    cliques_par_sommets = dico_correction["sommets_par_cliqs_avec_aretes"].copy();
    
    mode_select_noeuds_1 = dico_parametres_new["mode_select_noeuds_1s"][0];
    critere_select_pi1_pi2 = dico_parametres_new["critere_selection_pi1_pi2"];
    number_items_pi1_pi2 = dico_parametres_new["number_items_pi1_pi2"];
    DBG = dico_parametres_new["DBG"]
    
    dico_sommets_corriges = dict();
    
    id_sommet_1 = 0;
    while len(sommets_1)>0:
        sommet_1 = sommets_1.pop(id_sommet_1);
        
        dicos_p1_p2_ps = dict();
        dicos_p1_p2_ps = compression_sommet(id_sommet_1,
                            sommet_1 = sommet_1,
                            sommets_1 = sommets_1,
                            cliques_par_sommets = cliques_par_sommets,
                            cliques_couv_cor = cliques_couvertures_cor,
                            aretes_LG_k_alpha_cor = aretes_LG_k_alpha_cor,
                            mode_select_noeuds_1 = mode_select_noeuds_1,
                            critere_select_pi1_pi2 = critere_select_pi1_pi2,
                            number_items_pi1_pi2 = number_items_pi1_pi2,
                            DBG = DBG)
                                           
        
        dico_sol_p1_p2_ps, dico_compression = dict(), dict(); 
        min_c1 = 0; max_c2 = 0;
        min_c1, max_c2, dico_compression = critere_C2_C1_local(
                                            dicos_p1_p2_ps,
                                            mode_select_noeuds_1,
                                            critere_select_pi1_pi2,
                                            DBG = DBG)
        
        print("\n\n len(dico_compression)={}".format(len(dico_compression)))
        nbre_alea = random.randint(
                        0,
                        len(dico_compression)-1)
        dico_sol_p1_p2_ps = dico_compression[nbre_alea];
        
        if not dico_sol_p1_p2_ps:
            cout_T = {"aretes_ajoutees_p1": frozenset(),
                      "aretes_ajoutees_p2": frozenset(),
                      "aretes_supprimees": frozenset(),
                      "min_c1": min_c1,
                      "max_c2": max_c2};
            id_sommet_1 += 1;
            dico_sommets_corriges[("0_0", "0_0")] = {
                        "compression_p1": frozenset(),
                        "compression_p2": frozenset(),
                        "compression_ps": frozenset(),
                        "sommets_corriges": dict(), # voisins_corriges = {"id_voisin_ds_sommets_a_corriger":voisin}
                        "cout_T": cout_T
                        }
            sommets_1 = list();
        else:
            cliques_couv_new, \
            aretes_LG_k_alpha_cor_new, \
            sommets_LG_new, \
            cliques_par_nom_sommets_new, \
            sommets_1 = appliquer_correction(
                                dico_sol_p1_p2_ps,
                                sommets_1, 
                                DBG);
            
            # mise a jour variables
            cliques_par_sommets = cliques_par_nom_sommets_new.copy()
            cliques_couvertures_cor = cliques_couv_new.copy()
            aretes_LG_k_alpha_cor = aretes_LG_k_alpha_cor_new.copy();
            sommets_LG = sommets_LG_new.copy();
            cout_T = {
                "aretes_ajoutees_p1":dico_sol_p1_p2_ps["aretes_ajoutees_p1"],
                "aretes_ajoutees_p2":dico_sol_p1_p2_ps["aretes_ajoutees_p2"],
                "aretes_supprimees":dico_sol_p1_p2_ps["aretes_supprimees_ps"],
                "min_c1":min_c1, "max_c2":max_c2};
            id_sommet_1 += 1;
            dico_sommets_corriges[(id_sommet_1, 
                                   dico_sol_p1_p2_ps["sommet_1"])] = {
                        "compression_p1":dico_sol_p1_p2_ps["p1"],
                        "compression_p2":dico_sol_p1_p2_ps["p2"],
                        "compression_ps":dico_sol_p1_p2_ps["ps"],
                        "sommets_corriges":dico_sol_p1_p2_ps["sommets_corriges"], # voisins_corriges = {"id_voisin_ds_sommets_a_corriger":voisin}
                        "cout_T": cout_T
                        }
        
    # TODO A REVOIR JE DOIS RETOURNER dico_correction
    dico_correction["C"] = cliques_couvertures_cor;
    dico_correction["sommets_par_cliqs_avec_aretes"] = cliques_par_sommets;
    dico_correction["0"] = 0;
    return dico_correction, \
            dico_sommets_corriges;
#            sommets_LG, \
            
    pass
###############################################################################
#               correction aleatoire des sommets_1 => fin
###############################################################################

###############################################################################
#               algorithme de correction pour k entier naturel => debut
###############################################################################
def correction_cliques_k(dico_correction, 
                         aretes_matE_k_alpha,
                         dico_gamma_noeud, 
                         dico_parametres_new):
    """
    corriger les cliques selon le mode de selection des sommets a -1
    """
    aretes_LG_k_alpha_cor = aretes_matE_k_alpha.copy();
    
    sommets_a_corriger = dico_correction["sommets_a_corriger"];
    dico_sommets_corriges = dict();
    
    if dico_correction["critere_selection_pi1_pi2"] == "voisins_corriges" \
        and dico_correction["mode_select_noeuds_1s"] == "aleatoire":
        dico_correction, dico_sommets_corriges = correction_aleatoire(
                                                    sommets_a_corriger, 
                                                    aretes_LG_k_alpha_cor,
                                                    dico_correction,
                                                    dico_parametres_new);
        
    elif dico_correction["critere_selection_pi1_pi2"] == "nombre_aretes_corrigees" \
        and dico_correction["mode_select_noeuds_1s"] == "aleatoire":
        dico_correction, dico_sommets_corriges = correction_aleatoire(
                                                    sommets_a_corriger, 
                                                    aretes_LG_k_alpha_cor,
                                                    dico_correction,
                                                    dico_parametres_new);
        
    elif dico_correction["critere_selection_pi1_pi2"] == "voisins_nombre_aretes_corrigees" \
        and dico_correction["mode_select_noeuds_1s"] == "aleatoire":
        dico_correction, dico_sommets_corriges = correction_aleatoire(
                                                    sommets_a_corriger, 
                                                    aretes_LG_k_alpha_cor,
                                                    dico_correction,
                                                    dico_parametres_new);
                
    return dico_correction, dico_sommets_corriges;
###############################################################################
#               algorithme de correction pour k entier naturel => fin
###############################################################################