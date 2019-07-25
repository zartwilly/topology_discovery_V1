#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 16:06:00 2019

@author: willy
"""
import math;
import logging;

import numpy as np;
import itertools as it;

import fonctions_auxiliaires as fct_aux;
###############################################################################
#                fonctions de bases pour la correction => debut
###############################################################################
def mise_a_jour_aretes_cliques(C_new, 
                               aretes_Ec, 
                               aretes_ps, 
                               sommets_a_corriger, 
                               dico_sommets_par_cliqs):
    """ mettre a jour les sommets par cliques puis 
        verifier les sommets couverts par plus de deux cliques.
        
    """
#    print("PS: aretes_ps={}, cpt_aretes_Ec={}".format(aretes_ps, len(aretes_Ec)))
    
    # suppression des aretes_ps dans aretes_Ec
    aretes_ps = set(map(tuple, aretes_ps));
    aretes_ps_new = set();
    for arete in aretes_ps:
        aretes_ps_new.add(arete)
        aretes_ps_new.add( (arete[1],arete[0]) )
    aretes_Ec.difference_update(aretes_ps_new);
#    print("PS new cpt_aretes_Ec={}".format(len(aretes_Ec)))
    
    # suppression cliques dont on a supprime des aretes_ps
    C_nouvelle = C_new.copy();
    for c in C_new:
        for arete_ps in aretes_ps:
            if set(arete_ps).issubset(c):
                C_nouvelle.difference_update({c});
                temp_aretes = [arete_ps, (arete_ps[1],arete_ps[0])]
                aretes_Ec.difference_update(set(temp_aretes));             
                break;
                
#    print("PS cpt_C_new={},cpt_C_nouvelle={}".format(len(C_new), len(C_nouvelle)))
    dico_sommets_par_cliqs_new = fct_aux.couverture_par_sommets(
                                    sommets_matE = dico_sommets_par_cliqs.keys(),
                                    C = C_nouvelle);
    dico_sommets_corriges = dict(); dico_sommets_non_corriges = dict();
    
    for id_sommet, sommet_a_corriger in enumerate(sommets_a_corriger):
        cliques_sommet_a_corr = dico_sommets_par_cliqs_new[sommet_a_corriger];
        gamma_sommet_a_corriger = set(fct_aux.gamma(aretes_Ec,
                                                    sommet_a_corriger));
        if len(cliques_sommet_a_corr) == 0:
            dico_sommets_non_corriges[id_sommet] = sommet_a_corriger;
            
        elif len(cliques_sommet_a_corr) == 1 and \
            cliques_sommet_a_corr[0] == gamma_sommet_a_corriger:
            dico_sommets_corriges[id_sommet] = sommet_a_corriger;
            
        elif len(cliques_sommet_a_corr) == 1 and \
            cliques_sommet_a_corr[0] != gamma_sommet_a_corriger:
            dico_sommets_non_corriges[id_sommet] = sommet_a_corriger;
            
        elif len(cliques_sommet_a_corr) == 2:
            dico_sommets_corriges[id_sommet] = sommet_a_corriger;
            
        elif len(cliques_sommet_a_corr) > 2:
            dico_sommets_non_corriges[id_sommet] = sommet_a_corriger;
#    print("??? dico_sommets_corriges={},dico_sommets_non_corriges={}, dico_sommets_par_cliqs_new={}".\
#          format(dico_sommets_corriges, dico_sommets_non_corriges, dico_sommets_par_cliqs_new))
            
    return C_nouvelle, \
            aretes_Ec,\
            dico_sommets_corriges, \
            dico_sommets_non_corriges, \
            dico_sommets_par_cliqs_new;
    

###############################################################################
#                fonctions de bases pour la correction => fin
###############################################################################

###############################################################################
#               cliques contractables et S_sommet => debut
###############################################################################
def S_sommet(sommet_z, gamma_z, aretes_Ec, C, aretes_cliques):
    """ 
    voisins v de sommet_z tels que 
        * {v, sommet_z} est une clique de C
        * {v, sommet_z} de aretes_Ec n'est couverte par aucune clique de C.
        
    """
    logger = logging.getLogger('S_sommet');
    S_z = list();
    for voisin_z in gamma_z :
        if {voisin_z, sommet_z} in C :
            S_z.append(voisin_z);
        elif ((voisin_z, sommet_z) in aretes_Ec or \
            (sommet_z, voisin_z) in aretes_Ec) and \
            ((voisin_z, sommet_z) not in aretes_cliques and \
            (sommet_z, voisin_z) not in aretes_cliques):
            S_z.append(voisin_z);
    logger.debug(" * S_z: {}".format(S_z))
    return S_z;
    
def is_contractable(clique1, clique2, aretes_Ec, aretes_cliques, C):
    """ 
    determine si deux cliques sont contractables. 
    
    if true : cliques 1 et 2 sont contractables
    if false : sinon.
    """
    boolean_contractable = True; 
    cliques_suppl_contractables = set()
    sommet_inter = clique1.intersection(clique2)
    for noeud1, noeud2 in it.product(clique1-sommet_inter, 
                                     clique2-sommet_inter):    
        if noeud1 != noeud2 and \
           ((noeud1,noeud2) in aretes_Ec or (noeud2,noeud2) in aretes_Ec) and \
           ( frozenset((noeud1,noeud2)) not in C) :
               boolean_contractable = False; 
               cliques_suppl_contractables = set()
               break;
        if frozenset((noeud1,noeud2)) in C:
           cliques_suppl_contractables.add( frozenset((noeud1,noeud2)) ) 
    return boolean_contractable, cliques_suppl_contractables;
    
def cliques_contractables(sommet_z, 
                          aretes_Ec, 
                          aretes_cliques, 
                          cliques_sommet_z, 
                          C):
    """ retourne la liste des cliques contractables autour du sommet_z. """
    
    logger = logging.getLogger('cliques_contractables');
    cliq_contractables = [];
    """
    TODO a poser la question a un forum.
    how to do itertools with Nonetype inside a list.
    for c1, c2 in it.combinations(cliques_sommet_z.append({}),2):
        if not c1 or not c2 :
            cliq_contractables.append((c1, c2));
        else:
            if is_contractable(c1, c2, aretes_Ec, aretes_cliques):
                cliq_contractables.append((c1, c2))
    return cliq_contractables;
    """
    for c1, c2 in it.combinations(cliques_sommet_z, 2):
        boolean_contractable, cliques_suppl_contractables = is_contractable(
                                                            c1, 
                                                            c2, 
                                                            aretes_Ec, 
                                                            aretes_cliques, 
                                                            C)
        if boolean_contractable:
            cliq_contractables.append((c1, c2, cliques_suppl_contractables))
    for c in cliques_sommet_z:
        cliq_contractables.append((c,frozenset(), set()))
        
    logger.debug(" * cliq_contractables: {}".format( len(cliq_contractables)) )
    return cliq_contractables;
###############################################################################
#               cliques contractables et S_sommet => fin
###############################################################################

###############################################################################
#           voisinage sommet, dependance et augmentation  => debut
###############################################################################
def voisine_sommet(sommet_z, cliques_sommet_z, cliques_not_sommet_z, s_z):
    """ cliques voisines du sommet_z. """
#    cliques_sommet_z = [frozenset(c) for c in cliques_sommet_z]
    logger = logging.getLogger('voisine_sommet');
    cliques_voisines = [c for c in cliques_not_sommet_z 
                        if len(c.intersection(set(s_z))) >=1]
    logger.debug(" * voisine_sommet: {}".format(len(cliques_voisines)))
    return cliques_voisines;

def dependance_sommet(sommet_z, gamma_z, cliques_sommet_z, clique_voisine):
    """ retourner les cliques dependantes d une clique clique_voisine. """
    return [cliq for cliq in cliques_sommet_z \
            if len(cliq.intersection(clique_voisine.intersection(gamma_z))) != 0]
    
def augmentation(sommet_z, 
                 gamma_z, 
                 cliques_sommet_z, 
                 s_z, 
                 args):
    """ retourne les augmentations possibles autour du sommet sommet_z. """
    
    logger = logging.getLogger('augmentation');
    # cliques_sommet_z = cliques_sommet(sommet_z, args["dico_sommets_par_cliqs"]);
    cpt = 0;
    dico_cliques_augmentante = dict();
    cliques_not_sommet_z = list();
    cliques_not_sommet_z = [ c for c in args["C"] 
                             if not c.intersection({sommet_z})]
    
    cliques_voisine = voisine_sommet(sommet_z, 
                                     cliques_sommet_z, 
                                     cliques_not_sommet_z, 
                                     s_z);
    for clique_voisine in cliques_voisine:
        cliques_dependante = list();
        cliques_dependante = dependance_sommet(sommet_z, 
                                               gamma_z, 
                                               cliques_sommet_z, 
                                               clique_voisine);
#        print("!!!clique_voisine={},\n cliques_dependante={}".\
#              format(clique_voisine,cliques_dependante))
        if not cliques_dependante:
#                dico_cliques_augmentante[cpt] = {"cliq":cliq, 
#                                                "voisine":clique_voisine,
#                                                "dependante":frozenset(),
#                                                "sommet_z":sommet_z}
            dico_cliques_augmentante[(cpt, clique_voisine,\
                                     frozenset(), frozenset())] = {
                                      "voisine":clique_voisine,
                                      "dependante":frozenset(),
                                      "cliques_suppl_contractables": frozenset(),
                                      "sommet_z":sommet_z
                                      };                                  
            cpt += 1;
        else:
            for clique_dependante in cliques_dependante:
                boolean_contractable, cliques_suppl_contractables = \
                        is_contractable(clique_voisine, 
                                        clique_dependante, 
                                        args["aretes_Ec"], 
                                        args["aretes_cliques"],
                                        args["C"]);
                if boolean_contractable:
                    dico_cliques_augmentante[(cpt, clique_voisine,\
                        clique_dependante, \
                        frozenset(cliques_suppl_contractables))] = {
                        "voisine": clique_voisine,
                        "dependante": clique_dependante,
                        "cliques_suppl_contractables": \
                            frozenset(cliques_suppl_contractables),
                        "sommet_z": sommet_z
                        }                        
                    cpt += 1;
                    
#    print("??? dico_cliques_augmentante={}".format(dico_cliques_augmentante));
    logger.debug(" * augmentation: {}".format(len(dico_cliques_augmentante)))
    return dico_cliques_augmentante;
###############################################################################
#          voisinage sommet, dependance et augmentation => fin
###############################################################################
    
###############################################################################
#               compression d'un sommet => debut
###############################################################################
def compression_sommet(id_sommet_z, 
                       sommet_z, 
                       sommets_a_corriger, 
                       cliques_sommet_z, 
                       args):
    """ retourne la compression d'un sommet sommet_z. 
    
    la compression est le triplet (pi1, pi2, ps) dans lequel 
        * pi1, pi2 sont des cliques qui fusionnent 
            - des cliques augmentantes C1, C2 ou 
            - des cliques contractables C1, C2 ou 
            - un ensemble S1 tel que S1 n'est contenu par aucune clique C1 ou C2
        * pi1, pi2 sont des augmentations
        * ps est un ensemble de sommets u tel que (z,u) doit etre supprime de aretes_Ec
        
    """
    logger = logging.getLogger('compression_sommet');
#    print("X cliques_sommet_z={}".format(cliques_sommet_z))
    s_z = S_sommet(sommet_z, 
                   args["dico_gamma_sommets"][sommet_z][1], 
                   args["aretes_Ec"], 
                   args["C"], 
                   args["aretes_cliques"]);
    
    # determination de C1 = (C_1,C_2) avec C_1, C_2 contratables
    dico_C1_C2_S1 = dict(); cpt = 0;
    for C1, C2, cliques_suppl_contractables in cliques_contractables(
                                       sommet_z, 
                                       args["aretes_Ec"], 
                                       args["aretes_cliques"], 
                                       cliques_sommet_z, 
                                       args["C"]):
        S1 = C1.union(C2) - C1.union(C2).intersection(s_z);
        bool_sommet_a_exclu = True; S1_new = frozenset();
        for sommet_S1 in S1:
            for s1_, c1_c2 in it.product(frozenset({sommet_S1}), C1.union(C2)):
                if frozenset({s1_, c1_c2}) not in args["C"]:
                    bool_sommet_a_exclu = False;
                    break;
            if bool_sommet_a_exclu:
                S1_new.union({s1_})
                
        dico_C1_C2_S1[(cpt, C1, C2, S1_new)] =\
                {
                "cliques_contratables":(C1,C2),
                "cliques_suppl_contractables":cliques_suppl_contractables,
                "S1":S1,
                "clique_possible": 
                        C1.union(C2.union(S1_new.union(frozenset({sommet_z}))))
                                            }
        cpt += 1;
    
    
    # determination de pi1_pi2_ps
    dico_cliques_augmentante = dict();
    dico_cliques_augmentante = augmentation(
                                    sommet_z,
                                    args["dico_gamma_sommets"][sommet_z][1], 
                                    cliques_sommet_z, 
                                    s_z, 
                                    args);
            
    nb_prod_cartesien = pow(len(dico_C1_C2_S1), len(dico_cliques_augmentante)) \
                        if len(dico_C1_C2_S1) >= len(dico_cliques_augmentante) \
                        else pow(len(dico_cliques_augmentante),len(dico_C1_C2_S1))
    nbre_elts_pi1_pi2 = math.ceil( nb_prod_cartesien * 
                                  args["number_items_pi1_pi2"]);
    nbre_elts_pi1_pi2 = 200 if nbre_elts_pi1_pi2 > 1024 else nbre_elts_pi1_pi2;
                           
    cpt_prod_cartesien = 0;
    dico_p1_p2_ps = dict();
    print(" sommet_z ={}, ".format(sommet_z) \
         +" dico_C1_C2_S1:{}, ".format(len(dico_C1_C2_S1)) \
         +" dico_cliques_augmentante:{}, ".format(len(dico_cliques_augmentante)) \
         +" nbre_elts_pi1_pi2:{}".format(nbre_elts_pi1_pi2))
    logger.debug(" * compression_sommet : " \
            +" sommet_z ={}, ".format(sommet_z) \
            +" dico_C1_C2_S1:{}, ".format(len(dico_C1_C2_S1)) \
            +" dico_cliques_augmentante:{}, ".format(len(dico_cliques_augmentante)) \
            +" nbre_elts_pi1_pi2:{}".format(nbre_elts_pi1_pi2)
            )
    
    #TODO NOK: A refaire car pas de melange de solution """
    ##################33 test combinaision de dico
#    """
#    dico_C1_C2_S1[(cpt, C1, C2, S1_new)] = {
#                      "cliques_contratables":,"S1":,"clique_possible": }
#    dico_cliques_augmentante[(cpt, clique_voisine,\
#                             clique_dependante)] = {
#                              "cliq":, "voisine":,
#                              "dependante":,"sommet_z":}  
#    """                                        
#    for dico_c1c2s1_augm in it.islice(map(dict, it.product(dico_C1_C2_S1.items(), 
#                                         dico_cliques_augmentante.items())),
#                              nbre_elts_pi1_pi2):  
    ##################33 test combinaision de dico
    if not dico_C1_C2_S1 and not dico_cliques_augmentante :
        logger.debug(" * compression_sommet : *** dico_C1_C2_S1, " \
                     +"dico_cliques_augmentante VIDES ")
        dico_sommets_non_corriges = dict();
        dico_sommets_corriges = dict();
        for id_sommet_1, sommet_1 in enumerate(sommets_a_corriger):
            dico_sommets_non_corriges[id_sommet_1] = sommet_1;
                
        dico_p1_p2_ps[cpt_prod_cartesien] = {
                    "id_sommet_1": id_sommet_z,
                    "sommet_1": sommet_z,
                    "p1": frozenset(),
                    "p2": frozenset(),
                    "ps": frozenset(),
                    "voisine": frozenset(),
                    "dependante": frozenset(),
                    "contractable1": frozenset(),
                    "contractable2": frozenset(),
                    "S1": frozenset(),
                    "S_z": s_z,
                    "aretes_ajoutees_p1": frozenset(),
                    "aretes_ajoutees_p2": frozenset(),
                    "aretes_supprimees_ps": frozenset(),
                    "aretes_Ec_new": args["aretes_Ec"],
                    "C_new": args["C"],
                    "sommets_corriges": dico_sommets_corriges,
                    "sommets_non_corriges": dico_sommets_non_corriges,
                    "sommets_par_cliqs_avec_aretes_new": 
                        args["sommets_par_cliqs_avec_aretes"]
                    }
        
    elif not dico_C1_C2_S1 and dico_cliques_augmentante :
        logger.debug(" * compression_sommet : *** dico_C1_C2_S1 VIDE, " \
                     + "dico_cliques_augmentante NON VIDE ")
        for k_cpt_vois_depend, val_cpt_vois_depend in dico_cliques_augmentante.items():
            cpt_prod_cartesien += 1;
            p1 = val_cpt_vois_depend["voisine"].union(
                                val_cpt_vois_depend["dependante"].union(
                                frozenset({sommet_z})));
            p2 = frozenset();
            gamma_z = args["dico_gamma_sommets"][sommet_z][1];
            ps = gamma_z - p1.intersection(gamma_z);
            aretes_ps = set( frozenset((sommet_z, sommet_ps)) 
                             for sommet_ps in ps);
            aretes_p1 = fct_aux.determiner_aretes_cliques(p1);
            aretes_ajoutees_p1 = fct_aux.aretes_differente(args["aretes_Ec"], 
                                                   aretes_p1);
            aretes_Ec_new = set(args["aretes_Ec"]).union(aretes_ajoutees_p1);
            
            C_new = set(args["C"].copy());
            ens_cliq_a_supprimer = set();
            ens_cliq_a_supprimer.add(val_cpt_vois_depend["voisine"]);
            ens_cliq_a_supprimer.add(val_cpt_vois_depend["dependante"]);
            for cliq in val_cpt_vois_depend["cliques_suppl_contractables"]: 
                ens_cliq_a_supprimer.add(cliq);
            C_new.difference_update(ens_cliq_a_supprimer);
            
            C_new.add( p1 );
            dico_sommets_corriges = dict();
            dico_sommets_non_corriges = dict();
            dico_sommets_par_cliqs_new = dict();
            C_new, aretes_Ec,\
            dico_sommets_corriges, \
            dico_sommets_non_corriges, \
            dico_sommets_par_cliqs_new = \
                mise_a_jour_aretes_cliques(
                        C_new.copy(), 
                        aretes_Ec_new.copy(),
                        aretes_ps,
                        sommets_a_corriger.copy(),
                        args["sommets_par_cliqs_avec_aretes"])
                
            dico_p1_p2_ps[cpt_prod_cartesien] = {
                        "id_sommet_1": id_sommet_z,
                        "sommet_1": sommet_z,
                        "p1": p1,
                        "p2": p2,
                        "ps": ps,
                        "voisine": val_cpt_vois_depend["voisine"],
                        "dependante": val_cpt_vois_depend["dependante"],
                        "contractable1": frozenset(),
                        "contractable2": frozenset(),
                        "S1": frozenset(),
                        "S_z": s_z,
                        "aretes_ajoutees_p1": aretes_ajoutees_p1,
                        "aretes_ajoutees_p2": list(),
                        "aretes_supprimees_ps": aretes_ps,
                        "aretes_Ec_new": aretes_Ec_new,
                        "C_new": C_new,
                        "sommets_corriges": dico_sommets_corriges,
                        "sommets_non_corriges": dico_sommets_non_corriges,
                        "sommets_par_cliqs_avec_aretes_new": 
                            dico_sommets_par_cliqs_new
                        }
            
    elif dico_C1_C2_S1 and not dico_cliques_augmentante:
        logger.debug(" * compression_sommet : *** dico_C1_C2_S1 NON VIDE, " \
                     + "dico_cliques_augmentante VIDE ")
        for k_c1_c2_s1, val_cpt_c1_c2_s1 in dico_C1_C2_S1.items():
            cpt_prod_cartesien += 1;
            p1 = val_cpt_c1_c2_s1["clique_possible"];
            p2 = frozenset();
            gamma_z = args["dico_gamma_sommets"][sommet_z][1];
            ps = gamma_z - p1.intersection(gamma_z);
            aretes_ps = set( frozenset((sommet_z, sommet_ps)) 
                             for sommet_ps in ps);
            aretes_p1 = fct_aux.determiner_aretes_cliques(p1);
            aretes_ajoutees_p1 = fct_aux.aretes_differente(args["aretes_Ec"], 
                                                 aretes_p1);
            
            aretes_p2 = fct_aux.determiner_aretes_cliques(p2);
            aretes_ajoutees_p2 = fct_aux.aretes_differente(args["aretes_Ec"], 
                                                 aretes_p2);
                                                
            aretes_Ec_new = set(args["aretes_Ec"]).union(
                            aretes_ajoutees_p1.union(aretes_ajoutees_p2));
            
            C_new = set(args["C"].copy());
            
            ens_cliq_a_supprimer = set();
            ens_cliq_a_supprimer.add(
                    val_cpt_c1_c2_s1["cliques_contratables"][0]);
            ens_cliq_a_supprimer.add(
                    val_cpt_c1_c2_s1["cliques_contratables"][1]);
            for cliq in val_cpt_c1_c2_s1["cliques_suppl_contractables"] :
                ens_cliq_a_supprimer.add(cliq);
            C_new.difference_update(ens_cliq_a_supprimer);
            
            C_new.add( p1 );
            dico_sommets_corriges = dict();
            dico_sommets_non_corriges = dict();
            dico_sommets_par_cliqs_new = dict();
            C_new, aretes_Ec,\
            dico_sommets_corriges, \
            dico_sommets_non_corriges, \
            dico_sommets_par_cliqs_new = \
                mise_a_jour_aretes_cliques(
                        C_new.copy(), 
                        aretes_Ec_new.copy(), 
                        aretes_ps,\
                        sommets_a_corriger.copy(), \
                        args["sommets_par_cliqs_avec_aretes"])
                
            dico_p1_p2_ps[cpt_prod_cartesien] = {
                        "id_sommet_1": id_sommet_z,
                        "sommet_1": sommet_z,
                        "p1": val_cpt_c1_c2_s1["clique_possible"],
                        "p2": frozenset(),
                        "ps": ps,
                        "voisine": frozenset(),
                        "dependante": frozenset(),
                        "contractable1": 
                            val_cpt_c1_c2_s1["cliques_contratables"][0],
                        "contractable2": 
                            val_cpt_c1_c2_s1["cliques_contratables"][1],
                        "S1": val_cpt_c1_c2_s1["S1"],
                        "S_z": s_z,
                        "aretes_ajoutees_p1": aretes_ajoutees_p1,
                        "aretes_ajoutees_p2": aretes_ajoutees_p2,
                        "aretes_supprimees_ps": aretes_ps,
                        "aretes_Ec_new": aretes_Ec_new,
                        "C_new": C_new,
                        "sommets_corriges": dico_sommets_corriges,
                        "sommets_non_corriges": dico_sommets_non_corriges,
                        "sommets_par_cliqs_avec_aretes_new": 
                            dico_sommets_par_cliqs_new
                        }
    else:
        logger.debug(" * compression_sommet : *** dico_C1_C2_S1, " \
                     + "dico_cliques_augmentante NON VIDE ")
        for k_c1_c2_s1, val_cpt_c1_c2_s1 in dico_C1_C2_S1.items():
            for k_cpt_vois_depend, val_cpt_vois_depend in dico_cliques_augmentante.items():                                        
                cpt_prod_cartesien += 1;
                
                inter_p1_p2 = val_cpt_c1_c2_s1["clique_possible"].intersection(
                                k_cpt_vois_depend[1].union(k_cpt_vois_depend[2])
                                )
                logger.debug(" * compression_sommet : *** "+\
                             "cpt_prod_cart={},".format(cpt_prod_cartesien)+\
                             "inter_p1_p2={},".format(len(inter_p1_p2)))
#                print(" ***{} inter_p1_p2={},".format(
#                                            cpt_prod_cartesien,inter_p1_p2) \
#                      +"cliq_possible={},".format(
#                                          val_cpt_c1_c2_s1["clique_possible"]) \
#                      +"vois_dep={}".format(k_cpt_vois_depend[1].union(
#                                                  k_cpt_vois_depend[2]
#                                                  )
#                      ))
                
                if len(inter_p1_p2) <= 1 \
                    and inter_p1_p2 == frozenset({sommet_z}):
                    
                    p1 = val_cpt_c1_c2_s1["clique_possible"];
                    p2 = val_cpt_vois_depend["voisine"].union(
                                val_cpt_vois_depend["dependante"].union(
                                frozenset({sommet_z})))
                    gamma_z = args["dico_gamma_sommets"][sommet_z][1];
                    ps = gamma_z - p1.intersection(gamma_z).union(
                                    p2.intersection(gamma_z)
                                    );
                    aretes_ps = set( frozenset((sommet_z, sommet_ps)) 
                                     for sommet_ps in ps);
                    
                    aretes_p1 = fct_aux.determiner_aretes_cliques(p1);
                    aretes_ajoutees_p1 = fct_aux.aretes_differente(
                                                args["aretes_Ec"], 
                                                aretes_p1);
                    
                    aretes_p2 = fct_aux.determiner_aretes_cliques(p2);
                    aretes_ajoutees_p2 = fct_aux.aretes_differente(
                                                args["aretes_Ec"], 
                                                aretes_p2);
                                                        
                    aretes_Ec_new = set(args["aretes_Ec"]).union(
                                    aretes_ajoutees_p1.union(aretes_ajoutees_p2));
                    
                    C_new = set(args["C"].copy());

                    ens_cliq_a_supprimer = set();
                    ens_cliq_a_supprimer.add(
                            val_cpt_c1_c2_s1["cliques_contratables"][0]);
                    ens_cliq_a_supprimer.add(
                            val_cpt_c1_c2_s1["cliques_contratables"][1]);
                    ens_cliq_a_supprimer.add(
                            val_cpt_vois_depend["voisine"]);
                    ens_cliq_a_supprimer.add(
                            val_cpt_vois_depend["dependante"]);
                    for cliq in [item for subitem \
                                 in [val_cpt_c1_c2_s1["cliques_suppl_contractables"], 
                                 val_cpt_vois_depend["cliques_suppl_contractables"]] \
                                 for item in subitem]:
                        ens_cliq_a_supprimer.add(cliq);
                    C_new.difference_update(ens_cliq_a_supprimer);
                    
                    C_new.add( p1 );
                    C_new.add( p2 );
                    dico_sommets_corriges = dict();
                    dico_sommets_non_corriges = dict();
                    dico_sommets_par_cliqs_new = dict();
                    C_new, aretes_Ec,\
                    dico_sommets_corriges, \
                    dico_sommets_non_corriges, \
                    dico_sommets_par_cliqs_new = \
                        mise_a_jour_aretes_cliques(
                                C_new.copy(), 
                                aretes_Ec_new.copy(), 
                                aretes_ps,
                                sommets_a_corriger.copy(), 
                                args["sommets_par_cliqs_avec_aretes"])
                    
                    dico_p1_p2_ps[cpt_prod_cartesien] = {
                        "id_sommet_1": id_sommet_z,
                        "sommet_1": sommet_z,
                        "p1": val_cpt_c1_c2_s1["clique_possible"],
                        "p2": val_cpt_vois_depend["voisine"].union(
                                val_cpt_vois_depend["dependante"].union(
                                frozenset({sommet_z}))),
                        "ps": ps,
                        "voisine": val_cpt_vois_depend["voisine"],
                        "dependante": val_cpt_vois_depend["dependante"],
                        "contractable1": 
                            val_cpt_c1_c2_s1["cliques_contratables"][0],
                        "contractable2": 
                            val_cpt_c1_c2_s1["cliques_contratables"][1],
                        "S1": val_cpt_c1_c2_s1["S1"],
                        "S_z": s_z,
                        "aretes_ajoutees_p1": aretes_ajoutees_p1,
                        "aretes_ajoutees_p2": aretes_ajoutees_p2,
                        "aretes_supprimees_ps": aretes_ps,
                        "aretes_Ec_new": aretes_Ec_new,
                        "C_new": C_new,
                        "sommets_corriges": dico_sommets_corriges,
                        "sommets_non_corriges": dico_sommets_non_corriges,
                        "sommets_par_cliqs_avec_aretes_new": 
                            dico_sommets_par_cliqs_new
                        }
                else:
                    # TODO traiter le cas ou il ya  inter_p1_p2 != frozenset({sommet_z})
#                    p1 = val_cpt_c1_c2_s1["clique_possible"];
#                    p2 = val_cpt_vois_depend["voisine"].union(
#                                val_cpt_vois_depend["dependante"].union(
#                                frozenset({sommet_z})))
#                    gamma_z = args["dico_gamma_sommets"][sommet_z][1];
#                    ps = gamma_z - p1.intersection(gamma_z).union(
#                                    p2.intersection(gamma_z));
#                    aretes_ps = set( frozenset((sommet_z, sommet_ps)) 
#                                            for sommet_ps in ps)
                    pass
                if cpt_prod_cartesien >= nbre_elts_pi1_pi2:
                    break;
            if cpt_prod_cartesien >= nbre_elts_pi1_pi2:
                break;
#    print("@@@cpt_prod_cartesien={}, dico_C1_C2_S1={}, dico_cliques_augmentante={}".\
#          format(cpt_prod_cartesien, len(dico_C1_C2_S1), len(dico_cliques_augmentante)))   
    logger.debug(" * compression_sommet : *** fin compression_sommet, "+\
                 " dico_p1_p2_ps:{}, ".format(len(dico_p1_p2_ps)))
    return dico_p1_p2_ps;
###############################################################################
#               compression d'un sommet => fin
###############################################################################

###############################################################################
#                critere selection  compression => debut
###############################################################################
def rechercher_min_max(liste_tuples, critere):
    """ 
    retourne la tuple (min, max)
    """
    min_c1 = np.inf;
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
        
    
def critere_C2_C1_local(dico_compression, args):
    """ 
    selectionner le dico selon C2 puis C1 parmi les compressions possibles 
        sommet a corriger (sommet a -1)
    
    C2 : le maximum de sommets corriges
        * choisir le sommet a -1 qui corrige le max de sommets a -1 possibles
    C1 : le minimum d'aretes corriges
    
    dico_compression : dictionnaire contenant les compression (p1,p2,ps) du 
                        sommet sommet_z
    """
    
    max_c2 = 0;
    min_c1 = np.inf;
    dico_c1_c2 = dict();
    
    print("dico_compression={}".format( len(dico_compression) ))
    if not dico_compression :
        return min_c1, max_c2, [];
        
    # definition de C2
    if args["critere_selection_compression"] == "voisins_corriges":             # C2
        for cpt_prod_cartesien, dico_p1_p2_ps in dico_compression.items():
            if len(dico_p1_p2_ps["sommets_corriges"]) >= max_c2:
                max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
                nbre_aretes_corriges = \
                                len(dico_p1_p2_ps["aretes_ajoutees_p1"]) \
                                + len(dico_p1_p2_ps["aretes_ajoutees_p2"]) \
                                + len(dico_p1_p2_ps["aretes_supprimees_ps"]);
                min_c1 = nbre_aretes_corriges if min_c1 >= nbre_aretes_corriges \
                                                else min_c1;
                if (min_c1,max_c2) not in dico_c1_c2:
                    dico_c1_c2[(min_c1,max_c2)] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2[(min_c1,max_c2)].append(dico_p1_p2_ps);
        print("@@C2 max_c2={}, min_c1={}, dico_c1_c2={}".format(
              max_c2, min_c1,len(dico_c1_c2[(min_c1,max_c2)])))
        
    # definition de C1
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees":     # C1
        for cpt_prod_cartesien, dico_p1_p2_ps in dico_compression.items():
            nbre_aretes_corriges = len(dico_p1_p2_ps["aretes_ajoutees_p1"]) \
                                   + len(dico_p1_p2_ps["aretes_ajoutees_p2"]) \
                                   + len(dico_p1_p2_ps["aretes_supprimees_ps"]);
            min_c1 = nbre_aretes_corriges if min_c1 >= nbre_aretes_corriges \
                                          else min_c1;
            max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
            if (min_c1,max_c2) not in dico_c1_c2:
                dico_c1_c2[(min_c1,max_c2)] = [dico_p1_p2_ps];
            else:
                dico_c1_c2[(min_c1,max_c2)].append(dico_p1_p2_ps);
        print("@@C1 min_c1={}, max_c2={}, dico_c1_c2={}".format(
              min_c1, max_c2, len(dico_c1_c2[(min_c1,max_c2)])))
        
    # definition de C2 puis de C1
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees": # C2_C1
        for cpt_prod_cartesien, dico_p1_p2_ps in dico_compression.items():
            if len(dico_p1_p2_ps["sommets_corriges"]) >= max_c2:
                max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
                nbre_aretes_corriges = \
                                    len(dico_p1_p2_ps["aretes_ajoutees_p1"]) \
                                    + len(dico_p1_p2_ps["aretes_ajoutees_p2"]) \
                                    + len(dico_p1_p2_ps["aretes_supprimees_ps"]);
                min_c1 = nbre_aretes_corriges if min_c1 >= nbre_aretes_corriges \
                                                else min_c1;
                if (min_c1,max_c2) not in dico_c1_c2:
                    dico_c1_c2[(min_c1,max_c2)] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2[(min_c1,max_c2)].append(dico_p1_p2_ps);
        print("@@C1_C2 min_c1={}, max_c2={}, dico_c1_c2={}".format(min_c1, 
              max_c2, len(dico_c1_c2[(min_c1,max_c2)])))
        
    if not dico_c1_c2:
        return min_c1, max_c2, [];
    else:
        return min_c1, max_c2, dico_c1_c2[(min_c1,max_c2)];
        
def critere_C2_C1_global(dico_compression, args):
    """ 
    recherche la compression optimale parmi tous les sommets a corriger 
        selon les criteres C1 et C2.
        
    C2 : le maximum de sommets corriges
        * choisir le sommet a -1 qui corrige le max de sommets a -1 possibles
    C1 : le minimum d'aretes corriges    
    """
    
    max_c2_global = 0;
    min_c1_global = np.inf;
    dico_c1_c2_global = dict();
    cle_min_max_c2 = None;
    
#    if not dico_compression:
#        return min_c1_global, \
#                max_c2_global, \
#                dico_c1_c2_global[cle_min_max_c2][numero_sol_c1_c2];
    
    # selection C2
    # je cherche le min local de c1 pour tous les sommets a corriger
    # parmi les min locaux, je cherche le max global de c2
    # une fois la liste des (min_global,max_global), je prends le 1er element.
    if args["critere_selection_compression"] == "voisins_corriges":            # C2
        for id_sommet_z, dicos_p1_p2_ps in dico_compression.items():
            # selection de dico selon C1
            min_c1_local = dicos_p1_p2_ps[0];
            max_c2_local = dicos_p1_p2_ps[1]
            for dico_p1_p2_ps in dicos_p1_p2_ps[2] :
                nbre_aretes_corriges = \
                                    len(dico_p1_p2_ps["aretes_ajoutees_p1"]) \
                                    + len(dico_p1_p2_ps["aretes_ajoutees_p2"]) \
                                    + len(dico_p1_p2_ps["aretes_supprimees_ps"]);
                min_c1_local = nbre_aretes_corriges \
                                if min_c1_local >= nbre_aretes_corriges \
                                else min_c1_local;
                if (min_c1_local, max_c2_local) not in dico_c1_c2_global:
                    dico_c1_c2_global[(min_c1_local, 
                                       max_c2_local)] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2_global[(min_c1_local, 
                                       max_c2_local)].append(dico_p1_p2_ps);
                                     
        # selection selon C2
        cle_min_max_c2 = rechercher_min_max(dico_c1_c2_global.keys(), "C2");
        
    
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees":   # C1
        for id_sommet_z, dicos_p1_p2_ps in dico_compression.items():
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
        print("---> min_max_s={}".format(dico_c1_c2_global.keys()))
        cle_min_max_c2 = rechercher_min_max(dico_c1_c2_global.keys(), "C1");
        
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees": # C2_C1
        for id_sommet_z, dicos_p1_p2_ps in dico_compression.items():
            min_c1_local = dicos_p1_p2_ps[0]; #np.inf
            max_c2_local = dicos_p1_p2_ps[1]; #0
            for dico_p1_p2_ps in dicos_p1_p2_ps[2] :
                nbre_sommets_corriges = len(dico_p1_p2_ps["sommets_corriges"]);
                max_c2_local = nbre_sommets_corriges \
                                if nbre_sommets_corriges > max_c2_local \
                                else max_c2_local;
                nbre_aretes_corriges = len(dico_p1_p2_ps["aretes_ajoutees_p1"]) \
                                      + len(dico_p1_p2_ps["aretes_ajoutees_p2"]) \
                                      + len(dico_p1_p2_ps["aretes_supprimees_ps"]);
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
    return min_c1_global, \
            max_c2_global, \
            dico_c1_c2_global[cle_min_max_c2][numero_sol_c1_c2];
###############################################################################
#               critere selection  compression  => fin
###############################################################################

###############################################################################
#               application de la compression => debut
###############################################################################
def appliquer_correction(dico_sol_C2_C1, sommets_a_corriger, args):
    """ 
    appliquer la compression choisie dans le graphe.
    """
    C = list();
    C = dico_sol_C2_C1["C_new"];
    aretes_Ec = list();
    aretes_Ec = dico_sol_C2_C1["aretes_Ec_new"];
    
    id_sommets_1 = list(dico_sol_C2_C1["sommets_corriges"].keys());
#    print("****1 id_sommets_1={}".format(id_sommets_1))
    
    id_sommets_1.append(dico_sol_C2_C1["id_sommet_1"]);
#    print("****2 id_sommets_1={}".format(id_sommets_1))
    
    sommets_corriges = dico_sol_C2_C1["sommets_corriges"].values();
#    print("**** sommets_corriges={},sommet_1={}".format(sommets_corriges,
#          dico_sol_C2_C1["sommet_1"]))
#    print("****1 avant supp sommets_a_corriger={}".format(sommets_a_corriger))
    
    sommets_a_corriger = np.delete(sommets_a_corriger, id_sommets_1).tolist();
#    print("****2 apres supp sommets_a_corriger={}".format(sommets_a_corriger))
    
    if set(sommets_a_corriger).intersection(set(sommets_corriges)) :
        print("---ERROR : sommets {} suppression : NOK -----".
              format(sommets_corriges))

    return C, aretes_Ec, sommets_a_corriger;
###############################################################################
#               application de la compression => fin
###############################################################################

###############################################################################
#               correction des sommets_1 avec la selection avec critere
#                            => debut
###############################################################################
def correction_avec_critere(sommets_a_corriger, 
                        dico_correction,
                        dico_parametres_new):
    """
    corrige un graphe de correlation en ajoutant ou supprimant des aretes
    avec la selection specifique des sommets a corriger.
    La selection depend du critere critere_selection_compression
    
    # TODO a corriger : presente des erreurs au niveau des criteres globaux
                        et locaux avec une cle a min.
    """
    logger = logging.getLogger('correction_aleatoire_graphe_correlation');
    dico_sommets_corriges = dict();
    
    logger.debug(" * critere_selection_compression : {}".\
                     format(dico_parametres_new["critere_selection_compression"]))
    logger.debug(" * mode_correction : {}".\
                     format(dico_parametres_new["mode_select_noeuds_1"]))
    
    cpt_noeud = 0;
    dico_correction["number_items_pi1_pi2"] = \
                                dico_parametres_new["number_items_pi1_pi2"];
    
    print("sommets_a_corriger={} = {}".format(
            len(sommets_a_corriger), sommets_a_corriger));
            
    while(sommets_a_corriger):
        dico_compression = dict();
        
        for id_sommet_1, sommet_1 in enumerate(sommets_a_corriger):
            cliques_sommet_1 = \
                    dico_correction["sommets_par_cliqs_avec_aretes"][sommet_1];
            print(" \n* corr_graphe sommet_1:{}, ".format(sommet_1) \
                    +" cliques_sommets_1:{}".format(len(cliques_sommet_1))
                    )
            logger.debug(" * corr_graphe sommet_1:{}, ".format(sommet_1) \
                    +" cliques_sommets_1:{}".format(len(cliques_sommet_1))
                        )
            
            dico_p1_p2_ps = dict();
            dico_p1_p2_ps = compression_sommet(id_sommet_1,
                                               sommet_1,
                                               sommets_a_corriger,
                                               cliques_sommet_1,
                                               dico_correction);
                                               
            dico_compression[(id_sommet_1,sommet_1)] = critere_C2_C1_local(
                                                        dico_p1_p2_ps,
                                                        dico_parametres_new);
                    
            print(" ** cal_p1_p2_ps sommet_1:{},".format(sommet_1) \
            + " id_sommet_1:{},".format(id_sommet_1) \
            + " min_c1:{},".format(dico_compression[(id_sommet_1,sommet_1)][0]) \
            + " max_c2:{},".format(dico_compression[(id_sommet_1,sommet_1)][1]) \
            + " dico_c1_c2:{}".format(
                        len(dico_compression[(id_sommet_1,sommet_1)][2]))
            )
            
            logger.debug(" * cal_p1_p2_ps sommet_1:{},".format(sommet_1) \
             +" id_sommet_1:{},".format(id_sommet_1) \
             +" min_c1:{},".format(dico_compression[(id_sommet_1,sommet_1)][0]) \
             +" max_c2:{},".format(dico_compression[(id_sommet_1,sommet_1)][1]) \
             +" dico_c1_c2:{}".format(len(dico_compression[(id_sommet_1,sommet_1)][2]))
             )
            pass # end for id_sommet_1, sommet_1
            
        dico_sol_C2_C1 = dict();
        min_c1 = 0; max_c2 = 0;
        min_c1, max_c2, dico_sol_C2_C1 = critere_C2_C1_global(
                                        dico_compression,
                                        dico_parametres_new)               # C2 : nombre maximum de voisins corriges par un sommet, C1 : nombre minimum d'aretes a corriger au voisinage d'un sommet  
        
        logger.debug(
        " * choix_sommet : sommet_1:{}".format(dico_sol_C2_C1["sommet_1"]) \
        +" min_c1:{}".format(min_c1) \
        +" max_c2:{}".format(max_c2) \
        +" aretes_ajoutees_p1:{}".format(dico_sol_C2_C1["aretes_ajoutees_p1"]) \
        +" aretes_ajoutees_p2:{}".format(dico_sol_C2_C1["aretes_ajoutees_p2"]) \
        +" aretes_supprimees_ps:{}".format(dico_sol_C2_C1["aretes_supprimees_ps"]) \
        +" sommets_corriges:{}".format(dico_sol_C2_C1["sommets_corriges"]) \
        +" sommets_non_corriges:{}".format(dico_sol_C2_C1["sommets_non_corriges"]) 
                     )
        
        print("AVANT => sommets_a_corriger={}={}".format(
                len(sommets_a_corriger), sommets_a_corriger) \
            +" sommet_1_choisi = {}".format(dico_sol_C2_C1["sommet_1"])
            )
        C, aretes_Ec, sommets_a_corriger = appliquer_correction(
                                            dico_sol_C2_C1,
                                            sommets_a_corriger,
                                            dico_correction)
        print("APRES => sommets_a_corriger={}={}".format(
                len(sommets_a_corriger), sommets_a_corriger));
        
        logger.debug(
        " * appli_correction: "
        +" C_old:{}".format(len(dico_correction["C"])) \
        +" C:{}".format(len(C)) \
        +" aretes_Ec_old:{}".format(len(dico_correction["aretes_Ec"])) \
        +" aretes_Ec:{}".format(len(aretes_Ec)) \
        +" sommets_a_corriger={}".format(len(sommets_a_corriger))
                     )
        
        for sommet, cliques in \
            dico_sol_C2_C1["sommets_par_cliqs_avec_aretes_new"].items():
            logger.debug(" * appli_correction: "
                         +"sommet_par_cliques {}={}".format(
                                 sommet, cliques));
        
        dico_correction["C"] = C;
        dico_correction["aretes_cliques"] = \
                                    fct_aux.determiner_aretes_cliques(C);
        dico_correction["aretes_Ec"] = aretes_Ec;
        cout_T = {"aretes_ajoutees_p1": dico_sol_C2_C1["aretes_ajoutees_p1"],
                  "aretes_ajoutees_p2": dico_sol_C2_C1["aretes_ajoutees_p2"],
                  "aretes_supprimees": dico_sol_C2_C1["aretes_supprimees_ps"],
                  "min_c1": min_c1, 
                  "max_c2": max_c2};
        cpt_noeud += 1;
        dico_sommets_corriges[(cpt_noeud, dico_sol_C2_C1["sommet_1"])] = {
                    "compression_p1": dico_sol_C2_C1["p1"],
                    "compression_p2": dico_sol_C2_C1["p2"],
                    "compression_ps": dico_sol_C2_C1["ps"],
                    "sommets_corriges": dico_sol_C2_C1["sommets_corriges"], # voisins_corriges = {"id_voisin_ds_sommets_a_corriger":voisin}
                    "cout_T": cout_T
                    }
        
        # mettre a jour les cliques couvrants les sommets.
        dico_correction["sommets_par_cliqs_avec_aretes"] = \
                    dico_sol_C2_C1["sommets_par_cliqs_avec_aretes_new"];
        pass # end while sommets_a_corriger
        
    return dico_correction, dico_sommets_corriges;
    pass
###############################################################################
#               correction des sommets_1 avec la selection avec critere  
#                               => fin
###############################################################################
    
###############################################################################
#               correction aleatoire des sommets_1 => debut
###############################################################################
def correction_aleatoire(sommets_a_corriger, 
                        dico_correction,
                        dico_parametres_new):
    """
    corrige un graphe de correlation en ajoutant ou supprimant des aretes
    avec la selection aleatoire des sommets a corriger.
    """
    logger = logging.getLogger('correction_aleatoire_graphe_correlation');
    dico_sommets_corriges = dict();
    
    logger.debug(" * critere_selection_compression : {}".\
                     format(dico_parametres_new["critere_selection_compression"]))
    logger.debug(" * mode_correction : {}".\
                     format(dico_parametres_new["mode_select_noeuds_1"]))
    
    cpt_noeud = 0;
    dico_correction["number_items_pi1_pi2"] = \
                                dico_parametres_new["number_items_pi1_pi2"];
    
    print("sommets_a_corriger={} = {}".format(
            len(sommets_a_corriger), sommets_a_corriger));
            
    dico_compression = dict();
    
    for id_sommet_1, sommet_1 in enumerate(sommets_a_corriger):
        cliques_sommet_1 = \
                dico_correction["sommets_par_cliqs_avec_aretes"][sommet_1];
        print(" \n* corr_graphe sommet_1:{}, ".format(sommet_1) \
                +" cliques_sommets_1:{}".format(len(cliques_sommet_1))
                )
        logger.debug(" * corr_graphe sommet_1:{}, ".format(sommet_1) \
                +" cliques_sommets_1:{}".format(len(cliques_sommet_1))
                    )
        
        dico_p1_p2_ps = dict();
        dico_p1_p2_ps = compression_sommet(id_sommet_1,
                                           sommet_1,
                                           sommets_a_corriger,
                                           cliques_sommet_1,
                                           dico_correction);
                                           
        dico_compression[(id_sommet_1,sommet_1)] = critere_C2_C1_local(
                                                    dico_p1_p2_ps,
                                                    dico_parametres_new);
                
        print(" ** cal_p1_p2_ps sommet_1:{},".format(sommet_1) \
        + " id_sommet_1:{},".format(id_sommet_1) \
        + " min_c1:{},".format(dico_compression[(id_sommet_1,sommet_1)][0]) \
        + " max_c2:{},".format(dico_compression[(id_sommet_1,sommet_1)][1]) \
        + " dico_c1_c2:{}".format(
                    len(dico_compression[(id_sommet_1,sommet_1)][2]))
        )
        
        logger.debug(" * cal_p1_p2_ps sommet_1:{},".format(sommet_1) \
         +" id_sommet_1:{},".format(id_sommet_1) \
         +" min_c1:{},".format(dico_compression[(id_sommet_1,sommet_1)][0]) \
         +" max_c2:{},".format(dico_compression[(id_sommet_1,sommet_1)][1]) \
         +" dico_c1_c2:{}".format(len(dico_compression[(id_sommet_1,sommet_1)][2]))
         )
        pass # end for id_sommet_1, sommet_1
        
        dico_sol_C2_C1 = dict();
        min_c1 = 0; max_c2 = 0;
        min_c1, max_c2, dico_sol_C2_C1 = critere_C2_C1_global(
                                        dico_compression,
                                        dico_parametres_new)                    # C2 : nombre maximum de voisins corriges par un sommet, C1 : nombre minimum d'aretes a corriger au voisinage d'un sommet  
        
        logger.debug(
        " * choix_sommet : sommet_1:{}".format(dico_sol_C2_C1["sommet_1"]) \
        +" min_c1:{}".format(min_c1) \
        +" max_c2:{}".format(max_c2) \
        +" aretes_ajoutees_p1:{}".format(dico_sol_C2_C1["aretes_ajoutees_p1"]) \
        +" aretes_ajoutees_p2:{}".format(dico_sol_C2_C1["aretes_ajoutees_p2"]) \
        +" aretes_supprimees_ps:{}".format(dico_sol_C2_C1["aretes_supprimees_ps"]) \
        +" sommets_corriges:{}".format(dico_sol_C2_C1["sommets_corriges"]) \
        +" sommets_non_corriges:{}".format(dico_sol_C2_C1["sommets_non_corriges"]) 
                     )
        
        print("AVANT => sommets_a_corriger={}={}".format(
                len(sommets_a_corriger), sommets_a_corriger) \
            +" sommet_1_choisi = {}".format(dico_sol_C2_C1["sommet_1"])
            )
        C, aretes_Ec, sommets_a_corriger = appliquer_correction(
                                            dico_sol_C2_C1,
                                            sommets_a_corriger,
                                            dico_correction)
        print("APRES => sommets_a_corriger={}={}".format(
                len(sommets_a_corriger), sommets_a_corriger));
        
        logger.debug(
        " * appli_correction: "
        +" C_old:{}".format(len(dico_correction["C"])) \
        +" C:{}".format(len(C)) \
        +" aretes_Ec_old:{}".format(len(dico_correction["aretes_Ec"])) \
        +" aretes_Ec:{}".format(len(aretes_Ec)) \
        +" sommets_a_corriger={}".format(len(sommets_a_corriger))
                     )
        
        for sommet, cliques in \
            dico_sol_C2_C1["sommets_par_cliqs_avec_aretes_new"].items():
            logger.debug(" * appli_correction: "
                         +"sommet_par_cliques {}={}".format(
                                 sommet, cliques));
        
        dico_correction["C"] = C;
        dico_correction["aretes_cliques"] = \
                                    fct_aux.determiner_aretes_cliques(C);
        dico_correction["aretes_Ec"] = aretes_Ec;
        cout_T = {"aretes_ajoutees_p1": dico_sol_C2_C1["aretes_ajoutees_p1"],
                  "aretes_ajoutees_p2": dico_sol_C2_C1["aretes_ajoutees_p2"],
                  "aretes_supprimees": dico_sol_C2_C1["aretes_supprimees_ps"],
                  "min_c1": min_c1, 
                  "max_c2": max_c2};
        cpt_noeud += 1;
        dico_sommets_corriges[(cpt_noeud, dico_sol_C2_C1["sommet_1"])] = {
                    "compression_p1": dico_sol_C2_C1["p1"],
                    "compression_p2": dico_sol_C2_C1["p2"],
                    "compression_ps": dico_sol_C2_C1["ps"],
                    "sommets_corriges": dico_sol_C2_C1["sommets_corriges"], # voisins_corriges = {"id_voisin_ds_sommets_a_corriger":voisin}
                    "cout_T": cout_T
                    };
        
        # mettre a jour les cliques couvrants les sommets.
        dico_correction["sommets_par_cliqs_avec_aretes"] = \
                    dico_sol_C2_C1["sommets_par_cliqs_avec_aretes_new"];
        
    return dico_correction, dico_sommets_corriges;
###############################################################################
#               correction aleatoire des sommets_1  => fin
###############################################################################

###############################################################################
#               algorithme de correction pour k entier naturel => debut
###############################################################################
def correction_cliques_k(dico_correction,
                         dico_parametres_new):
    """
    corriger les cliques selon le mode de selection des sommets a -1
    """
    
    sommets_a_corriger = dico_correction["sommets_a_corriger"];
    dico_sommets_corriges = dict();
    
    dico_correction["aretes_cliques"] = fct_aux.determiner_aretes_cliques(
                                            dico_correction["C"]
                                            )
    print("########### critere_selection_compression={}, mode_select_noeuds_1={}"\
          .format(dico_parametres_new["critere_selection_compression"],
                  dico_parametres_new["mode_select_noeuds_1"]))
    
    if dico_parametres_new["critere_selection_compression"] == "voisins_corriges" \
        and dico_parametres_new["mode_select_noeuds_1"] == "aleatoire":
        dico_correction, dico_sommets_corriges = correction_aleatoire(
                                                    sommets_a_corriger, 
                                                    dico_correction,
                                                    dico_parametres_new);
        
    elif dico_parametres_new["critere_selection_compression"] == \
        "nombre_aretes_corrigees" \
        and dico_parametres_new["mode_select_noeuds_1"] == "aleatoire":
        dico_correction, dico_sommets_corriges = correction_aleatoire(
                                                    sommets_a_corriger, 
                                                    dico_correction,
                                                    dico_parametres_new);
        # TODO presente des erreurs
#        dico_correction, dico_sommets_corriges = correction_avec_critere(
#                                                    sommets_a_corriger, 
#                                                    dico_correction,
#                                                    dico_parametres_new);
        
    elif dico_parametres_new["critere_selection_compression"] == \
        "voisins_nombre_aretes_corrigees" \
        and dico_parametres_new["mode_select_noeuds_1"] == "aleatoire":
        dico_correction, dico_sommets_corriges = correction_aleatoire(
                                                    sommets_a_corriger, 
                                                    dico_correction,
                                                    dico_parametres_new);
                
    return dico_correction, dico_sommets_corriges;
###############################################################################
#               algorithme de correction pour k entier naturel => fin
###############################################################################



























####################################
##COMMENTAIRE ==> PARTIE INACHEVE, TOUT EST COMMENTE
####################################
################################################################################
##               calcul de la compression d un sommet => debut
################################################################################
#def compression_sommet(id_sommet_1,
#                       nom_sommet_1,
#                       noms_sommets_1,
##                       cliques_sommet_1,
#                       cliques_par_nom_sommets,
#                       cliques_couvertures,
#                       aretes_LG_k_alpha_cor,
#                       sommets_LG,
#                       mode_correction,
#                       critere_correction,
#                       number_items_pi1_pi2,
#                       DBG
#                       
#                       id_sommet_1,
#                       sommet_1,
#                       sommets_1,
#                       cliques_par_sommets,
#                       cliques_couv_cor = cliques_couvertures_cor,
#                            aretes_LG_k_alpha_cor = aretes_LG_k_alpha_cor,
#                            mode_select_noeuds_1 = mode_select_noeuds_1,
#                            critere_select_pi1_pi2 = critere_select_pi1_pi2,
#                            number_items_pi1_pi2 = number_items_pi1_pi2,
#                            DBG = DBG):
#    """ retourne la compression d'un sommet sommet_z. 
#    
#    la compression est le triplet (pi1, pi2, ps) dans lequel 
#        * pi1, pi2 sont des cliques qui fusionnent 
#            - des cliques augmentantes C1, C2 ou 
#            - des cliques contractables C1, C2 ou 
#            - un ensemble S1 tel que S1 n'est contenu par aucune clique C1 ou C2
#        * pi1, pi2 sont des augmentations
#        * ps est un ensemble de sommets u tel que (z,u) doit etre supprime de aretes_LG_k_alpha
#        
#    """
#    aretes_cliques = fct_aux.edges_in_cliques(cliques_couvertures);
#    s_z = S_sommet(sommet_z = nom_sommet_1,
#                   gamma_z = sommets_LG[nom_sommet_1].voisins,
#                   aretes_LG_k_alpha = aretes_LG_k_alpha_cor,
#                   cliques_couvertures = cliques_couvertures,
#                  aretes_cliques = aretes_cliques);
#    
#    # determination de C1 = (C_1,C_2) avec C_1, C_2 contratables
#    dico_C1_C2_S1 = dict(); cpt = 0;
#    
##    cliques_sommet_1 = cliques_par_nom_sommets[nom_sommet_1];
#    cliques_sommet_1 = set(cliques_par_nom_sommets[nom_sommet_1]);
#    cliques_voisines_sommet_1 = clique_voisine_sommet_z(
#                                    sommet_z = nom_sommet_1,
#                                    C = cliques_couvertures,
#                                    cliques_sommet_z = cliques_sommet_1)
#    
#    cliques_contractables_s = cliques_contractables(
#                                nom_sommet_z = nom_sommet_1, 
#                                aretes_LG_k_alpha = aretes_LG_k_alpha_cor, 
#                                aretes_cliques = aretes_cliques, 
#                                cliques_sommet_z = cliques_sommet_1, 
#                                cliques_voisines_z = cliques_voisines_sommet_1,
#                                C = cliques_couvertures,
#                                DBG = DBG)
#
#    random.shuffle(cliques_contractables_s)
#    
#    logger.debug(
#        "****** compres =>" \
#        +" nom_sommet_z={}".format(nom_sommet_1) \
#        +" cliques_contractables_s={}".format(len(cliques_contractables_s)) \
#        +" cliques_voisines_sommet_1 = {}".format(len(cliques_voisines_sommet_1))
#        ) if DBG else None;
#    
#    for clique_C1_C2_Cx in cliques_contractables_s[:MAX_NUMBER]:
#        # construction de dico_C1_C2_S1
#        #        dico_C1_C2_S1[(cpt, (C1,C2,...), (C3,C4,...))] = {
#        #          "cliques_contratables_1":(C1, C2),
#        #          "cliques_contratables_2":(C3, C4),
#        #          "clique_possible_1": ,
#        #          "clique_possible_2": ,
#        #                   }
#        dico_C1_C2_S1[(cpt, clique_C1_C2_Cx, frozenset())] = {
#                             "cliques_contractables_1" : clique_C1_C2_Cx,
#                             "cliques_contractables_2" : frozenset(),
#                             "clique_possible_1" : \
#                                 frozenset.union(
#                                            *clique_C1_C2_Cx).union(
#                                                    frozenset({nom_sommet_1})
#                                                                ),
#                             "clique_possible_2" : frozenset()
#                            }
#        cpt += 1;
#        
#    ## *chercher les paires de cliques contractables tel que 
#    ## *  |contr1 \cap contr2 |= 1
#    logger.debug("****** compres => Avant " \
#                 +" nom_sommet_z={}".format(nom_sommet_1) \
#                 +" dico_C1_C2_S1={}".format(len(dico_C1_C2_S1))
#                ) if DBG else None;      
#          
#    for clique_p1_p2 in it.combinations(cliques_contractables_s[:MAX_NUMBER], 2):
#        clique_p1 = frozenset.union(*clique_p1_p2[0]);
#        clique_p2 = frozenset.union(*clique_p1_p2[1]);
#        if cpt > MAX_NUMBER:
#            break;
#        if len(clique_p1.intersection(clique_p2)) == 1 and \
#            clique_p1.intersection(clique_p2) == frozenset({nom_sommet_1}):
#            cpt += 1;
#            dico_C1_C2_S1[(cpt, clique_p1, clique_p2)] = {
#                            "cliques_contractables_1" : clique_p1_p2[0],
#                            "cliques_contractables_2" : clique_p1_p2[1],
#                            "clique_possible_1" : frozenset.union(
#                                                    clique_p1).union(
#                                                    frozenset({nom_sommet_1})
#                                                                ),
#                            "clique_possible_2" : frozenset.union(
#                                                    clique_p2).union(
#                                                    frozenset({nom_sommet_1})
#                                                                )
#                            }
#    print("****** compres => Avant " \
#                 +" nom_sommet_z={}".format(nom_sommet_1) \
#                 +" dico_C1_C2_S1={}".format(len(dico_C1_C2_S1))
#                ) if DBG else None;  
#    logger.debug("****** compres => Avant " \
#                 +" nom_sommet_z={}".format(nom_sommet_1) \
#                 +" dico_C1_C2_S1={}".format(len(dico_C1_C2_S1))
#                ) if DBG else None;   
#            
#    
#    
#    # determination de pi1_pi2_ps
#    nb_prod_cartesien = len(dico_C1_C2_S1);
#    nbre_elts_pi1_pi2 = math.ceil( nb_prod_cartesien * number_items_pi1_pi2);
#    cpt_prod_cartesien = 0;
#    dico_p1_p2_ps = dict();
#    
#    if not dico_C1_C2_S1:
#        ens_cliq_a_supprimer, aretes_ps = set(), set();
#        dico_sommets_corriges, dico_sommets_non_corriges = dict(), dict();
#        cliques_par_nom_sommets_new = dict();
#        
#        cliqs_couv_new, \
#        aretes_LG_k_alpha_new, \
#        dico_sommets_corriges, \
#        dico_sommets_non_corriges, \
#        cliques_par_nom_sommets_new, \
#        sommets_LG_new = \
#            mise_a_jour_aretes_cliques(
#                    nom_sommet_z = nom_sommet_1,
#                    cliques_couvertures_new = set(cliques_couvertures).copy(), 
#                    aretes_LG_k_alpha_new = set(aretes_LG_k_alpha_cor).copy(), 
#                    aretes_ps = aretes_ps,
#                    noms_sommets_1 = noms_sommets_1.copy(),
#                    sommets_LG = sommets_LG,
#                    cliques_par_nom_sommets = cliques_par_nom_sommets.copy())
#        
#        dico_p1_p2_ps[cpt_prod_cartesien] = {
#                    "id_sommet_1": id_sommet_1,
#                    "nom_sommet_1": nom_sommet_1,
#                    "p1": frozenset(),
#                    "p2": frozenset(),
#                    "ps": frozenset(),
#                    "S_z": s_z,
#                    "aretes_ajoutees_p1": frozenset(),
#                    "nbre_aretes_ajoutees_p1": np.inf,
#                    "aretes_ajoutees_p2": frozenset(),
#                    "nbre_aretes_ajoutees_p2": np.inf,
#                    "aretes_supprimees_ps": frozenset(),
#                    "nbre_aretes_supprimees_ps": np.inf,
#                    "aretes_LG_k_alpha_new": aretes_LG_k_alpha_new,
#                    "cliques_couvertures_new": cliqs_couv_new,
#                    "sommets_LG_new": sommets_LG_new,
#                    "sommets_corriges": dico_sommets_corriges,
#                    "sommets_non_corriges": dico_sommets_non_corriges,
#                    "cliques_par_nom_sommets_new": cliques_par_nom_sommets_new,
#                    "cliques_supprimees": ens_cliq_a_supprimer,
#                    "cliques_contractables_1": frozenset(),
#                    "cliques_contractables_2": frozenset()
#                            }
#    else:
#        for k_c1_c2_s1, val_cpt_c1_c2_s1 in dico_C1_C2_S1.items():
#            cpt_prod_cartesien += 1;
#            p1 = None; p2 = None;
#            if cpt_prod_cartesien > nbre_elts_pi1_pi2:
#                break;
#                
#            if val_cpt_c1_c2_s1["cliques_contractables_1"] and \
#                not val_cpt_c1_c2_s1["cliques_contractables_2"]:
#                p1 = val_cpt_c1_c2_s1["clique_possible_1"];
#                p2 = frozenset();
#            elif val_cpt_c1_c2_s1["cliques_contractables_1"] and \
#                val_cpt_c1_c2_s1["cliques_contractables_2"]:
#                p1 = val_cpt_c1_c2_s1["clique_possible_1"];
#                p2 = val_cpt_c1_c2_s1["clique_possible_2"];
#            else :
#                print("IMPOSSIBLE cliques_contr_1 ={}, cliques_contr_2={}".format(
#                      len(val_cpt_c1_c2_s1["cliques_contractables_1"]),
#                      len(val_cpt_c1_c2_s1["cliques_contractables_2"])))
#            
#            if p1 is not None and p2 is not None :
#                gamma_1 = sommets_LG[nom_sommet_1].voisins
#                ps = gamma_1 \
#                     - val_cpt_c1_c2_s1["clique_possible_1"].intersection(gamma_1) \
#                     - val_cpt_c1_c2_s1["clique_possible_2"].intersection(gamma_1);
#                
#                aretes_ps = set( frozenset((nom_sommet_1, sommet_ps)) 
#                                    for sommet_ps in ps
#                                )
#                aretes_p1 = set( map(frozenset, it.combinations(p1,2)) )
#                aretes_ajoutees_p1 = aretes_differente(
#                                        aretes_LG_k_alpha_cor, 
#                                        aretes_p1);
#            
#                aretes_p2 = set( map(frozenset, it.combinations(p2,2)) )
#                aretes_ajoutees_p2 = aretes_differente(
#                                        aretes_LG_k_alpha_cor, 
#                                        aretes_p2);
#                                                
#                aretes_LG_k_alpha_new = set(aretes_LG_k_alpha_cor).union(
#                                                aretes_ajoutees_p1.union(
#                                                    aretes_ajoutees_p2
#                                                    )
#                                                );
#                                                
#                cliques_couvertures_new = set(cliques_couvertures.copy());
#                ens_cliq_a_supprimer = set();                                       
#                for cliq_a_supps in [val_cpt_c1_c2_s1["cliques_contractables_1"],
#                                     val_cpt_c1_c2_s1["cliques_contractables_2"]]:
#                    for cliq_a_supp in cliq_a_supps:
#                        ens_cliq_a_supprimer.add(cliq_a_supp);
#                                           
#                for cliq_couv_new in cliques_couvertures_new :
#                    if cliq_couv_new.issubset(val_cpt_c1_c2_s1["clique_possible_1"]) or \
#                        cliq_couv_new.issubset(val_cpt_c1_c2_s1["clique_possible_2"]) :
#                        ens_cliq_a_supprimer.add(cliq_couv_new);
#               
#                cliques_couvertures_new.difference_update(ens_cliq_a_supprimer);
#            
#                cliques_couvertures_new.add( 
#                                        val_cpt_c1_c2_s1["clique_possible_1"] );
#                cliques_couvertures_new.add( 
#                                        val_cpt_c1_c2_s1["clique_possible_2"] ) \
#                          if val_cpt_c1_c2_s1["clique_possible_2"] else None;
#        
#        
#                dico_sommets_corriges, dico_sommets_non_corriges = dict(), dict();
#                cliques_par_nom_sommets_new = dict();
#                
#                cliqs_couv_new, \
#                aretes_LG_k_alpha_new, \
#                dico_sommets_corriges, \
#                dico_sommets_non_corriges, \
#                cliques_par_nom_sommets_new, \
#                sommets_LG_new = \
#                    mise_a_jour_aretes_cliques(
#                        nom_sommet_z = nom_sommet_1,
#                        cliques_couvertures_new = set(cliques_couvertures_new).copy(), 
#                        aretes_LG_k_alpha_new = set(aretes_LG_k_alpha_new).copy(), 
#                        aretes_ps = aretes_ps,
#                        noms_sommets_1 = noms_sommets_1.copy(),
#                        sommets_LG = sommets_LG,
#                        cliques_par_nom_sommets = cliques_par_nom_sommets.copy())
#                
#                dico_p1_p2_ps[cpt_prod_cartesien] = {
#                        "id_sommet_1": id_sommet_1,
#                        "nom_sommet_1": nom_sommet_1,
#                        "p1": val_cpt_c1_c2_s1["clique_possible_1"],
#                        "p2": val_cpt_c1_c2_s1["clique_possible_2"],
#                        "ps": ps,
#                        "S_z": s_z,
#                        "aretes_ajoutees_p1": aretes_ajoutees_p1,
#                        "nbre_aretes_ajoutees_p1": len(aretes_ajoutees_p1),
#                        "aretes_ajoutees_p2": aretes_ajoutees_p2,
#                        "nbre_aretes_ajoutees_p2": len(aretes_ajoutees_p2),
#                        "aretes_supprimees_ps": aretes_ps,
#                        "nbre_aretes_supprimees_ps": len(aretes_ps),
#                        "aretes_LG_k_alpha_new": aretes_LG_k_alpha_new.copy(),
#                        "cliques_couvertures_new": cliqs_couv_new,
#                        "sommets_LG_new": sommets_LG_new,
#                        "sommets_corriges": dico_sommets_corriges,
#                        "sommets_non_corriges": dico_sommets_non_corriges,
#                        "cliques_par_nom_sommets_new": cliques_par_nom_sommets_new,
#                        "cliques_supprimees" : ens_cliq_a_supprimer,
#                        "cliques_contractables_1": set(val_cpt_c1_c2_s1["cliques_contractables_1"]),
#                        "cliques_contractables_2": set(val_cpt_c1_c2_s1["cliques_contractables_2"])
#                        } 
#                
#            pass # end for
#        pass # end else
#        
#    logger.debug("****** compres ===> Fin compression " \
#                 +" sommet_z : {}, ".format(nom_sommet_1) \
#                 +" nbre_elts_pi1_pi2:{}, ".format(nbre_elts_pi1_pi2) \
#                 +" dico_C1_C2_S1:{}, ".format(len(dico_C1_C2_S1)) \
#                 +" dico_p1_p2_ps:{}".format(len(dico_p1_p2_ps))
#          )  
#    return dico_p1_p2_ps;
#    pass
################################################################################
##               calcul de la compression d un sommet => fin
################################################################################
#
################################################################################
##      critere selection compression: critere local et global => debut
################################################################################
#def critere_C2_C1_local(dico_compression, 
#                        mode_correction,
#                        critere_correction, 
#                        DBG):
#    """ 
#    selectionner le dico selon C2 puis C1 parmi les compressions possibles 
#        sommet a corriger (sommet a -1)
#    
#    C2 : le maximum de sommets corriges
#        * choisir le sommet a -1 qui corrige le max de sommets a -1 possibles
#    C1 : le minimum d'aretes corriges
#    
#    dico_compression : dictionnaire contenant les compressions (p1,p2,ps) du 
#                        sommet sommet_z
#    """
#    max_c2 = 0;
#    min_c1 = np.inf;
#    dico_c1_c2 = dict();
#    critere = "";
#    
#    if not dico_compression :
#        #print("@@CritereLocal: dico_compression={}".format( len(dico_compression) ))
#        return min_c1, max_c2, [];
#    
#    # definition de C2
#    if critere_correction == "voisins_corriges":                               # C2
#        critere = "C2";
#        for cpt_prod_cartesien, dico_p1_p2_ps in dico_compression.items():
#            if len(dico_p1_p2_ps["sommets_corriges"]) >= max_c2:
#                max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
#                nbre_aretes_corriges = \
#                            dico_p1_p2_ps["nbre_aretes_ajoutees_p1"] + \
#                            dico_p1_p2_ps["nbre_aretes_ajoutees_p2"] + \
#                            dico_p1_p2_ps["nbre_aretes_supprimees_ps"];
#                min_c1 = nbre_aretes_corriges if min_c1 >= nbre_aretes_corriges \
#                                                else min_c1;
#                if (min_c1,max_c2) not in dico_c1_c2:
#                    dico_c1_c2[(min_c1,max_c2)] = [dico_p1_p2_ps];
#                else:
#                    dico_c1_c2[(min_c1,max_c2)].append(dico_p1_p2_ps);
#    
#    # definition de C1
#    elif critere_correction == "nombre_aretes_corrigees":                      # C1
#        critere = "C1";
#        for cpt_prod_cartesien, dico_p1_p2_ps in dico_compression.items():
#            nbre_aretes_corriges = \
#                        dico_p1_p2_ps["nbre_aretes_ajoutees_p1"] + \
#                        dico_p1_p2_ps["nbre_aretes_ajoutees_p2"] + \
#                        dico_p1_p2_ps["nbre_aretes_supprimees_ps"];
#            min_c1 = nbre_aretes_corriges if min_c1 >= nbre_aretes_corriges \
#                                          else min_c1;
#            max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
#            if (min_c1,max_c2) not in dico_c1_c2:
#                dico_c1_c2[(min_c1,max_c2)] = [dico_p1_p2_ps];
#            else:
#                dico_c1_c2[(min_c1,max_c2)].append(dico_p1_p2_ps);
#    
#    # definition de C2 puis de C1
#    elif critere_correction == "voisins_nombre_aretes_corrigees":              # C2_C1
#        critere = "C2_C1"
#        for cpt_prod_cartesien, dico_p1_p2_ps in dico_compression.items() :
#            if len(dico_p1_p2_ps["sommets_corriges"]) >= max_c2 :
#                max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
#                nbre_aretes_corriges = \
#                            dico_p1_p2_ps["nbre_aretes_ajoutees_p1"] + \
#                            dico_p1_p2_ps["nbre_aretes_ajoutees_p2"] + \
#                            dico_p1_p2_ps["nbre_aretes_supprimees_ps"];
#                min_c1 = nbre_aretes_corriges if min_c1 >= nbre_aretes_corriges \
#                                                else min_c1;
#                if (min_c1,max_c2) not in dico_c1_c2:
#                    dico_c1_c2[(min_c1,max_c2)] = [dico_p1_p2_ps];
#                else:
#                    dico_c1_c2[(min_c1,max_c2)].append(dico_p1_p2_ps);
#                    
#    logger.debug("@@CritereLocal: critere={},".format(critere) \
#                 +" min_c1={},".format(min_c1) \
#                 +" max_c2={},".format(max_c2) \
#                 +" cles_dico_c1_c2={},".format(set(dico_c1_c2.keys())) \
#                 +" dico_c1_c2={},".format(len(dico_c1_c2[(min_c1,max_c2)])) \
#                 +" dico_compression={}".format(len(dico_compression))) \
#                 if DBG else None
#                          
#    if not dico_c1_c2:
#        return min_c1, max_c2, [];
#    else:
#        return min_c1, max_c2, dico_c1_c2[(min_c1,max_c2)];
#    pass
#
#def rechercher_min_max(liste_tuples, critere):
#    """ retourne la tuple (min, max)
#    """
#    min_c1 = np.inf;
#    if len(liste_tuples) == 0:
#        return np.inf
#    
#    #max_c2 = 0;
#    if critere == "C1":
#        return min(liste_tuples)
#    elif critere == "C2":
#        return max(liste_tuples)
#    elif critere == "C2_C1":
#        liste_intermediaires = [];
#        min_c1 = min(liste_tuples)[0]
#        for tuple_ in liste_tuples:
#           if tuple_[0] == min_c1:
#               liste_intermediaires.append(tuple_)
#        return max(liste_intermediaires)
#
#def critere_C2_C1_global(dico_compression,
#                         mode_correction,
#                         critere_correction,
#                         DBG):
#    """ recherche la compression optimale parmi tous les sommets a corriger 
#        selon les criteres C1 et C2.
#        
#    C2 : le maximum de sommets corriges
#        * choisir le sommet a -1 qui corrige le max de sommets a -1 possibles
#    C1 : le minimum d'aretes corriges 
#    
#    methode de selection C2
#        je cherche le min local de c1 pour tous les sommets a corriger
#        parmi les min locaux, je cherche le max global de c2
#        une fois la liste des (min_global,max_global), je prends le 1er element.
#    """
#    
#    max_c2_global = 0;
#    min_c1_global = np.inf;
#    dico_c1_c2_global = dict();
#    cle_min_max_c2 = None;
#    
#    critere = ""
#    
#    if critere_correction == "voisins_corriges":                               # C2
#        critere = "C2"
#        for id_sommet_z, dicos_p1_p2_ps in dico_compression.items():
#            # selection de dico selon C1
#            min_c1_local = dicos_p1_p2_ps[0];
#            max_c2_local = dicos_p1_p2_ps[1];
#            
#            for dico_p1_p2_ps in dicos_p1_p2_ps[2]:
#                nbre_aretes_corriges = \
#                                dico_p1_p2_ps["nbre_aretes_ajoutees_p1"] \
#                                + dico_p1_p2_ps["nbre_aretes_ajoutees_p2"] \
#                                + dico_p1_p2_ps["nbre_aretes_supprimees_ps"];
#                min_c1_local = nbre_aretes_corriges \
#                                if min_c1_local >= nbre_aretes_corriges \
#                                else min_c1_local;
#                if (min_c1_local,max_c2_local) not in dico_c1_c2_global :
#                    dico_c1_c2_global[(min_c1_local, 
#                                       max_c2_local)] = [dico_p1_p2_ps];
#                else:
#                    dico_c1_c2_global[(min_c1_local, 
#                                       max_c2_local)].append(dico_p1_p2_ps);
#                           
#        # selection selon C2
#        cle_min_max_c2 = rechercher_min_max(dico_c1_c2_global.keys(), "C2");
#    
#    elif critere_correction == "nombre_aretes_corrigees":                      # C1
#        critere = "C1"
#        for id_sommet_z, dicos_p1_p2_ps in dico_compression.items() :
#            # selection de dico selon C2
#            max_c2_local = dicos_p1_p2_ps[1];
#            min_c1_local = dicos_p1_p2_ps[0];
#            
#            for dico_p1_p2_ps in dicos_p1_p2_ps[2] :
#                nbre_sommets_corriges = len(dico_p1_p2_ps["sommets_corriges"]);
#                max_c2_local = nbre_sommets_corriges \
#                                if nbre_sommets_corriges > max_c2_local \
#                                else max_c2_local;
#                if (min_c1_local,max_c2_local) not in dico_c1_c2_global:
#                    dico_c1_c2_global[(min_c1_local, 
#                                       max_c2_local)] = [dico_p1_p2_ps];
#                else:
#                    dico_c1_c2_global[(min_c1_local, 
#                                       max_c2_local)].append(dico_p1_p2_ps);
#
#        # selection selon C1
#        cle_min_max_c2 = rechercher_min_max(dico_c1_c2_global.keys(), "C1");
#        
#    elif critere_correction == "voisins_nombre_aretes_corrigees":              # C2_C1
#        critere = "C2_C1"
#        for id_sommet_z, dicos_p1_p2_ps in dico_compression.items():
#            min_c1_local = dicos_p1_p2_ps[0]; #np.inf
#            max_c2_local = dicos_p1_p2_ps[1]; #0
#            
#            for dico_p1_p2_ps in dicos_p1_p2_ps[2] :
#                nbre_sommets_corriges = len(dico_p1_p2_ps["sommets_corriges"]);
#                max_c2_local = nbre_sommets_corriges \
#                                if nbre_sommets_corriges > max_c2_local \
#                                else max_c2_local;
#                nbre_aretes_corriges = \
#                                dico_p1_p2_ps["nbre_aretes_ajoutees_p1"] \
#                                + dico_p1_p2_ps["nbre_aretes_ajoutees_p2"] \
#                                + dico_p1_p2_ps["nbre_aretes_supprimees_ps"];
#                min_c1_local = nbre_aretes_corriges \
#                                if min_c1_local >= nbre_aretes_corriges \
#                                else min_c1_local;
#                if (min_c1_local,max_c2_local) not in dico_c1_c2_global:
#                    dico_c1_c2_global[(min_c1_local, 
#                                       max_c2_local)] = [dico_p1_p2_ps];
#                else:
#                    dico_c1_c2_global[(min_c1_local, 
#                                       max_c2_local)].append(dico_p1_p2_ps);
#                
#        # selection selon C2_C1
#        cle_min_max_c2 = rechercher_min_max(dico_c1_c2_global.keys(), "C2_C1");
#    
#    numero_sol_c1_c2 = np.random.randint(
#                        low=0, 
#                        high=len(dico_c1_c2_global[cle_min_max_c2])
#                        )
#    
#    min_c1_global = cle_min_max_c2[0];
#    max_c2_global = cle_min_max_c2[1];
#    
#    logger.debug("@@CritereGlobal critere={}, ".format(critere) \
#                 +" dico_compression={}, ".format(len(dico_compression)) \
#                 +" cle_globale={}, ".format(set(dico_c1_c2_global.keys())) \
#                 +" cle_min_max_c2={}, ".format(cle_min_max_c2) \
#                 +" nbre_sol_c1_c2={}".format(
#                         len(dico_c1_c2_global[cle_min_max_c2]))
#                 ) if DBG else None;
#                 
#    return min_c1_global, \
#            max_c2_global, \
#            dico_c1_c2_global[cle_min_max_c2][numero_sol_c1_c2];
################################################################################
##      critere selection compression: critere local et global => fin
################################################################################
#
################################################################################
##                application de la correction => debut
################################################################################
#def appliquer_correction(dico_sol_C2_C1,
#                         sommets_1,
#                         DBG):
#    """ appliquer la compression choisie dans le graphe.
#    """
#    """
#    {
#    "id_sommet_1": id_sommet_1,
#    "sommet_1": nom_sommet_1,
#    "p1": val_cpt_c1_c2_s1["clique_possible_1"],
#    "p2": val_cpt_c1_c2_s1["clique_possible_2"],
#    "ps": ps,
#    "S_z": s_z,
#    "aretes_ajoutees_p1": aretes_ajoutees_p1,
#    "nbre_aretes_ajoutees_p1": len(aretes_ajoutees_p1),
#    "aretes_ajoutees_p2": aretes_ajoutees_p2,
#    "nbre_aretes_ajoutees_p2": len(aretes_ajoutees_p2),
#    "aretes_supprimees_ps": aretes_ps,
#    "nbre_aretes_supprimees_ps": len(aretes_ps),
#    "aretes_LG_k_alpha_new": aretes_LG_k_alpha_new.copy(),
#    "cliques_couvertures_new": cliqs_couv_new,
#    "sommets_LG_new": sommets_LG_new,
#    "sommets_corriges": dico_sommets_corriges,
#    "sommets_non_corriges": dico_sommets_non_corriges,
#    "cliques_par_nom_sommets_new": cliques_par_nom_sommets_new,
#    "cliques_supprimees" : ens_cliq_a_supprimer,
#    "cliques_contractables_1" : set(val_cpt_c1_c2_s1["cliques_contractables_1"]),
#    "cliques_contractables_2" : set(val_cpt_c1_c2_s1["cliques_contractables_2"])
#    } 
#    """
#    
#    cliques_couvertures = set();
#    cliques_couvertures = dico_sol_C2_C1['cliques_couvertures_new'];
#    aretes_LG_k_alpha = dico_sol_C2_C1['aretes_LG_k_alpha_new'];
#    sommets_LG = dico_sol_C2_C1['sommets_LG_new'];
#    cliques_par_nom_sommets = dico_sol_C2_C1['cliques_par_nom_sommets_new'];
#    
#    id_sommets_1 = set(dico_sol_C2_C1["sommets_corriges"].keys());
#    id_sommets_1.add(dico_sol_C2_C1["id_sommet_1"]);
#    sommets_corriges = dico_sol_C2_C1["sommets_corriges"].values();
#    
#    logger.debug(
#            "*** Avant correction : id_sommets_1:{}, ".format(id_sommets_1) \
#            +" sommets_corriges={}, ".format(sommets_corriges) \
#            +" sommet_1={}".format(dico_sol_C2_C1["nom_sommet_1"])) \
#        if DBG else None;
#                
#    sommets_1 = np.delete(sommets_1, list(id_sommets_1)).tolist();
#    logger.debug("*** Apres correction : "
#                  +"sommets_1 restants = {}".format(sommets_1)) \
#        if DBG else None;
#                     
#    if set(sommets_1).intersection(set(sommets_corriges)) :
#        print("---ERROR : sommets {} suppression : NOK -----".
#              format(sommets_corriges))
#                     
#    """
#    cliques_couv_new,\
#            aretes_LG_k_alpha_cor_new,\
#            sommets_LG,\
#            noms_sommets_1
#    """
#    return cliques_couvertures,\
#            aretes_LG_k_alpha,\
#            sommets_LG,\
#            cliques_par_nom_sommets,\
#            sommets_1;
################################################################################
##               application de la correction => fin
################################################################################
#
################################################################################
##                correction aleatoire des sommets_1 => debut
################################################################################
#def correction_aleatoire(sommets_1, 
#                         aretes_LG_k_alpha_cor,
#                         dico_correction,
#                         dico_parametres_new):
#    """
#    correction des sommets a corriger par une selection aleatoire de ceux-ci.
#    """
#    cliques_couvertures_cor = dico_correction["C"].copy();
#    cliques_par_sommets = dico_correction["sommets_par_cliqs_avec_aretes"].copy();
#    
#    mode_select_noeuds_1 = dico_parametres_new["mode_select_noeuds_1s"][0];
#    critere_select_pi1_pi2 = dico_parametres_new["critere_selection_pi1_pi2"];
#    number_items_pi1_pi2 = dico_parametres_new["number_items_pi1_pi2"];
#    DBG = dico_parametres_new["DBG"]
#    
#    dico_sommets_corriges = dict();
#    
#    id_sommet_1 = 0;
#    while len(sommets_1)>0:
#        sommet_1 = sommets_1.pop(id_sommet_1);
#        
#        dicos_p1_p2_ps = dict();
#        dicos_p1_p2_ps = compression_sommet(id_sommet_1,
#                            sommet_1 = sommet_1,
#                            sommets_1 = sommets_1,
#                            cliques_par_sommets = cliques_par_sommets,
#                            cliques_couv_cor = cliques_couvertures_cor,
#                            aretes_LG_k_alpha_cor = aretes_LG_k_alpha_cor,
#                            mode_select_noeuds_1 = mode_select_noeuds_1,
#                            critere_select_pi1_pi2 = critere_select_pi1_pi2,
#                            number_items_pi1_pi2 = number_items_pi1_pi2,
#                            DBG = DBG)
#                                           
#        
#        dico_sol_p1_p2_ps, dico_compression = dict(), dict(); 
#        min_c1 = 0; max_c2 = 0;
#        min_c1, max_c2, dico_compression = critere_C2_C1_local(
#                                            dicos_p1_p2_ps,
#                                            mode_select_noeuds_1,
#                                            critere_select_pi1_pi2,
#                                            DBG = DBG)
#        
#        print("\n\n len(dico_compression)={}".format(len(dico_compression)))
#        nbre_alea = random.randint(
#                        0,
#                        len(dico_compression)-1)
#        dico_sol_p1_p2_ps = dico_compression[nbre_alea];
#        
#        if not dico_sol_p1_p2_ps:
#            cout_T = {"aretes_ajoutees_p1": frozenset(),
#                      "aretes_ajoutees_p2": frozenset(),
#                      "aretes_supprimees": frozenset(),
#                      "min_c1": min_c1,
#                      "max_c2": max_c2};
#            id_sommet_1 += 1;
#            dico_sommets_corriges[("0_0", "0_0")] = {
#                        "compression_p1": frozenset(),
#                        "compression_p2": frozenset(),
#                        "compression_ps": frozenset(),
#                        "sommets_corriges": dict(), # voisins_corriges = {"id_voisin_ds_sommets_a_corriger":voisin}
#                        "cout_T": cout_T
#                        }
#            sommets_1 = list();
#        else:
#            cliques_couv_new, \
#            aretes_LG_k_alpha_cor_new, \
#            sommets_LG_new, \
#            cliques_par_nom_sommets_new, \
#            sommets_1 = appliquer_correction(
#                                dico_sol_p1_p2_ps,
#                                sommets_1, 
#                                DBG);
#            
#            # mise a jour variables
#            cliques_par_sommets = cliques_par_nom_sommets_new.copy()
#            cliques_couvertures_cor = cliques_couv_new.copy()
#            aretes_LG_k_alpha_cor = aretes_LG_k_alpha_cor_new.copy();
#            sommets_LG = sommets_LG_new.copy();
#            cout_T = {
#                "aretes_ajoutees_p1":dico_sol_p1_p2_ps["aretes_ajoutees_p1"],
#                "aretes_ajoutees_p2":dico_sol_p1_p2_ps["aretes_ajoutees_p2"],
#                "aretes_supprimees":dico_sol_p1_p2_ps["aretes_supprimees_ps"],
#                "min_c1":min_c1, "max_c2":max_c2};
#            id_sommet_1 += 1;
#            dico_sommets_corriges[(id_sommet_1, 
#                                   dico_sol_p1_p2_ps["sommet_1"])] = {
#                        "compression_p1":dico_sol_p1_p2_ps["p1"],
#                        "compression_p2":dico_sol_p1_p2_ps["p2"],
#                        "compression_ps":dico_sol_p1_p2_ps["ps"],
#                        "sommets_corriges":dico_sol_p1_p2_ps["sommets_corriges"], # voisins_corriges = {"id_voisin_ds_sommets_a_corriger":voisin}
#                        "cout_T": cout_T
#                        }
#        
#    # TODO A REVOIR JE DOIS RETOURNER dico_correction
#    dico_correction["C"] = cliques_couvertures_cor;
#    dico_correction["sommets_par_cliqs_avec_aretes"] = cliques_par_sommets;
#    dico_correction["0"] = 0;
#    return dico_correction, \
#            dico_sommets_corriges;
##            sommets_LG, \
#            
#    pass
################################################################################
##               correction aleatoire des sommets_1 => fin
################################################################################
#
################################################################################
##               algorithme de correction pour k entier naturel => debut
################################################################################
#def correction_cliques_k(dico_correction, 
#                         aretes_matE_k_alpha,
#                         dico_gamma_noeud, 
#                         dico_parametres_new):
#    """
#    corriger les cliques selon le mode de selection des sommets a -1
#    """
#    aretes_LG_k_alpha_cor = aretes_matE_k_alpha.copy();
#    
#    sommets_a_corriger = dico_correction["sommets_a_corriger"];
#    dico_sommets_corriges = dict();
#    
#    if dico_correction["critere_selection_pi1_pi2"] == "voisins_corriges" \
#        and dico_correction["mode_select_noeuds_1s"] == "aleatoire":
#        dico_correction, dico_sommets_corriges = correction_aleatoire(
#                                                    sommets_a_corriger, 
#                                                    aretes_LG_k_alpha_cor,
#                                                    dico_correction,
#                                                    dico_parametres_new);
#        
#    elif dico_correction["critere_selection_pi1_pi2"] == "nombre_aretes_corrigees" \
#        and dico_correction["mode_select_noeuds_1s"] == "aleatoire":
#        dico_correction, dico_sommets_corriges = correction_aleatoire(
#                                                    sommets_a_corriger, 
#                                                    aretes_LG_k_alpha_cor,
#                                                    dico_correction,
#                                                    dico_parametres_new);
#        
#    elif dico_correction["critere_selection_pi1_pi2"] == "voisins_nombre_aretes_corrigees" \
#        and dico_correction["mode_select_noeuds_1s"] == "aleatoire":
#        dico_correction, dico_sommets_corriges = correction_aleatoire(
#                                                    sommets_a_corriger, 
#                                                    aretes_LG_k_alpha_cor,
#                                                    dico_correction,
#                                                    dico_parametres_new);
#                
#    return dico_correction, dico_sommets_corriges;
################################################################################
##               algorithme de correction pour k entier naturel => fin
################################################################################