#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 15:32:30 2019

@author: willy
"""
import math;
import random;
import logging;

import numpy as np;
import itertools as it;

import fonctions_auxiliaires as fct_aux;

###############################################################################
#                          cliques_contractables => debut
###############################################################################
def cliques_contractables(sommet_a_corr,
                            clique_sommet_a_corrs,
                            cliques_couvertures,
                            aretes_Ec,
                            aretes_cliques):
  
    """
    retourne la liste de cliques contractables
    Chaque clique contractable est un couple d'elements (x1, x2) tel que
        x1 = {c1,c2} : clique formee par l'union de C1 et C2
        x2 = [c1,c2] : la liste des cliques contractables
    exple [({c1,c2},[c1, c2]), ...]
     
    """    
    cliques_contractables = list();
    aretes_cliques = list( map(set, aretes_cliques) );
    aretes_Ec = list( map(set, aretes_Ec) );
    for tuple_cliq in it.combinations(clique_sommet_a_corrs, 2):
        booleen = True;
        
        f = lambda arete : set(arete) in aretes_Ec \
                            and set(arete) in aretes_cliques \
                            and set(arete) not in aretes_cliques;
        for arete in it.product(tuple_cliq[0]-{sommet_a_corr}, 
                                tuple_cliq[1]-{sommet_a_corr}):
            booleen = False if list(filter(f, arete)) else True;
            if booleen:
                cliques_contractables.append(
                        (frozenset.union(*tuple_cliq), list(tuple_cliq))
                        );
    
    #ajout des cliques contractables dont unee clique est l'ensemble vide
    for tuple_cliq in it.combinations(clique_sommet_a_corrs, 1):
        cliques_contractables.append( (tuple_cliq[0], [tuple_cliq[0]]) );
        
    return cliques_contractables;
#def cliques_contractables(noeud_z, C_z, ens_C_i, liste_aretes_E_C_i, liste_aretes_ens_C_i):
#    l_cliqs_contractables = list()
##    #print("C_z = ",C_z)
#    #for tuple_cliq in it.combinations(ens_C_i, 2): =====> A Effacer
#    for tuple_cliq in it.combinations(C_z, 2):
#        #tu0 = tuple_cliq[0]; tu1 = tuple_cliq[1];
#        booleen = True
#        for arete in it.product( tuple_cliq[0] - {noeud_z} , tuple_cliq[1] - {noeud_z} ):
#            if ( (arete[0],arete[1]) in liste_aretes_E_C_i or (arete[1],arete[0]) in liste_aretes_E_C_i) \
#                and ( (arete[0],arete[1]) in liste_aretes_ens_C_i or (arete[1],arete[0]) in liste_aretes_ens_C_i )\
#                and set(arete) not in ens_C_i:
#                booleen = False; #print("**** not contractable: ", tuple_cliq, " arete_problem: ", arete)
#                break;
#        if booleen == True:
#            l_cliqs_contractables.append( (tuple_cliq[0].union(tuple_cliq[1]), [tuple_cliq[0],tuple_cliq[1]] ))
#    
#    #ajout des cliques contractables dont unee clique est l'ensemble vide
#    for tuple_cliq in it.combinations(C_z, 1):
#        l_cliqs_contractables.append( (tuple_cliq[0], [tuple_cliq[0]]) )
#    return l_cliqs_contractables 
###############################################################################
#                          cliques_contractables => fin
###############################################################################

###############################################################################
#                          S_z => debut
###############################################################################
def S_sommet(sommet_a_corr,
               voisins_sommet_a_corr,
               cliques_couvertures,
               aretes_cliques):
    """
    gamma_z tel que 
        - {noeud_z, voisin_z} in cliques_couvertures ou
        - (noeud_z, voisin_z) not in aretes_cliques and 
            (voisin_z, noeud_z) not in aretes_cliques
    return la liste des voisin_z
    """
    s_z = list();
    for voisin_som_a_corr in voisins_sommet_a_corr:
        if {sommet_a_corr, voisin_som_a_corr} in cliques_couvertures \
            or {sommet_a_corr, voisin_som_a_corr} not in aretes_cliques:
            s_z.append(voisin_som_a_corr);
    return s_z;

#def S_z(noeud_z, gamma_z, ens_C_i, liste_aretes_ens_C_i):
#    """
#    gamma_z tel que 
#        - {noeud_z, voisin_z} in ens_C_i ou
#        - (noeud_z, voisin_z) not in liste_aretes_ens_C_i and 
#            (voisin_z, noeud_z) not in liste_aretes_ens_C_i
#    return la liste des voisin_z
#    """
#    s_z = list()
#    for voisin_z in gamma_z:
#        if {noeud_z, voisin_z} in ens_C_i or \
#            ((noeud_z, voisin_z) not in liste_aretes_ens_C_i and \
#             (voisin_z, noeud_z) not in liste_aretes_ens_C_i):
#            s_z.append(voisin_z)
#    return s_z;
###############################################################################
#                          S_z => fin
###############################################################################

###############################################################################
#         voisine_sommet et D_sommet_c d'un sommet a corriger => debut
###############################################################################
def voisine_sommet(sommet_a_corr, 
                   clique,
                   clique_sommet_a_corrs,
                   s_sommet_a_corr, 
                   n = 2):
    """
    clique "clique" n'appartenant pas a clique_sommet_a_corrs
        et card("clique" \cup "s_sommet_a_corr") > n
    return un booleen tel que 
        si False alors la clique "clique" n'est pas voisine de sommet_a_corr
        si True alors la clique "clique" est voisine de sommet_a_corr
    """
    if clique in clique_sommet_a_corrs:
        return False;
    if len( clique.intersection(s_sommet_a_corr) ) > n-2:
        return True;

def D_sommet_c(clique, 
               voisins_sommet_a_corr, 
               clique_sommet_a_corrs):
    """
    determine si une clique clique_sommet_a_corr contenant "sommet_a_corr_z" 
            est dependance avec "clique"
            si oui retourne les listes des cliques dependante avec "clique"
    return
        list_clique: liste de cliques dependantes de "clique" 
        
    clique_sommet_a_corrs : liste des cliques contenant noeud_z = sommet_a_corr 
    voisins_sommet_a_corr : un ens de voisins de noeud_z = sommet_a_corr 
                            cad il est un set
    
    NB : "clique" est une clique voisine de noeud_z
    """
    list_clique = list()
    list_clique = [cliq_z for cliq_z in clique_sommet_a_corrs 
                   if len(clique.intersection(
                           cliq_z.intersection(voisins_sommet_a_corr))) != 0]
    return list_clique;
    
#def voisine_z(noeud_z, C, C_z, S_z, n = 2):
#    """
#    clique C n'appartenant pas a C_z et card(C \cup S_Z) > n
#    return un booleen tel que 
#        si False alors la clique C n'est pas voisine de noeud_z
#        si True alors la clique C est voisine de noeud_z
#    """
#    if C in C_z:
#        return False
#    #if len( C.intersection(S_z) ) > n-1:   # == >= 1
#    if len( C.intersection(S_z) ) > n-2:   # == >= 1
#        return True
#
#def D_z_c( C, gamma_z, C_z):
#    """
#    determine si une clique c_z contenant "noeud_z" est dependance avec C
#            si oui retourne les listes des cliques dependante avec C
#    return
#        list_clique: liste de cliques dependantes de C 
#        
#    C_z: liste des cliques contenant noeud_z = z 
#    gamma_z: un ens de voisins de noeud_z = z cad il est un set
#    
#    NB : C est une clique voisine de noeud_z
#    """
#    list_clique = list()
#    list_clique = [cliq_z for cliq_z in C_z if len(C.intersection(cliq_z.intersection(gamma_z))) != 0]
#    return list_clique
###############################################################################
#           voisine_sommet et D_sommet_c d'un sommet a corriger => fin
###############################################################################
    
###############################################################################
#                    liste de cliques contractables => debut
###############################################################################
def liste_cliques_contractables(clique, 
                                liste_cliques,
                                sommet_a_corr, 
                                clique_sommet_a_corrs, 
                                cliques_couvertures, 
                                aretes_cliques):
    """
    verifie si une paire de cliques de liste_cliques est contractable.
    Si cette paire est contractable alors
        * on ajoute [ paire[0], paire[1], C] a l'ensemble des listes augmentantes
    """
    
    f = lambda arete: set(arete) in aretes_cliques \
                        and set(arete) not in cliques_couvertures;
    
    if len(liste_cliques) == 0:
        return [clique, set()];
    elif len(liste_cliques) == 1:
        return [liste_cliques[0], clique];
    else:
        res_liste_cliques = list();
        for tu_contr_possib in it.combinations(liste_cliques, 2):
            set_inter = set.intersection(*tu_contr_possib);
            set_inter = set_inter - {sommet_a_corr};
            tu0 = set(tu_contr_possib[0]) - set_inter;
            tu1 = set(tu_contr_possib[1]) - set_inter;
            
            booleen = True;
            booleen = False if list(filter(f, it.product(tu0, tu1))) else True;
            if booleen:
                res_liste_cliques.append(
                        (tu_contr_possib[0], tu_contr_possib[1], clique)
                        );
    return res_liste_cliques;
    
#def ens_cliques_contractables( C, liste_cliques, noeud_z, C_z, ens_C_i, liste_aretes_ens_C_i):
#    
#    if len(liste_cliques) == 0:
#        return [(C, set())]
#    if len(liste_cliques) == 1:
#        return [( liste_cliques.pop(), C)]
#    else:
#        res_liste_cliques = list();
#        for tu_contr_possib in it.combinations(liste_cliques, 2):
#            set_inter = set(tu_contr_possib[0]).intersection(set(tu_contr_possib[1]));
#            set_inter = set_inter - {noeud_z};
#            tu0 = set(tu_contr_possib[0]) - set_inter ; tu1 = set(tu_contr_possib[1]) - set_inter;
#        
#            booleen = True;
#            for arete in it.product( tu0 , tu1 ):
#                if ( (arete[0],arete[1]) in liste_aretes_ens_C_i \
#                    or (arete[1],arete[0]) in liste_aretes_ens_C_i) \
#                    and set(arete) not in ens_C_i:
#                        booleen = False; break;
#            if booleen == True:
#                res_liste_cliques.append( (tu_contr_possib[0], tu_contr_possib[1], C) )
#                
#        return res_liste_cliques   
###############################################################################
#                    liste de cliques contractables => fin
###############################################################################
    
###############################################################################
#                    augmentation d'un sommet a corriger => debut
###############################################################################
def augmentation_sommet(sommet_a_corr, 
                        voisins_sommet_a_corr, 
                        clique_sommet_a_corrs, 
                        s_sommet_a_corr, 
                        cliques_couvertures, 
                        aretes_cliques, 
                        n=2):
    """
    determine l augmentation d'un sommet a corriger et retourne les cliques 
    selon le format 
        augm_sommet_a_corrs: [ ({1,2,3,sommet_a_corr},[{1,2,3}]), ... ]
    """
    augm_sommet_a_corrs = list()
    for clique in cliques_couvertures:
        if voisine_sommet(sommet_a_corr, 
                          clique,
                          clique_sommet_a_corrs,
                          s_sommet_a_corr, 
                          n):
            liste_cliques = D_sommet_c(clique, 
                                       voisins_sommet_a_corr, 
                                       clique_sommet_a_corrs);
            liste_cliques = liste_cliques_contractables( 
                                clique, 
                                liste_cliques,
                                sommet_a_corr, 
                                clique_sommet_a_corrs, 
                                cliques_couvertures, 
                                aretes_cliques);
            while len(liste_cliques) != 0:
                cliq_contratables = liste_cliques.pop();
                augm_sommet_a_corrs.append( 
                        (set.union( *cliq_contratables).union({sommet_a_corr}), 
                         list(cliq_contratables)) 
                        );
                        
    return augm_sommet_a_corrs;
#def augmentation_z(noeud_z, gamma_z, C_z, S_z, ens_C_i, liste_aretes_ens_C_i, n=2):
#    """
#    format l_augm_z: [ ({1,2,3,z},[{1,2,3}]), ... ]
#    """
#    l_augm_z = list()
#    for C in ens_C_i:
#        liste_cliques = list()
#        if voisine_z(noeud_z, C, C_z, S_z, n) == True:
#            liste_cliques = D_z_c( C, gamma_z, C_z)
#            liste_cliques = ens_cliques_contractables( C, liste_cliques, \
#                                                       noeud_z, C_z, ens_C_i, 
#                                                       liste_aretes_ens_C_i)
#            while len(liste_cliques) != 0:
#                cliq_contratables = liste_cliques.pop();
#                l_augm_z.append( (set().union( *cliq_contratables).union({noeud_z}), 
#                                  list(cliq_contratables)) )
#    return l_augm_z
###############################################################################
#                    augmentation d'un sommet a corriger => fin
###############################################################################

###############################################################################
#                       determiner C1_S1  => debut
###############################################################################
def all_subsets(ss):
    #print("all ss ", ss)
    return it.chain(*map(lambda x: it.combinations(ss, x), range(1, len(ss)+1)))

def combinaison_possible_C1_S1(C1, s1):
    """
    p1 = c1[0]; p2= s1.pop(), p3 = s1.pop(), .....
    resultat p1p2p3, p1p2, p1p3, p2p3
    
    NB: TOUTES LES COMBINAISONS POSSIBLE DE GROUPE DE 1, 2 ,... LEN(N) 
    AVEC N LA TAILLE DE LA LISTE DE [C1, s1] 
    """
#    #print("C1 = ", C1)
#    #print("s1 = ", s1)
    l = list(); l_r = list()
    l.append(C1[0]); l.extend(s1);
    for subset in all_subsets(l):
        if C1[0] in subset:
            l_r.append( (set().union(*subset),C1[1]) )
    return l_r
    
def C1_S1(sommet_a_corr, C1, s_sommet_a_corr, cliques_couvertures):
    """
    sommets v de s_sommet_a_corr (S_z) 
    n'appartenant a aucune clique de C1 tel que:
        v in S_z, x in C1, !C' in ens_C_i dont len(C') >2 et {v,x} convert par C' ===> "1"
        
    C1: clique contractable de la forme ({C1,C2}, [C1,C2])
    ens_C_i : cliques_couvertures
    return liste de sous ensemble verifiant "1"
    
    exple: 
        C1 = ({1,2,3},[{1,2,3}]);s1= [{4},{5}];
        ===> return [ ({1, 2, 3}, [{1, 2, 3}]), ({1, 2, 3, 4}, [{1, 2, 3}]),
         ({1, 2, 3, 5}, [{1, 2, 3}]), ({1, 2, 3, 4, 5}, [{1, 2, 3}])]
    """
    l_sous_ensemble_C1_S1 = list();
    if len(C1) == 0:
        for i in range(1,len(s_sommet_a_corr)+1):
            for tu in it.combinations(s_sommet_a_corr, i):
                l_sous_ensemble_C1_S1.append( 
                        (set().union( *[{sommet_a_corr}, tu] ), []) 
                        );                                                       #[{'z'},{1,2,3}]
    else:
        S1_tmp = list(); s1 = list();
        S1_tmp = [s for s in s_sommet_a_corr if s not in C1[0]] 
        ens_C_i_sup2 = [Cliq for Cliq in cliques_couvertures if len(Cliq) > 2]
#        l_aretes_cliques = aretes_ens_C_i(ens_C_i_sup2);
        l_aretes_cliques = fct_aux.determiner_aretes_cliques(ens_C_i_sup2)
        
        while len(S1_tmp) != 0:
            booleen = True; noeud_S1_tmp = S1_tmp.pop()
            for arete in it.product({noeud_S1_tmp},C1[0]):
                if (arete[0], arete[1]) in l_aretes_cliques \
                    or (arete[1], arete[0]) in l_aretes_cliques:
                    booleen = False; 
                    break;
            if booleen == True:
                s1.append({noeud_S1_tmp})
        
        # ici faire combinaison de C1[0] et S1 A DEPLACER
        l = combinaison_possible_C1_S1(C1, s1);
        l_sous_ensemble_C1_S1.extend(l)
    return l_sous_ensemble_C1_S1;

#def C1_S1(noeud_z, C1, S_z, ens_C_i):
#    """
#    sommets v de S_z n'appartenant a aucune clique de C1 tel que:
#        v in S_z, x in C1, !C' in ens_C_i dont len(C') >2 et {v,x} convert par C' ===> 1
#        
#    C1: clique contractable de la forme ({C1,C2}, [C1,C2])
#    return liste de sous ensemble verifiant 1
#    
#    exple: 
#        C1 = ({1,2,3},[{1,2,3}]);s1= [{4},{5}];
#        ===> return [ ({1, 2, 3}, [{1, 2, 3}]), ({1, 2, 3, 4}, [{1, 2, 3}]),
#         ({1, 2, 3, 5}, [{1, 2, 3}]), ({1, 2, 3, 4, 5}, [{1, 2, 3}])]
#    """
#    #TOD mettre une entree dans fichier log "run fonction S1 ===> debut"
#    l_sous_ensemble_C1_S1 = list()
#    if len(C1) == 0:
#        for i in range(1,len(S_z)+1):
#            for tu in it.combinations(S_z, i):
#                l_sous_ensemble_C1_S1.append( (set().union( *[{noeud_z},tu] ),[]) )   #[{'z'},{1,2,3}]
#    else:
#        S1_tmp = list(); s1 = list()
#        S1_tmp = [s for s in S_z if s not in C1[0]] 
#        ens_C_i_sup2 = [Cliq for Cliq in ens_C_i if len(Cliq) > 2]
#        l_aretes_cliques = aretes_ens_C_i(ens_C_i_sup2)
#        
#        while len(S1_tmp) != 0:
#            booleen = True; noeud_S1_tmp = S1_tmp.pop()
#            for arete in it.product({noeud_S1_tmp},C1[0]):
#                if (arete[0], arete[1]) in l_aretes_cliques or (arete[1], arete[0]) in l_aretes_cliques:
#                    booleen = False; 
#                    # TOD mettre une entree dans fichier log
#                    break;
#            if booleen == True:
#                s1.append({noeud_S1_tmp})
#                # TOD ici une entre dans le log
#        
#        # ici faire combinaison de C1[0] et S1 A DEPLACER
#        l = combinaison_possible_C1_S1(C1, s1);
#        l_sous_ensemble_C1_S1.extend(l)
#    #TOD mettre une entree dans fichier log "run fonction S1 ===> fin"
#    return l_sous_ensemble_C1_S1
###############################################################################
#                       determiner C1_S1  => fin
###############################################################################
    
###############################################################################
#                       determiner C1 et S1  => debut
###############################################################################
def determiner_C1_S1(sommet_a_corr, 
                     clique_contractables, 
                     s_sommet_a_corr, 
                     cliques_couvertures):
    """
    return tous les C1 \cup S1
    format return 
    [({'E', 'H', 'I'}, [{'E', 'H', 'I'}]), 
     ({'F', 'E'}, []), 
     ({'D', 'E', 'G'}, [{'D', 'E', 'G'}]) ]
    """
    C1_S1_s = list()
    if len(clique_contractables) == 0:
        C1_S1_s.extend( C1_S1(sommet_a_corr, (), 
                              s_sommet_a_corr, cliques_couvertures) )
    for C1 in clique_contractables:
        #print("cliq contractables = ", C1)
        sub_set_S1 = C1_S1(sommet_a_corr, C1, 
                           s_sommet_a_corr, cliques_couvertures);
        C1_S1_s.extend(sub_set_S1)
    
    # a ajouter toutes combinaison possibles dans l_C1_S1
    return C1_S1_s;

#def set_C1_S1(noeud_z, cliqs_contract, S_z, ens_C_i):
#    l_C1_S1 = list()
#    if len(cliqs_contract) == 0:
#        l_C1_S1.extend( C1_S1(noeud_z, (), S_z, ens_C_i) )
#    for C1 in cliqs_contract:
#        #print("cliq contractables = ", C1)
#        sub_set_S1 = C1_S1(noeud_z, C1, S_z, ens_C_i);
#        l_C1_S1.extend(sub_set_S1)
##    #print("l_C1_S1 = ", l_C1_S1)
#    
#    # a ajouter toutes combinaison possibles dans l_C1_S1
#    return l_C1_S1
###############################################################################
#                       determiner C1 et S1  => fin
###############################################################################

###############################################################################
#                  determiner la compression P1 et P2  => debut
###############################################################################
def pi1_pi2(sommet_a_corr, 
            augm_sommet_a_corrs, 
            C1_S1_s, 
            number_items_pi1_pi2):
    """
    format C1_S1_s : [ ({1, 2, 3}, [{1, 2, 3}]), ({1, 2, 3, 4}, [{1, 2, 3}]), ....]
    format augm_sommet_a_corrs : [ ({1,2,3,z},[{1,2,3}]), ... ]
    * on fusionne les 2 listes C1_S1_s et augm_sommet_a_corrs, 
        on a fusion_liste_C1S1_augmz
    * ensuite on fait une combinaison par 2 de la liste 
        fusion_liste_C1S1_augmz tel que 
        - dans chaque tuple, l intersection des 2 items du tuple == {sommet_a_corr} 
        
    number_items_pi1_pi2: "pourcentage" d'elements de l_pi1_pi2 a recuperer ( <=1 )
    """
    fusion_liste_C1S1_augmz = C1_S1_s + augm_sommet_a_corrs;
    l_pi1_pi2 = list()

    # test debut
    if len(fusion_liste_C1S1_augmz) == 1:
        tuple_possible = fusion_liste_C1S1_augmz.pop()
        # tuple_possible = ({'D', 'F', 'A', 'E'}, [{'D', 'F', 'A', 'E'}])
        pi1 = tuple_possible[0]
        l_cont = tuple_possible[1]
        l_pi1_pi2.append( (pi1, set(), l_cont) )
        return l_pi1_pi2;
    
    #for tuple_possible in it.combinations(fusion_liste_C1S1_augmz, 2):
    number_items = math.ceil(len(fusion_liste_C1S1_augmz) \
                                 * number_items_pi1_pi2)
    for tuple_possible in it.islice(
                            it.combinations(fusion_liste_C1S1_augmz, 2), 
                            number_items):
        if tuple_possible[0][0].intersection(tuple_possible[1][0]) == \
            {sommet_a_corr}:
            #tmp_pi1_pi2 = (tuple_possible[0][0],  tuple_possible[1][0], tuple_possible[0][1]+tuple_possible[1][1])
            pi1 = tuple_possible[0][0]; 
            pi2 = tuple_possible[1][0];
            l_cont = tuple_possible[0][1] + tuple_possible[1][1];
            l_pi1_pi2.append( (pi1, pi2, l_cont) );
        else:
            pi1 = tuple_possible[0][0].intersection(tuple_possible[1][0]);
            pi2 = set().union( tuple_possible[0][0].union( tuple_possible[1][0] ) \
                      - pi1).union({sommet_a_corr});
            l_cont = [pi1];
            l_pi1_pi2.append( (pi1, pi2, l_cont) );
            
    #print("l_pi1_pi2: ", len(l_pi1_pi2))
    return l_pi1_pi2

#def pi1_pi2(noeud_z, l_augm_z, l_C1_S1, number_items_pi1_pi2):
#    """
#    format l_C1_S1 : [ ({1, 2, 3}, [{1, 2, 3}]), ({1, 2, 3, 4}, [{1, 2, 3}]), ....]
#    format l_augm_z: [ ({1,2,3,z},[{1,2,3}]), ... ]
#    * on fusionne les 2 listes l_C1_S1 et l_augm_z, on a fusion_liste_C1S1_augmz
#    * ensuite on fait une combinaison par 2 de la liste fusion_liste_C1S1_augmz tel que 
#        - dans chaque tuple, l intersection des 2 items du tuple == {noeud_z} 
#        
#    number_items_pi1_pi2: "pourcentage" d'elements de l_pi1_pi2 a recuperer ( <=1 )
#    """
#    fusion_liste_C1S1_augmz = l_C1_S1 + l_augm_z;
#    #print("noeud traite: ", noeud_z," l_C1_S1: ",len(l_C1_S1)," l_augm_z: ",len(l_augm_z)," fusion_l_C1S1_augm_z: ", len(fusion_liste_C1S1_augmz))
#    l_pi1_pi2 = list()
#
#    # test debut
#    if len(fusion_liste_C1S1_augmz) == 1:
#        tuple_possible = fusion_liste_C1S1_augmz.pop()
#        # tuple_possible = ({'D', 'F', 'A', 'E'}, [{'D', 'F', 'A', 'E'}])
#        pi1 = tuple_possible[0]
#        l_cont = tuple_possible[1]
#        l_pi1_pi2.append( (pi1, set(), l_cont) )
#        return l_pi1_pi2
#    #for tuple_possible in it.combinations(fusion_liste_C1S1_augmz, 2):
#    number_items = math.ceil(len(fusion_liste_C1S1_augmz) * number_items_pi1_pi2)
#    for tuple_possible in it.islice(it.combinations(fusion_liste_C1S1_augmz, 2), number_items):
#        ##print("tuple_possible= ", tuple_possible)
#        if tuple_possible[0][0].intersection(tuple_possible[1][0]) == {noeud_z}:
#            #tmp_pi1_pi2 = (tuple_possible[0][0],  tuple_possible[1][0], tuple_possible[0][1]+tuple_possible[1][1])
#            pi1 = tuple_possible[0][0]; 
#            pi2 = tuple_possible[1][0];
#            l_cont = tuple_possible[0][1]+tuple_possible[1][1];
#            l_pi1_pi2.append( (pi1, pi2, l_cont) )
#        else:
#            pi1 = tuple_possible[0][0].intersection(tuple_possible[1][0]);
#            pi2 = set().union( tuple_possible[0][0].union( tuple_possible[1][0] ) - pi1).union({noeud_z})
#            l_cont = [pi1]
#            l_pi1_pi2.append( (pi1, pi2, l_cont) )
#    #print("l_pi1_pi2: ", len(l_pi1_pi2))
#    return l_pi1_pi2
#    # test fin
###############################################################################
#                  determiner la compression P1 et P2  => fin
###############################################################################

###############################################################################
#                  determiner la compression Pis  => debut
###############################################################################
def pis(sommet_a_corr, 
        voisins_sommet_a_corr, 
        pi1_pi2_s, 
        cliques_couvertures):
    """
    retourne le n-uplet (pi1,pi2,pis,[cliques a supprimer])
    l_pi1_pi2: liste des pi1_pi2 de la forme 
                [({"a","b"}, {"c","d"},[{'A','B','C'}]), ... ]
    """
    l_4uplet_pi1_pi2_pis_cliqsASupp = list();
    for p1_p2 in pi1_pi2_s:
        pi1 = p1_p2[0];                                                        # pi1 de la forme pi1 = ({1,2,3},[{1,2},{2,3}])
        pi2 = p1_p2[1];
        gammaZ_pi1 = set(); gammaZ_pi2 = set(); Y = set(); pi_s = set();
        gammaZ_pi1 = set(pi1).intersection( voisins_sommet_a_corr );
        gammaZ_pi2 = set(pi2).intersection( voisins_sommet_a_corr );
        Y = gammaZ_pi1.union(gammaZ_pi2);
        pi_s = set(voisins_sommet_a_corr) - Y;
        
        l_cliques_ASupp = list();
        # pour eviter les doublons de cliste a supprimer 
        # et aussi les cliques a suppirmer 
        # n'etant pas dans ens_C_i
        l_cliques_ASupp = list()
        for cliq_a_suppr in p1_p2[2]:
            if cliq_a_suppr in cliques_couvertures \
                and cliq_a_suppr not in l_cliques_ASupp:
                l_cliques_ASupp.append(cliq_a_suppr);
                
        l_4uplet_pi1_pi2_pis_cliqsASupp.append( 
                (pi1, pi2, pi_s, l_cliques_ASupp) 
                );
        
    ##print("l_4uplet_pi1_pi2_pis_cliqsASupp = ", l_4uplet_pi1_pi2_pis_cliqsASupp)
    return l_4uplet_pi1_pi2_pis_cliqsASupp;

#def pis(noeud_z, gamma_z, l_pi1_pi2, ens_C_i):
#    """
#    retourne le n-uplet (pi1,pi2,pis,[cliques a supprimer])
#    l_pi1_pi2: liste des pi1_pi2 de la forme 
#                [({"a","b"}, {"c","d"},[{'A','B','C'}]), ... ]
#    """
#    l_4uplet_pi1_pi2_pis_cliqsASupp = list()
##    #print("===== l_pi1_pi2 = ", l_pi1_pi2)
#    for p1_p2 in l_pi1_pi2:
##        #print("p1_p2 = ", p1_p2)
#        pi1 = p1_p2[0] # pi1 de la forme pi1 = ({1,2,3},[{1,2},{2,3}])
#        pi2 = p1_p2[1]
#        gammaZ_pi1 = set(); gammaZ_pi2 = set(); Y = set(); pi_s = set()
#        gammaZ_pi1 = set(pi1).intersection(gamma_z)
#        gammaZ_pi2 = set(pi2).intersection(gamma_z)
#        Y = gammaZ_pi1.union(gammaZ_pi2)
#        pi_s = set(gamma_z) - Y
#        
#        l_cliques_ASupp = list()
#        # pour eviter les doublons de cliste a supprimer et aussi les cliques a suppirmer 
#        # n'etant pas dans ens_C_i
#        l_cliques_ASupp = list()
#        for cliq_a_suppr in p1_p2[2]:
#            if cliq_a_suppr in ens_C_i and cliq_a_suppr not in l_cliques_ASupp:
#                l_cliques_ASupp.append(cliq_a_suppr)
#        l_4uplet_pi1_pi2_pis_cliqsASupp.append( (pi1, pi2, pi_s, l_cliques_ASupp) )
#        
#    ##print("l_4uplet_pi1_pi2_pis_cliqsASupp = ", l_4uplet_pi1_pi2_pis_cliqsASupp)
#    return l_4uplet_pi1_pi2_pis_cliqsASupp
###############################################################################
#                  determiner la compression Pis  => fin
###############################################################################
    
###############################################################################
#                     supprimer les petites cliques
#                   contenues dans d'autres cliques  => debut
###############################################################################
def supprimer_cliq_couvert_min_solution(min_tuple_solution_0, 
                                        cliques_couvertures):
    """
    supprimer les petites cliques contenues dans d'autres cliques
    A REVOIR
    {1,2,3}.issubset({1,2}) ====> FALSE
    """
    # TODO si cliq < min_tuple_solution_0 => 
    #         cliq.issubset(min_tuple_solution_0) == FALSE  ===> PROBLEME
    tmp_ens_C = cliques_couvertures.copy()
    for cliq in tmp_ens_C:
        if cliq.issubset(min_tuple_solution_0):
            cliques_couvertures.remove(cliq);
    return cliques_couvertures; 
###############################################################################
#                     supprimer les petites cliques
#                   contenues dans d'autres cliques  => fin
###############################################################################

###############################################################################
#            correction d'une permutation de sommets a corriger => debut
###############################################################################
def corriger_une_permutation_sommets_a_corriger(
                    permut_sommets_a_corriger, 
                    dico_correction_cp,
                    dico_parametres_new, 
                    aretes_cliques):
    """
    corriger une permutation de sommets a corriger
    """
    dico_sommets_corriges = dict();
    som_cout_min = 0;
    aretes_cliques = list( map(set, aretes_cliques) );
    dico_correction_cp["aretes_Ec"] = list( map(
                                              set, 
                                              dico_correction_cp["aretes_Ec"]
                                              ) 
                                            );
    for num_sommet_a_corr, sommet_a_corr in \
        enumerate(permut_sommets_a_corriger):
            
        voisins_sommet_a_corr = dico_correction_cp["dico_gamma_sommets"]\
                                                    [sommet_a_corr];
        clique_sommet_a_corrs = dico_correction_cp\
                                    ["sommets_par_cliqs_avec_aretes"]\
                                    [sommet_a_corr];
        
        clique_contracts = cliques_contractables(
                            sommet_a_corr,
                            clique_sommet_a_corrs,
                            dico_correction_cp["C"],
                            dico_correction_cp["aretes_Ec"],
                            aretes_cliques);
                
        s_sommet_a_corr = S_sommet(sommet_a_corr,
                       voisins_sommet_a_corr,
                       dico_correction_cp["C"],
                       aretes_cliques);
                       
        augm_sommet_a_corrs = augmentation_sommet(
                                sommet_a_corr, 
                                voisins_sommet_a_corr, 
                                clique_sommet_a_corrs, 
                                s_sommet_a_corr, 
                                dico_correction_cp["C"], 
                                aretes_cliques, 
                                2);
        C1_S1_s = determiner_C1_S1(sommet_a_corr, 
                            clique_contracts, 
                            s_sommet_a_corr, 
                            dico_correction_cp["C"]);
        pi1_pi2_s = pi1_pi2(sommet_a_corr, 
                            augm_sommet_a_corrs, 
                            C1_S1_s, 
                            dico_parametres_new["number_items_pi1_pi2"]);
        pi1_pi2_pis_cliqsASupp_s = pis(sommet_a_corr, 
                                       voisins_sommet_a_corr, 
                                       pi1_pi2_s, 
                                       dico_correction_cp["C"]);
                                       
        min_tuple_solution = None; min_cout = 0;
        min_cout, \
        min_tuple_solution, \
        min_arete_a_ajouter_s, \
        min_arete_a_suppr_s, \
        min_cliques_a_suppr_de_C_s = cout_min(
                                        sommet_a_corr, 
                                        pi1_pi2_pis_cliqsASupp_s, 
                                        aretes_cliques, 
                                        dico_correction_cp["aretes_Ec"],
                                        dico_parametres_new[
                                            "critere_selection_compression"]
                                        );
        
        som_cout_min += min_cout if min_cout != None else 0;
        
        if min_tuple_solution != None:
            # supprimer aretes dans liste_arcs_E_C_i
            for arete in min_arete_a_suppr_s:
                if (arete[0],arete[1]) in aretes_cliques:
                    aretes_cliques.remove( (arete[0],arete[1]) )
                    for clique in dico_correction_cp["C"]:
                        if arete[0] in clique and arete[1] in clique:
                            dico_correction_cp["C"].remove(clique);
                if (arete[1],arete[0]) in aretes_cliques:
                    aretes_cliques.remove( (arete[1],arete[0]) )
                    for clique in dico_correction_cp["C"]:
                        if arete[0] in clique and arete[1] in clique:
                            dico_correction_cp["C"].remove(clique);
            
            aretes_cliques.extend(min_arete_a_ajouter_s);
            dico_correction_cp["etats_sommets"][sommet_a_corr] = 1;
        
            # supprimer les cliques de min_cliques_aSupprDe_ens_C dans ens_C ===> Ne pas oublier
            for clique_a_suppr_de_C in min_cliques_a_suppr_de_C_s:
                if  clique_a_suppr_de_C in dico_correction_cp["C"]:
                    dico_correction_cp["C"].remove(clique_a_suppr_de_C)
                    
            
            #ajouter nouvelle clique dans ens_c
            pi1 = min_tuple_solution[0];
            pi2 = min_tuple_solution[1];
            aretes_pi1 = fct_aux.determiner_aretes_cliques(pi1);
            aretes_pi2 = fct_aux.determiner_aretes_cliques(pi2);
            aretes_ajoutees_p1 = fct_aux.aretes_differente(
                                    dico_correction_cp["aretes_Ec"], 
                                    aretes_pi1);
            aretes_ajoutees_p2 = fct_aux.aretes_differente(
                                    dico_correction_cp["aretes_Ec"], 
                                    aretes_pi2);
            aretes_supprimees_ps = min_arete_a_suppr_s;
            if len(pi1) != 0:
                dico_correction_cp["C"]  = supprimer_cliq_couvert_min_solution(
                                            pi1, 
                                            dico_correction_cp["C"]);
                dico_correction_cp["C"].append(pi1);
            if len(pi2) != 0:
                dico_correction_cp["C"]  = supprimer_cliq_couvert_min_solution(
                                            pi2, 
                                            dico_correction_cp["C"]);
                dico_correction_cp["C"].append(pi2);
            
    
            dico_sommets_corriges[(num_sommet_a_corr, sommet_a_corr)] = \
                    {"aretes_ajoutees_p1": aretes_ajoutees_p1,
                     "aretes_ajoutees_p2": aretes_ajoutees_p2,
                     "aretes_supprimees": aretes_supprimees_ps};
            dico_sommets_corriges["som_cout_min"] = som_cout_min;
                  
    return dico_correction_cp, dico_sommets_corriges;
        
###############################################################################
#            correction d'une permutation de sommets a corriger => fin
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
    
    dico_correction["number_items_pi1_pi2"] = \
                                dico_parametres_new["number_items_pi1_pi2"];
    
    print("sommets_a_corriger={} = {}".format(
            len(sommets_a_corriger), sommets_a_corriger));
    
    # liste de combinaisons de sommets a traiter
    permut_sommets_a_corriger_s = [];
    for i in dico_parametres_new["number_permutations_nodes_1"]:
        list_i = random.shuffle(sommets_a_corriger);
        if list_i not in permut_sommets_a_corriger_s:
            permut_sommets_a_corriger_s.append(list_i);
            
    aretes_cliques = fct_aux.determiner_aretes_cliques(
                            dico_correction["C"])
    
    # correction de chaque combi_sommets_a_corriger
    dico_permut = dict();
    for permut_sommets_a_corriger in permut_sommets_a_corriger_s:
        dico_permut[frozenset(permut_sommets_a_corriger)] = \
            corriger_une_permutation_sommets_a_corriger(
                    permut_sommets_a_corriger, 
                    dico_correction.copy(),
                    dico_parametres_new, 
                    aretes_cliques);
                    
    # determination de la meilleure pemutation cad DC ou DH minimale
    dico_best_permut, dico_sommets_corriges = \
                best_permutation(dico_permut, 
                                 dico_correction["aretes_matE"], 
                                 dico_correction["aretes_Ec"]);                 # aretes_Ec = aretes_matE_k_alpha
                                        
    # retourner le resultat
    return dico_best_permut, dico_sommets_corriges;
    
###############################################################################
#               correction aleatoire des sommets_1 => fin
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