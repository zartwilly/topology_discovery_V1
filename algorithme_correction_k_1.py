#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 11:30:14 2019

@author: willy
"""
import sys;
sys.path.append('../');
import fonctions_auxiliaires as fct_aux;

###############################################################################
#           correction_pi1_pi2 sans suppression d aretes --> debut
###############################################################################
def correction_pi1_pi2_sans_ps(sommet_1, sommets_1, dico_correction):
    """
    corriger les cliques au voisinage de sommet_1.
    """
    C_ab = []; 
    C_x = []; 
    C_xb = []; 
    
    sommets_1 = set(sommets_1) - set({sommet_1});
    
    # identifier les C_ab, C_xb, C_x dans 'sommets_par_cliqs_avec_aretes'
    for clique in dico_correction['sommets_par_cliqs_avec_aretes'][sommet_1]:
        if len(clique) <= 2:
            C_xb = clique;
        elif len(clique) > 2 and len(clique.intersection(sommets_1)) == 0:
            C_x = clique;
        elif len(clique) > 2:
            C_ab = clique;
    print("C_xb={}, C_x={}, C_ab={}".format(C_xb, C_x, C_ab));
    
    # former les cliques pi_1, pi_2
    pi_1 = set(C_ab).union(C_xb); 
    pi_2 = set(C_x);
    print("pi_1={}, pi_2={}".format(pi_1, pi_2));
    
    # supprimer les aretes ou cliques, sous ensemble de pi_1, pi_2
    cliques_C = dico_correction["C"];
    cliques_a_del = [];
    for clique in cliques_C:
        if clique.issubset(pi_1) or pi_1.issubset(clique) or \
            clique.issubset(pi_2) or pi_2.issubset(clique):
            cliques_a_del.append(clique);
    
    cliques_a_del.append(C_x) if len(C_x) > 0 else None;
    cliques_a_del.append(C_ab) if len(C_ab) > 0 else None;
    cliques_a_del.append(C_xb) if len(C_xb) > 0 else None;
    
    print("cliques_a_del={}, cliques_C_avant_del={}".format(
            len(cliques_a_del),len(cliques_C)))
    for clique_a_del in cliques_a_del:
        id_cliq_a_del = cliques_C.index(clique_a_del) \
                        if clique_a_del in cliques_C else None;
        cliques_C.pop(id_cliq_a_del) if id_cliq_a_del != None else None;
    print("cliques_apres_del={}".format(len(cliques_C)));
        
    return cliques_C, pi_1, pi_2;
###############################################################################
#           correction_pi1_pi2 sans suppression d aretes --> fin
###############################################################################

###############################################################################
#                       rechercher_sommet_b_B_ab --> debut
###############################################################################
def rechercher_sommet_b_B_ab(B_ab, dico_correction):
    """
    rechercher le sommet b au sein de l'ensemble B_ab des sommets a -1.
    """
    sommet_b = None;
    id_sommet_b_s = {};
    for id_som_x, sommet_x in enumerate(B_ab):
        for id_som_y, sommet_y in enumerate(B_ab[id_som_x+1:]):
            cliques_x = dico_correction["sommets_par_cliqs_avec_aretes"][sommet_x];
            cliques_y = dico_correction["sommets_par_cliqs_avec_aretes"][sommet_y];
            fr_cliqs_x = set(map(frozenset, cliques_x));
            fr_cliqs_y = set(map(frozenset, cliques_y));
            print("fr_cliqs_x={}, fr_cliqs_y={}".format(fr_cliqs_x, fr_cliqs_y))
            print("REC1 INT={} \n".format(fr_cliqs_x.intersection(fr_cliqs_y)))
            if len(fr_cliqs_x.intersection(fr_cliqs_y)) == 0:
                id_sommet_b_s.add(id_som_x);
            print("REC2")
            
    print("REC3 id_sommet_b_s={}".format(id_sommet_b_s))
    sommet_b = B_ab[list(id_sommet_b_s)[0]] if id_sommet_b_s else None;
    B_ab = B_ab.pop(list(id_sommet_b_s)[0]) if id_sommet_b_s else B_ab;
    print("REC4")
    
    print('sommet_b={}'.format(sommet_b))
    return B_ab, sommet_b;
###############################################################################
#                       rechercher_sommet_b_B_ab --> fin
###############################################################################

###############################################################################
#                       correction_cliques --> debut
###############################################################################
def correction_cliques(dico_correction, dico_gamma_noeud, dico_parametres_new):
    """
    
    """
    print("cliques_sans_aretes ={}".format(len(dico_correction["C"])))
    dico_correction["C"] = dico_correction["C"] \
                            + list(map(set, dico_correction['aretes_restantes']));
    print("cliques_avec_aretes ={}".format(len(dico_correction["C"])))
    
    dico_correction["sommets_par_cliqs_avec_aretes"] = \
                fct_aux.couverture_par_sommets(
                    sommets_matE=list(dico_correction["etats_sommets"].keys()),
                    C=dico_correction["C"]);
                        
    # construction des ensembles B_a,b et U_a,b
    ## selection de sommets -1 
    sommets_1 = [k for k, v in dico_correction["etats_sommets"].items() 
                    if v == -1]
    B_ab, U_ab, not_B_U_ab = list(), list(), list();
    print("R1")
    for sommet_1 in sommets_1:
        if len(dico_correction['sommets_par_cliqs_avec_aretes'][sommet_1]) == 3:
            B_ab.append(sommet_1);
        elif len(dico_correction['sommets_par_cliqs_avec_aretes'][sommet_1]) == 2:
            U_ab.append(sommet_1);
        else:
            not_B_U_ab.append(sommet_1);
    
    print("R2")        
    B_ab, sommet_b = rechercher_sommet_b_B_ab(B_ab, dico_correction);
            
    print("sommets_1={}, B_ab={}, U_ab={}, not_B_U_ab={}".format(
            len(sommets_1), len(B_ab), len(U_ab), len(not_B_U_ab)))
    
    print("R3")
    if len(B_ab) == 0 and len(U_ab) == 1:
        print("NOT POSSIBLE TO FILL REMOVED EDGE: Graph is a LINEGRAPH")
    elif len(B_ab) == 1 and len(U_ab) == 0:
        print("R311")
        sommet_1 = B_ab[0];
        cliques_C, pi_1, pi_2 = correction_pi1_pi2_sans_ps(sommet_1, 
                                    sommets_1, dico_correction)
        
        #mise a jour cliques C dans dico_correction
        dico_correction["C"] = cliques_C;
        # ajout cliques pi_x dans C
        dico_correction["C"].extend([pi_1, pi_2]);
        print("R312")
        
    elif (len(B_ab) == 1 and len(U_ab) == 0) or \
         (len(B_ab) >= 1 and len(U_ab) >= 0):
        print("R321")
        sommet_1 = B_ab[0];
        cliques_C, pi_1, pi_2 = correction_pi1_pi2_sans_ps(sommet_1, 
                                    sommets_1, dico_correction)
        
        #mise a jour cliques C dans dico_correction
        dico_correction["C"] = cliques_C;
        # ajout cliques pi_x dans C
        dico_correction["C"].extend([pi_1, pi_2]);
        print("R322")
        
    elif (len(B_ab) == 1 and len(U_ab) == 1):
        print("R331")
        sommet_x = B_ab[0];
        sommet_z = U_ab[0];
        couv_x = dico_correction['sommets_par_cliqs_avec_aretes'][sommet_x];
        couv_b = dico_correction['sommets_par_cliqs_avec_aretes'][sommet_b];
        
        voisins_b = dico_gamma_noeud[sommet_b][1] \
                        - {sommet_x} \
                        - {sommet_z};
        voisins_x = dico_gamma_noeud[sommet_x][1] \
                        - {sommet_b} \
                        - {sommet_z}
        
        if len(voisins_b) == 0 and len(voisins_x) == 0:
            if len(couv_x.intersection(couv_b)) != 0:
                index_clique_xb = dico_correction["C"].index({sommet_x, sommet_b});
                dico_correction["C"].pop(index_clique_xb);
            
        elif len(voisins_b) == 0 and len(voisins_x) != 0:
            sommet_1 = B_ab[0];
            cliques_C, pi_1, pi_2 = correction_pi1_pi2_sans_ps(sommet_1, 
                                    sommets_1, dico_correction)
        
            #mise a jour cliques C dans dico_correction
            dico_correction["C"] = cliques_C;
            # ajout cliques pi_x dans C
            dico_correction["C"].extend([pi_1, pi_2]);
            
        elif len(voisins_b) != 0 and len(voisins_x) == 0:
            print("CORRECTION NOT POSSIBLE because VOIS_X = EMPTY"+
                  " and VOIS_b NOT EMPTY")
            
        elif len(voisins_b) != 0 and len(voisins_x) != 0:
            sommet_1 = B_ab[0];
            cliques_C, pi_1, pi_2 = correction_pi1_pi2_sans_ps(sommet_1, 
                                    sommets_1, dico_correction)
        
            #mise a jour cliques C dans dico_correction
            dico_correction["C"] = cliques_C;
            # ajout cliques pi_x dans C
            dico_correction["C"].extend([pi_1, pi_2]);
            
        print("R332")
    
    dico_correction["sommets_par_cliqs_avec_aretes"] = \
                fct_aux.couverture_par_sommets(
                    sommets_matE=list(dico_correction["etats_sommets"].keys()),
                    C=dico_correction["C"]);    
        
    return dico_correction;
    pass
###############################################################################
#                       correction_cliques --> fin
###############################################################################