#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 12:10:52 2019

@author: willy
"""

import networkx as nx;
import itertools as it;

import clique_max as clique;
import fonctions_auxiliaires as fct_aux;

from operator import itemgetter;

############################ Constantes ---> debut ############################
NBRE_CLIQUES_POSSIBLES = 100;
############################ Constantes ---> fin ##############################


###############################################################################
#                      fonctions annexes ==> debut
###############################################################################
################## ambiguite et graphes doubles DEBUT #########################
def is_isomorphe_graphe_double(aretes_matE_k_alpha):
    """
    but: determiner l'isomorphisme entre matP et 
         les graphes doubles predefinis dans cette fonction
    retourne :
        True s'il ya isomorphisme 
        False s'il n'y a pas d'isomorphisme
    """
    A = nx.Graph()
    B = nx.Graph()
    C = nx.Graph()
    D = nx.Graph()
    edgeA = [(1,2),(1,3),(2,3)]
    edgeB = [(1,2),(1,3),(2,3),(4,2),(4,3)]
    edgeC = [(3,1),(3,2),(3,4),(3,5),(1,2),(1,4),(5,4),(5,2)]
    edgeD = [(1,2),(1,3),(1,4),(1,5),(6,2),(6,3),(6,4),(6,5),(2,3),(2,5),(4,5)]
    
    A.add_edges_from(edgeA)
    B.add_edges_from(edgeB)
    C.add_edges_from(edgeC)
    D.add_edges_from(edgeD)
    
    ### liste arcs matP et conversion nx.Graph()
    matP = nx.Graph()
    matP.add_edges_from(aretes_matE_k_alpha)
    
    
    if nx.is_isomorphic(matP,A) == True or \
        nx.is_isomorphic(matP,B) == True or \
        nx.is_isomorphic(matP,C) == True or \
        nx.is_isomorphic(matP,D) == True :
            return True
        
    return False
################## ambiguite et graphes doubles FIN ###########################

###################### selectionner un sommet DEBUT ###########################
def selectionner_sommet(etats_sommets, dico_gamma_noeud):
    """
    selectionner un sommet dans lequel son etat est soit 0 ou 3 
    selon la strategie suivante:
        * sommet de degre minimum.
    """
    sommets_0 = [(sommet, dico_gamma_noeud[sommet][0]) \
                 for sommet, etat in etats_sommets.items() if etat == 0];
    sommets_3 = [(sommet, dico_gamma_noeud[sommet][0]) \
                 for sommet, etat in etats_sommets.items() if etat == 3];
    
    sommet_u = None;
    if sommets_0:
        sommet_u = max(sommets_0, key=itemgetter(1));
    elif sommets_3:
        sommet_u = max(sommets_3, key=itemgetter(1));
    return sommet_u[0];
###################### selectionner un sommet fin ###########################
   
######################## verifier de cliques DEBUT ############################
def verifier_cliques(cliques):
    """
    verifier si tous les sommets d'une clique de cliques concourent a un sommet.
    """
    cliques_possibles = [];
    for clique_ in cliques:
        aretes = list(map(lambda x:set(x.split("_")), clique_));
        sommet_commun = set.intersection(*aretes);
        
        if sommet_commun and len(sommet_commun) == 1:
            cliques_possibles.append(clique_);
            
    return cliques_possibles;
######################## verifier de cliques fin ##############################
    
######################## trouver les cliques coherentes DEBUT #################
def trouver_cliques_coherentes(cliques_possibles):
    """
    verifier si 2 cliques ont un sommet commun.
    retourne les cliques coherentes en une liste de tuple ordonnee 
    par taille decroissante.
    """
    cliques_possibles = map(set, cliques_possibles);
    
    cliques_coherentes = [];
    tuples_cliques = it.combinations(cliques_possibles, 2);
    boolean = True; cpt_treated_tuple = 0;
    
    while boolean:
        try :
            tuple_clique = next(tuples_cliques);
            sommets_commun = set.intersection(*tuple_clique);
        
            if len(sommets_commun) == 1:
                cliques_coherentes.append(tuple_clique);
            elif len(sommets_commun) == 2:
                # TODO
                # cliques fragmentees : nous sommmes dans le cas de sous-graphes:
                # cervolant ou len(tup[0])>=3 ou len(tup[1])=3 
                pass
            if cpt_treated_tuple > NBRE_CLIQUES_POSSIBLES:
                boolean = False;
            cpt_treated_tuple += 1;
        except StopIteration:
            boolean = False;
    
    # ordonner la liste des tuples selon la taille des elements de chaque tuple.        
    if cliques_coherentes:
        sorted(cliques_coherentes, 
               key=lambda tup: (len(tup[0]),len(tup[1])), reverse=True)
    
    return cliques_coherentes;
######################## trouver les cliques coherentes fin ###################

#################### supprimer les aretes de C1 et C2 debut ###################
def mise_a_jour_voisinage_sommets(dico_gamma_noeud,
                                  cliques, 
                                  matE_k_alpha, 
                                  aretes_matE_k_alpha):
    """
    supprimer les aretes de C1 et C2 en changeant le voisinage des sommets.
    """
    aretes_supp = [];
    for C in cliques:
        if C is not None:
            aretes = set(it.combinations(C, 2));
            aretes_supp.extend(aretes);
            for arete in aretes:
                matE_k_alpha.loc[arete[0],arete[1]] = 0;
                matE_k_alpha.loc[arete[1],arete[0]] = 0;
    aretes_matE_k_alpha = fct_aux.liste_arcs(matE_k_alpha);
    dico_gamma_noeud = fct_aux.gamma_noeud(matE_k_alpha, aretes_matE_k_alpha);
    
    return dico_gamma_noeud, list(aretes_matE_k_alpha);
###################### supprimer les aretes de C1 et C2 fin ###################
    
###############################################################################
#                      fonctions annexes ==> fin
###############################################################################

###############################################################################
#                      couverture en cliques ==> debut
###############################################################################
def couverture_en_cliques(etats_sommets, 
                         dico_gamma_noeud, 
                         aretes_matE_k_alpha, 
                         matE_k_alpha, 
                         dico_ver, 
                         arguments_MAJ):
    """
    algo de couverture en cliques cad couvrir chaque par une clique et 
    marquer les sommets du graphe matE_k_alpha.
    """
    C = list();
    ordre_noeuds_traites = list();
    dico_sommets_par_cliqs = dict();
    
    while {0,3}.intersection(etats_sommets.values()):
        print("Y1")
        sommet_u = selectionner_sommet(etats_sommets, dico_gamma_noeud);
        if sommet_u is None:
            print("ERREUR sommet_u = None, EXIST sommet ETATS 0,3");
            return {"C":C, "etats_sommets": etats_sommets, 
                    "aretes_restantes": aretes_matE_k_alpha, 
                    "ordre_noeuds_traites": ordre_noeuds_traites, 
                    "sommets_par_cliqs": dico_sommets_par_cliqs};
        
        print("Y2 sommet_u={}".format(sommet_u))
        voisins_u = dico_gamma_noeud[sommet_u][1]
        cliques = clique.find_clique(matE_k_alpha, 
                                     voisins_u.union({sommet_u}), 
                                     []);
        print("Y3")
        
        #### TODO test voisinage A EFFACER ---> debut
        item0 = sommet_u.split("_")[0]; item1 = sommet_u.split("_")[1];
        ens0, ens1 = set(), set();
        for vois in voisins_u:
            if item0 in vois.split("_"):
                ens0.add(vois);
            if item1 in vois.split("_"):
                ens1.add(vois);
        print("voisins_u={}, ens0={}, ens1={}".format(len(voisins_u), ens0, ens1));
        #### test voisinage A EFFACER ---> fin
        
        cliques_possibles = verifier_cliques(cliques);
        print("Y4 cliques_possibles={}".format(cliques_possibles))
        fr_cliques = map(frozenset, cliques);
        fr_cliques_possibles = map(frozenset, cliques_possibles);
        print("Y4bis cliques_not_possibles={}".format( set(fr_cliques) - set(fr_cliques_possibles) ));
        print("cliques ={}".format(cliques)) if not cliques_possibles else None;
        
        C1, C2 = None, None;
        if not cliques_possibles:
            etats_sommets[sommet_u] = -1;
        elif len(cliques_possibles) == 1:
            C1 = cliques_possibles[0];
        else:
            cliques_coherentes = trouver_cliques_coherentes(cliques_possibles);
            print("Y41 cliques_coherentes={}".format( len(cliques_coherentes) ))
            if not cliques_coherentes:
                etats_sommets[sommet_u] = -1;
            else:
                C1, C2 = cliques_coherentes[0][0], cliques_coherentes[0][1];
        print("Y5  C1={}, C2={}".format(C1, C2))
        
        # mettre les etats des sommets
        # TODO probleme quand C1,C2=None, None
        if C1 is None and C2 is None:
            etats_sommets[sommet_u] = -1;
        elif etats_sommets[sommet_u] == 3 and C2 is not None:
            etats_sommets[sommet_u] = -1;
        else:
            for sommet_w in voisins_u:
                voisins_w = dico_gamma_noeud[sommet_w][1] - {sommet_u};
                C1_C2 = C1.union(C2) if C2 is not None else C1;
                if C1_C2.union(voisins_w) - C1_C2:                              # il existe un sommet appartenant a C1.union(C2)
                    if etats_sommets[sommet_w] == 3:
                        etats_sommets[sommet_w] = -1;
                    elif etats_sommets[sommet_w] == 0:
                        etats_sommets[sommet_w] = 3;
                else:
                    if etats_sommets[sommet_w] == 3:
                        etats_sommets[sommet_w] = 2;
                    elif etats_sommets[sommet_w] == 0:
                        etats_sommets[sommet_w] = 1;
        print("Y6")
        if etats_sommets[sommet_u] == 0:
            if C2 is None:
                etats_sommets[sommet_u] = 1;
            else:
                etats_sommets[sommet_u] = 2;
        elif etats_sommets[sommet_u] != -1:
             etats_sommets[sommet_u] = 2; 
        
        # ajouter C1 et C2 dans C
        print("Y7")
        C.append(C1) if C1 is not None else None;
        C.append(C2) if C2 is not None else None;
        
        # mise a jour des voisinages des sommets.
        print("Y8")
        dico_gamma_noeud, aretes_matE_k_alpha = \
            mise_a_jour_voisinage_sommets(
                                dico_gamma_noeud,
                                [C1,C2], 
                                matE_k_alpha,
                                aretes_matE_k_alpha);
        
        print("Y9 aretes_matE_k_alpha={}".format(len(aretes_matE_k_alpha)))
        ordre_noeuds_traites.append(sommet_u);
        pass # end while
        
    dico_sommets_par_cliqs = fct_aux.couverture_par_sommets(
                                sommets_matE = etats_sommets.keys(),
                                C = C)
    return {"C":C, "etats_sommets": etats_sommets, 
            "aretes_restantes": aretes_matE_k_alpha, 
            "ordre_noeuds_traites": ordre_noeuds_traites, 
            "sommets_par_cliqs": dico_sommets_par_cliqs};

###############################################################################
#                      couverture en cliques ==> fin
###############################################################################
    
###############################################################################
#                      algorithme de decomposition ==> debut
###############################################################################
def algo_decomposition_en_cliques(matE_k_alpha, 
                          dico_sommet_arete, 
                          seuil_U=10, 
                          epsilon=0.75,
                          chemin_datasets="", 
                          chemin_matrices="",
                          ascendant_1=True, simulation=True, 
                          dico_proba_cases=dict(),
                          dico_parametres_new=dict()
                          ):
    """
    obtenir la decomposition en cliques des sommets du graphe matE_k_alpha
    un sommet doit appartenir au plus a 2 cliques.
    une arete doit appartenir a 1 clique.
    chaque sommet a 5 etats {0,1,2,3,-1}
    """
    #initialisation cliq et ver et C
    etats_sommets = dict(); dico_ver = dict();
    for sommet in matE_k_alpha.columns:   # nbre de noeuds dans le graphe
        etats_sommets[sommet] = 0; dico_ver[sommet] = 0;
    
    # fusion des datasets 
    liste_grandeurs = []; #fct_aux.liste_grandeurs(chemin_datasets)
    arguments_MAJ = {"dico_sommet_arete": dico_sommet_arete, 
                     "df_fusion": dict(), 
                     "seuil_U": seuil_U, "epsilon":epsilon, 
                     "chemin_dataset": chemin_datasets,
                     "simulation": simulation, "grandeurs": liste_grandeurs}
          
    # copy E0 <- Ec
    aretes_matE_k_alpha = fct_aux.liste_arcs(matE_k_alpha);
    dico_gamma_noeud = fct_aux.gamma_noeud(matE_k_alpha, aretes_matE_k_alpha) # {"2":[3,{"1","3","4"}],....}
    
    
    if is_isomorphe_graphe_double(aretes_matE_k_alpha):
        """
        QU'EST CE un GRAPHE DOUBLE : graphe avec 2 couvertures.
        """
        #print("le traiter avec Verif_correl ou ORACLE")
        return {"C":list(), "etats_sommets": etats_sommets, 
                "aretes_restantes": aretes_matE_k_alpha, 
                "ordre_noeuds_traites": list(), 
                "sommets_par_cliqs": dict()}
    else:
        dico_couverture = couverture_en_cliques(etats_sommets, 
                                                 dico_gamma_noeud, 
                                                 aretes_matE_k_alpha, 
                                                 matE_k_alpha.copy(), 
                                                 dico_ver, 
                                                 arguments_MAJ)
     
    return dico_couverture;
###############################################################################
#                      algorithme de decomposition ==> fin
###############################################################################
