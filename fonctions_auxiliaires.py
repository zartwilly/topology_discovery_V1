#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 12:22:22 2019

ce fichier contient toutes les fonctions secondaires utilises 
dans tous les fichiers .py

@author: willy
"""
import random;
import numpy as np;
import pandas as pd;
import itertools as it;


import genererMatA as geneMatA;
import generations_mesures as mesures;

from scipy.stats import truncnorm;

###############################################################################
#           identifier les arcs du graphe ==> debut
###############################################################################
def nommage_arcs(matA):
    """
    le but est de nommer les arcs entre 2 sommets de la maniere suivante:
        dico_dual_arc_sommet= {"A":("u","v"),"B":("u","w"),"C":("w","v"),"D":("z","w")}

    dans le cas ou dico n'est pas necessaire alors faire ceci 
    dual_arcs_sommets = list(map(lambda u_v: ("_").join([u_v[0],u_v[1]]),list_arcs_))
    """
    dico_dual_arc_sommet = dict();
    liste_arcs_ = []
    if type(matA) is pd.DataFrame:
        liste_arcs_ = liste_arcs(matA, oriented=True);
    elif type(matA) is list:
        liste_arcs_ = matA
    else:
        print("type matA UNKNOWN");
    for cpt, arc in enumerate(liste_arcs_):
        nom_arc = "_".join([arc[0], arc[1]]) \
            if int(arc[0]) <= int(arc[1]) else "_".join([arc[1], arc[0]]);
        dico_dual_arc_sommet[nom_arc] = arc;
    return dico_dual_arc_sommet;
###############################################################################
#           identifier les arcs du graphe ==> fin
###############################################################################

###############################################################################
#           liste des arcs/ aretes du graphe ==> debut
###############################################################################
def liste_arcs(mat, oriented=False, val=1):
    """ retourne la liste des 
        * arcs si oriented = True ou 
        * aretes si oriented = False 
        d'un graphe. 
        si val == 1: liste des arcs 
        si val == 0: liste des sommets non adjacents
    """
    res = list()
    if not oriented:
        res = list(filter(lambda row_col: mat.loc[row_col[0],row_col[1]]==val,
                 it.combinations(mat.columns.tolist(),2)))
    else:
        res = list(filter(
                lambda row_col: mat.loc[row_col[0],row_col[1]]==val or
                        mat.loc[row_col[1],row_col[0]]==val,
                 it.combinations(mat.columns.tolist(),2)))
    return res;
###############################################################################
#           liste des arcs/ aretes du graphe ==> fin
###############################################################################
   

    
"""
creation au graphe racine et de son linegraphe
"""
###############################################################################
#           creation du linegraphe  ==> debut
###############################################################################
def creer_mat_LG(arcs_or_aretes) :
    """ Methode qui determine la matrice du line graphe de mat_GR a partir de
        la liste des arcs ou aretes de mat_GR.
    """
    mat_LG = None;
    
    aretes = map(set, arcs_or_aretes)
    
    dico_graphe = dict()
    for (arete0,arete1) in it.combinations(aretes,2) :
        if arete0.intersection(arete1) :

            arete0 = list(arete0); arete1 = list(arete1);
            arete0 = "_".join([arete0[0], arete0[1]]) \
                if int(arete0[0]) <= int(arete0[1]) \
                else "_".join([arete0[1], arete0[0]])
            arete1 = "_".join([arete1[0], arete1[1]]) \
                if int(arete1[0]) <= int(arete1[1]) \
                else "_".join([arete1[1], arete1[0]])
            
            if arete0 not in dico_graphe and arete1 not in dico_graphe :
                dico_graphe[arete0] = [arete1];
                dico_graphe[arete1] = [arete0];
            elif arete0 not in dico_graphe and arete1 in dico_graphe :
                dico_graphe[arete1].append(arete0);
                dico_graphe[arete0] = [arete1];
            elif arete0 in dico_graphe and arete1 not in dico_graphe :
                dico_graphe[arete0].append(arete1);
                dico_graphe[arete1] = [arete0]
            elif arete0 in dico_graphe and arete1 in dico_graphe :
                dico_graphe[arete0].append(arete1);
                dico_graphe[arete1].append(arete0);
    mat_LG = pd.DataFrame(index = dico_graphe.keys(), 
                          columns = dico_graphe.keys());
    for k, vals in dico_graphe.items() :
        for v in vals:
            mat_LG.loc[k,v] = 1
            mat_LG.loc[v,k] = 1
    mat_LG.fillna(value=0, inplace=True);
    
    return mat_LG.astype(int);
###############################################################################
#           creation du linegraphe  ==> fin
###############################################################################

###############################################################################
#           creation du graphe racine et de son linegraphe  ==> debut
###############################################################################
def creer_GR_LG(dim_mat, 
                nbre_lien, 
                chemin_matrices, 
                chemin_datasets,
                nbre_ts,
                effet_joule,
                without_mesures):
    """
    creer un graphe racine et son line-graphe sans/avec generer 
    les mesures de flots dans le graphe selon la variable without_mesures
    si without_mesures=True alors PAS DE GENERATION DE MESURES
    si without_mesures=False alors GENERATION DE MESURES
    """
    matA_GR = geneMatA.genererMatriceA(dim_mat, nbre_lien);
    arcs_or_aretes = liste_arcs(matA_GR, oriented=True);
    matE_GR = creer_mat_LG(arcs_or_aretes);
    dico_sommet_arete = nommage_arcs(arcs_or_aretes);
    if not without_mesures:
        mesures.create_datasets(matA_GR, dico_sommet_arete,
                                chemin_datasets, nbre_ts, effet_joule);
        
    matA_GR.to_csv(chemin_matrices+"df_matA_generer.csv");
    matE_GR.to_csv(chemin_matrices+"matE_LG.csv");
    return matE_GR, matA_GR, dico_sommet_arete;
###############################################################################
#           creation du graphe racine et de son linegraphe  ==> fin
###############################################################################



"""
modification au graphe cree
"""
###############################################################################
#           selection et modification k cases ==> debut
###############################################################################
def delete_aretes(j, aretes_0_1, dico_deleted_add_edges, operation):
    """
    suppression de l'arete a l'indice j puis 
    ajout dans le dico_deleted_add_edges
    soit avec la cle "ajouter" ou soit avec la cle "supprimer".
    """
    if operation == 0:
        # case 0 -> 1: case from 0 to 1
        dico_deleted_add_edges["ajouter"].append(aretes_0_1[j])
        del aretes_0_1[j]   
    elif operation == 1:
        # case 1 -> 0: case from 1 to 0
        dico_deleted_add_edges["supprimer"].append(aretes_0_1[j])
        del aretes_0_1[j] 
    
    return aretes_0_1, dico_deleted_add_edges;  

def ajout_suppression_arete(proba, aretes_0, aretes_1, k_cases):
    """
    selection de cases et modification
    """
    dico_deleted_add_edges = dict();
    dico_deleted_add_edges["ajouter"] = []; 
    dico_deleted_add_edges["supprimer"] = [];
    for k in range(0, k_cases):
        if proba < 0.5:
            # case 0 -> 1
            i = random.choice(range(0, len(aretes_0)))
            aretes_1.append(aretes_0[i])
            aretes_0, dico_deleted_add_edges = delete_aretes(
                                                    i, 
                                                    aretes_0,
                                                    dico_deleted_add_edges, 
                                                    0)
        else:
            # case 1 -> 0
            j = random.choice(range(0, len(aretes_1)))
            aretes_0.append(aretes_1[j])
            aretes_1, dico_deleted_add_edges = delete_aretes(
                                                    j, 
                                                    aretes_1,
                                                    dico_deleted_add_edges, 
                                                    1)
    return aretes_0, aretes_1, dico_deleted_add_edges;
###############################################################################
#           selection et modification k cases ==> fin
###############################################################################
    
###############################################################################
#           modifier k aretes dans le graphe ==> debut
###############################################################################
def modif_k_cases(matE, k_cases, methode, p_correl_seuil):
    """
    changer k cases dans matE tel que les cases
        * 0 passent a 1 ( FAUX POSITIFS )
        * 1 passent a 0 ( FAUX NEGATIFS )
    le but est de reduire le biais cree par la suppression d aretes dans matE
    methode 0: avec proba
    methode 1: avec liste
    p_correl_seuil = 0.5, 0.75, 1 # si p_correl_seuil = 0 => on ajoute que des aretes
                                  # si p_correl_seuil = 1 => on supprime que des aretes
                                  # si p_correl_seuil = 0.5 => on ajoute et supprime des aretes
    explications:
        0           faux neg     faux pos      1
        ---------------------|----------------------------> p_correl
        0->0        1->0     |    0->1       1->1
                             |
         suppresion aretes   |     ajout aretes
                             |
                        p_correl_seuil
    
    test methode:
        dimMat = 5;nbre_lien = 5;nbre_ts = 10
        effet_joule = 0;epsilon = 0.75;test = "FINI";
        chemin_datasets = "data/G1/datasets/";chemin_matrices = "data/G1/matrices/";
        methode = 0
        matE, df_matA, dico_arc = simu50.matriceE(dimMat, nbre_lien, chemin_matrices,\
                                                  chemin_datasets, nbre_ts, epsilon, \
                                                  effet_joule, test)
        matE_, dico_deleted_add_edges= modif_k_cases(matE.copy(), 2, methode)
        DH, l_DH = simu50.distance_hamming(fct_aux.liste_arcs(matE),fct_aux.liste_arcs(matE_))
        print("DH = ", DH," l_dh = ", l_DH," dico_deleted_add_edges = ", dico_deleted_add_edges.values())
    """
    
    aretes_0 = liste_arcs(mat=matE, oriented=False, val=0);
    aretes_1 = liste_arcs(mat=matE, oriented=False, val=1);
    
    if methode == 0:
        # methode avec biais
        dico_deleted_add_edges = dict();
        dico_deleted_add_edges["ajouter"] = []; 
        dico_deleted_add_edges["supprimer"] = [];
        
        for k in range(0, k_cases):
            proba = random.random();
        
            m0 = len(aretes_0);
            m1 = pow(matE.shape[0],2) - m0;    
                     
            if proba < p_correl_seuil:
                # case 0 -> 1
                p0 = 1/m0;
                bool = True;
                while bool:
                    arete = random.choice(aretes_0)
                    if random.random() > p0 and \
                        arete not in dico_deleted_add_edges["ajouter"] and \
                        arete not in dico_deleted_add_edges["supprimer"] :
                        matE.loc[arete[0]][arete[1]] = 1; 
                        matE.loc[arete[1]][arete[0]] = 1;
                        dico_deleted_add_edges["ajouter"].append(arete) 
                        bool = False
                        break;
#                for arete in aretes_0 :
#                    if random.random() > p0 and arete not in dico_deleted_add_edges["ajouter"] \
#                        and arete not in dico_deleted_add_edges["supprimer"] :
#                        matE.loc[arete[0]][arete[1]] = 1; matE.loc[arete[1]][arete[0]] = 1;
#                        dico_deleted_add_edges["ajouter"].append(arete)   
#                        break;
            else:
                # case 1 -> 0
                p1 = 1/m1;
                bool = True;
                while bool:
                    arete = random.choice(aretes_1)
                    if random.random() > p1 and \
                        arete not in dico_deleted_add_edges["ajouter"] and \
                        arete not in dico_deleted_add_edges["supprimer"]:
                        matE.loc[arete[0]][arete[1]] = 0; 
                        matE.loc[arete[1]][arete[0]] = 0;
                        dico_deleted_add_edges["supprimer"].append(arete)   
                        bool = False;
                        break;   
#                for arete in aretes_1: 
#                    if random.random() > p1 and arete not in dico_deleted_add_edges["ajouter"] \
#                        and arete not in dico_deleted_add_edges["supprimer"]:
#                        matE.loc[arete[0]][arete[1]] = 0; matE.loc[arete[1]][arete[0]] = 0;
#                        dico_deleted_add_edges["supprimer"].append(arete)   
#                        break;
        return matE, dico_deleted_add_edges;
    else:
        # methode avec liste
        proba = random.random();
        dico_deleted_add_edges = dict();
        aretes_0, aretes_1, dico_deleted_add_edges = \
        ajout_suppression_arete(proba, aretes_0, aretes_1, k_cases)
        
        for key, aretes in dico_deleted_add_edges.items():
            if key == "ajouter":
                for arete in aretes:
                    matE.loc[arete[0]][arete[1]] = 1; 
                    matE.loc[arete[1]][arete[0]] = 1;
            elif key == "supprimer":
                for arete in aretes:
                    matE.loc[arete[0]][arete[1]] = 0; 
                    matE.loc[arete[1]][arete[0]] = 0; 
                    
        return matE, dico_deleted_add_edges;
###############################################################################
#           modifier k aretes dans le graphe ==> debut
###############################################################################

###############################################################################
#                ajouter la proba de chaque case de matE ====> debut     
###############################################################################
def loi_proba(x,s):
    return pow(x-s, 2);
def loi_de_probalibilite_old(debut_prob, fin_proba, nbre_correlations, correl_seuil):
    """
    genere des valeurs compris [debut_prob, fin_proba] et leur proba associe
    """
    correlations = np.linspace(debut_prob, fin_proba, nbre_correlations)
    proba_correlations =  [loi_proba(x,correl_seuil) for x in correlations]
    # Normalising to 1.0
    proba_correlations /= np.sum(proba_correlations)
    return round(np.random.choice(correlations, 1, p=proba_correlations)[0], 3)
    
def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)
    
def loi_de_probalibilite(debut_proba, fin_proba, nbre_correlations, correl_seuil):
    """
    generer des valeurs de correlations entre debut_proba et fin_proba.
    """
    mean = (fin_proba + debut_proba)/2; 
    sd = abs((fin_proba - debut_proba)/(2*4)) \
            if fin_proba != debut_proba else abs((fin_proba)/(2*4))
    correlations  = get_truncated_normal(mean*100, 
                                         sd*100, 
                                         low=1, 
                                         upp = nbre_correlations*100);
    return np.random.choice(correlations.rvs(nbre_correlations)/100);

def ajouter_proba_matE(matE, 
                       dico_deleted_add_edges, 
                       loi_stats, 
                       p_correl, 
                       correl_seuil=0.8):
    """
    le but est de definir des probabilites pour chaque case de matE pour simuler 
    la matrice de correlation obtenue apres calcul. les probas (correls) sont les suivantes:
        
        * case 0 -----> 0: proba entre 0 - 0.5 (0, correl_seuil-p_correl-0.01)           ==> vrai negatifs
        * case 1 -----> 0: proba entre 0.6 - 0.79 ( correl_seuil-p_correl, correl_seuil) ==> faux negatifs
        * case 0 -----> 1: proba entre p_correl et correl_seuil-0.001 (p_correl,correl_seuil-0.01)==> faux positifs
        * case 1 -----> 1: proba entre 0.8 - 1 (correl_seuil, 1)             ==> vrai positifs
    NB: X = correl_seuil
    dico_deleted_add_edges = {"ajouter":[(a,b),...],"supprimer":[(c,d),...]}
    loi_stats = loi uniforme, de poisson
    """
    dico_proba_case = dict(); nbre_correlations = 250;
    colonnes = matE.columns.tolist();
    for id_row, row in enumerate(colonnes):
        for col in colonnes[id_row+1:]:
            if matE.loc[row,col] == 0:
                if (row,col) in dico_deleted_add_edges["supprimer"] or \
                    (col,row) in dico_deleted_add_edges["supprimer"]:
                    # proba entre 0.6 et 0.79 (entre X-0.2 -- X-0.01) ===> faux negatifs (distributions uniforme)
                    if loi_stats == "poisson":
                        dico_proba_case[(row,col)] = \
                            round(random.choice(
                                    np.union1d(np.linspace(0.6,0.644,100), \
                                               np.linspace(0.72,0.79,100)))\
                                 ,3)
                    else:
                        correlation = loi_de_probalibilite(
                                        correl_seuil-p_correl, 
                                        p_correl,
                                        nbre_correlations, 
                                        correl_seuil)
                        dico_proba_case[(row,col)] = correlation;
                else:
                    # proba entre 0 - 0.5 ( 0 -- X-0.3)  ===> vrai negatifs (distributions asymetrique selon fonction loi_proba)
                    if loi_stats == "poisson":
                        dico_proba_case[(row,col)] = \
                            round(random.choice(
                                    np.union1d(np.linspace(0,0.117,100), \
                                               np.linspace(0.35,0.5,100)))\
                                 ,3)
                    else:
                        correlation = loi_de_probalibilite(
                                        0, 
                                        correl_seuil-p_correl-0.01,
                                        nbre_correlations, 
                                        correl_seuil);
                        dico_proba_case[(row,col)] =  correlation;              
            elif matE.loc[row,col] == 1:
                if (row,col) in dico_deleted_add_edges["ajouter"] or \
                    (col,row) in dico_deleted_add_edges["ajouter"]:
                    # proba entre 0.8 = X ===> faux positifs (distributions uniforme selon fonction loi_proba)
                    correlation = loi_de_probalibilite(
                                    p_correl, 
                                    correl_seuil-0.001,\
                                    nbre_correlations, 
                                    correl_seuil)
                    dico_proba_case[(row,col)] = correlation
                else:
                    # proba entre 0.8001 - 1 ( X -- 1 ) ===> vrai positifs (distributions asymetrique selon fonction loi_proba)
                    correlation = loi_de_probalibilite( 
                                    correl_seuil, 
                                    1,
                                    nbre_correlations, 
                                    correl_seuil)
                    dico_proba_case[(row,col)] = correlation
    return dico_proba_case;               

def ajouter_proba_matE_TROP_COMPLEX(matE, 
                       dico_deleted_add_edges, 
                       loi_stats, 
                       p_correl, 
                       correl_seuil=0.8):
    """
    le but est de definir des probabilites pour chaque case de matE pour simuler 
    la matrice de correlation obtenue apres calcul. les probas (correls) sont les suivantes:
        
        * case 0 -----> 0: proba entre 0 - 0.5 (0, correl_seuil-p_correl-0.01)           ==> vrai negatifs
        * case 1 -----> 0: proba entre 0.6 - 0.79 ( correl_seuil-p_correl, correl_seuil) ==> faux negatifs
        * case 0 -----> 1: proba entre p_correl et correl_seuil-0.001 (p_correl,correl_seuil-0.01)==> faux positifs
        * case 1 -----> 1: proba entre 0.8 - 1 (correl_seuil, 1)             ==> vrai positifs
    NB: X = correl_seuil
    dico_deleted_add_edges = {"ajouter":[(a,b),...],"supprimer":[(c,d),...]}
    loi_stats = loi uniforme, de poisson
    """
    cases_traites = list(); dico_proba_case = dict(); nbre_correlations = 250;
    colonnes = matE.columns.tolist();
    for id_row, row in enumerate(colonnes):
        for col in colonnes[id_row+1:]:
            if col != row and (row,col) not in cases_traites and \
                (row,col) not in cases_traites:
                if col != row and matE.loc[row,col] == 0:
                    if (row,col) in dico_deleted_add_edges["supprimer"] or \
                        (col,row) in dico_deleted_add_edges["supprimer"]:
                        # proba entre 0.6 et 0.79 (entre X-0.2 -- X-0.01) ===> faux negatifs (distributions uniforme)
                        if loi_stats == "poisson":
                            dico_proba_case[(row,col)] = \
                                round(random.choice(
                                        np.union1d(np.linspace(0.6,0.644,100), \
                                                   np.linspace(0.72,0.79,100)))\
                                     ,3)
                        else:
                            correlation = loi_de_probalibilite(
                                            correl_seuil-p_correl, 
                                            p_correl,
                                            nbre_correlations, 
                                            correl_seuil)
                            dico_proba_case[(row,col)] = correlation;
                        cases_traites.append( (row,col) )
                    else:
                        # proba entre 0 - 0.5 ( 0 -- X-0.3)  ===> vrai negatifs (distributions asymetrique selon fonction loi_proba)
                        if loi_stats == "poisson":
                            dico_proba_case[(row,col)] = \
                                round(random.choice(
                                        np.union1d(np.linspace(0,0.117,100), \
                                                   np.linspace(0.35,0.5,100)))\
                                     ,3)
                        else:
                            correlation = loi_de_probalibilite(
                                            0, 
                                            correl_seuil-p_correl-0.01,
                                            nbre_correlations, 
                                            correl_seuil);
                            dico_proba_case[(row,col)] =  correlation;              
                        cases_traites.append( (row,col) )
                elif col != row and matE.loc[row][col] == 1:
                    if (row,col) in dico_deleted_add_edges["ajouter"] or \
                        (col,row) in dico_deleted_add_edges["ajouter"]:
                        # proba entre 0.8 = X ===> faux positifs (distributions uniforme selon fonction loi_proba)
                        correlation = loi_de_probalibilite(
                                        p_correl, 
                                        correl_seuil-0.001,\
                                        nbre_correlations, 
                                        correl_seuil)
                        dico_proba_case[(row,col)] = correlation
                        cases_traites.append( (row,col) );
                    else:
                        # proba entre 0.8001 - 1 ( X -- 1 ) ===> vrai positifs (distributions asymetrique selon fonction loi_proba)
                        correlation = loi_de_probalibilite( 
                                        correl_seuil, 
                                        1,
                                        nbre_correlations, 
                                        correl_seuil)
                        dico_proba_case[(row,col)] = correlation
                        cases_traites.append( (row,col) );
                        
    return dico_proba_case;
###############################################################################
#                ajouter proba de chaque case de matE ====> fin     
###############################################################################
    

"""
application de l'algo de couverture au graphe modifie cad 
    (LG_k: LG modifie de k aretes) modifie 
"""
###############################################################################
#                determiner la voisinage de chaque sommet ====> debut     
###############################################################################
def gamma_noeud(matE, liste_aretes):
    """
    but: determine le nombre de voisin d'un noeud et la liste des noeuds de son voisinage
    return dico
    """
    dico = dict();
    if liste_aretes:
        for noeud in matE.columns:
            ens = set();
            for arc in liste_aretes:
                if noeud == arc[0]:
                    ens.add( arc[1] )
                if noeud == arc[1]:
                    ens.add( arc[0] )
            dico[noeud] = [len(ens), ens]
    else:
        for noeud in matE.columns:
            ens = frozenset([xj for xj in matE.columns \
                         if matE.loc[noeud][xj] == 1]);
            dico[noeud] = [len(ens), ens];
    return dico;

def couverture_par_sommets(sommets_matE, C):
    """ 
    retourne les cliques couvrants tous les sommets d'un graphe. 
    """ 
    dico_sommets_par_cliqs = dict();
    for cliq in C:
        for sommet in cliq:
            if sommet not in dico_sommets_par_cliqs.keys():
                dico_sommets_par_cliqs[sommet] = [cliq];
            else:
                dico_sommets_par_cliqs[sommet].append(cliq);
    sommets_not_in_cliq = set(sommets_matE) - set(dico_sommets_par_cliqs.keys());
    for sommet in sommets_not_in_cliq:
        dico_sommets_par_cliqs[sommet] = [];
    return dico_sommets_par_cliqs;
###############################################################################
#                determiner la voisinage de chaque sommet ====> fin     
###############################################################################
   
    

"""
application de l'algorithme de correction
"""
###############################################################################
#          determiner les aretes de toutes les cliques ====> debut    
###############################################################################
def gamma(liste_arcs, noeud):
    """
    recherche le voisinage de "noeud"
    cad pour chaque arc si noeud est une extremite de cet arc
    """        
    voisins = list();        
    voisins = list(filter(lambda arete: arete[0] == noeud \
                                                or arete[1] == noeud, 
                                liste_arcs));
    return voisins;

def aretes_dans_cliques(C):
    """ 
    retourne les aretes de tous les cliques. 
    """

    f_subset = lambda elt: type(elt) == list \
                            or type(elt) == set \
                            or type(elt) == frozenset;
                            
    aretes_cliques = list();
    boolean_subset = True if list(filter(f_subset, C)) else False;
    
    if boolean_subset:
        aretes_cliques = [item for sublist in [list(it.combinations(c,2)) 
                                            for c in C] 
                        for item in sublist]
    else:
        aretes_cliques = list(it.combinations(C,2));
    return aretes_cliques;

def determiner_aretes_cliques(cliques):
    """
    former un ensemble contenant les aretes des cliques.
    """
    
    f_subset = lambda elt: type(elt) == list \
                            or type(elt) == set \
                            or type(elt) == frozenset;
    aretes_cliques = []
    boolean_subset = False;
    boolean_subset = True if list(filter(f_subset, cliques)) else False;
            
    if boolean_subset:
        for Cu in cliques:
            aretes_cliques.extend( it.combinations(Cu,2) );
    else:
        aretes_cliques = list(it.combinations(cliques,2));
    return aretes_cliques;

def aretes_differente(aretes_Ec, aretes_cible):
    """ retourner le nombre d'aretes differente entre aretes_Ec, aretes_cible. """
    res = set()
    for arete in aretes_cible:
        if (arete[0], arete[1]) not in aretes_Ec and \
            (arete[1], arete[0]) not in aretes_Ec:
            res.add((arete[0], arete[1]))
#    res = aretes_Ec.union(aretes_cible) - aretes_Ec.intersection(aretes_cible)         
    return res;
###############################################################################
#         determiner les aretes de toutes les cliques ====> fin     
###############################################################################