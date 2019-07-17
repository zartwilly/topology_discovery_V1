#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 15:15:20 2019

@author: willy
"""

import fonctions_auxiliaires as fct_aux;
import algorithme_correction_k_1 as alg_corr_k_1;
 import algorithme_correction_k as alg_corr_k_sup_1;

###############################################################################
#                   corriger les sommets a -1 ===> debut
###############################################################################
def correction_cliques(dico_correction, 
                       aretes_matE_k_alpha,
                       dico_gamma_noeud, 
                       dico_parametres_new):
    """
    corriger les sommets  couverts par 3 ou plusieurs cliques
    """
    print("cliques_sans_aretes ={}".format(len(dico_correction["C"])))
    dico_correction["C"] = dico_correction["C"] \
                            + list(map(set, dico_correction['aretes_restantes']));
    print("cliques_avec_aretes ={}".format(len(dico_correction["C"])))
    
    dico_correction["sommets_par_cliqs_avec_aretes"] = \
                fct_aux.couverture_par_sommets(
                    sommets_matE=list(dico_correction["etats_sommets"].keys()),
                    C=dico_correction["C"]);
                        
    if dico_correction["k_erreur"] == 1:
        dico_correction = alg_corr_k_1.correction_cliques_k_1(
                                dico_correction, 
                                dico_gamma_noeud, 
                                dico_parametres_new);
    elif dico_correction["k_erreur"] > 1:
        dico_correction = alg_corr_k_sup_1.correction_cliques_k(
                                dico_correction, 
                                aretes_matE_k_alpha,
                                dico_gamma_noeud, 
                                dico_parametres_new);
                
    return dico_correction;
###############################################################################
#                   corriger les sommets a -1 ===> fin
###############################################################################