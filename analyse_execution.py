#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 11:29:15 2019

@author: willy
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 11:00:50 2019

@author: willy
"""

import os, time;
import numpy as np;
import pandas as pd;
import itertools as it;
import multiprocessing as mp;

import simulation_graphe as simi_graph_k;

from pathlib import Path    # http://stackoverflow.com/questions/6004073/how-can-i-create-directories-recursively

from multiprocessing import Pool;

import smtplib;
import os.path as op;
from email.mime.multipart import MIMEMultipart;
from email.mime.text import MIMEText;
from email.mime.base import MIMEBase;
from email.utils import COMMASPACE, formatdate;
from email import encoders;

import fonctions_auxiliaires as fct_aux;


###############################################################################
#                    envoie mail et envoie mail avec piece jointe
#                               ===> debut
###############################################################################
MY_ADDRESS = "zartwilly@gmail.com"
PASSWORD = "willis38"

def envoie_mail(message, subject):
    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.starttls()
    server.login(MY_ADDRESS, PASSWORD)
    
    msg = MIMEMultipart();
    msg["From"] = MY_ADDRESS;
    msg["To"] = MY_ADDRESS;
    msg["Subject"] = subject;
    
    msg.attach(MIMEText(message, 'plain'))
        
    server.send_message(msg)
    server.quit()
    
    del msg;
    
def envoie_mail_with_pieces_jointes(message, subject, path_files):
    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.starttls()
    server.login(MY_ADDRESS, PASSWORD)
    
    msg = MIMEMultipart();
    msg["From"] = MY_ADDRESS;
    msg["To"] = MY_ADDRESS;
    msg["Subject"] = subject;
    
    msg.attach(MIMEText(message, 'plain'))
    
    files = os.listdir(path_files)
    for path in files:
        part = MIMEBase('application', "octet-stream")
        with open(path_files+"/"+path, 'rb') as file:
            part.set_payload(file.read())
        encoders.encode_base64(part)
        part.add_header('Content-Disposition',
                        'attachment; filename="{}"'.format(op.basename(path)))
        msg.attach(part)
    server.send_message(msg)
    server.quit()
    
    del msg;
###############################################################################
#                    envoie mail et envoie mail avec piece jointe
#                               ===> fin 
###############################################################################
    
###############################################################################
#                   load graphes, modified k edges and 
#                       define parametres ===> debut
###############################################################################
def define_parametres(dico_parametres, dbg_mesures):
    """
    charger les graphes, modifier k aretes et 
    definir les parametres pour chaque execution
    """
    ## initialisation valeurs flots des graphes
    nbre_ts = 10; effet_joule = 0#0.1    
    
    if dico_parametres["DBG"] :
            print("modes={}, p_correls={}, k_erreurs={}, nbre_graphes={}"\
                  .format(len(dico_parametres["mode_select_noeuds_1s"]), 
                          len(dico_parametres["p_correls"]), 
                          len(dico_parametres["k_erreurs"]), 
                          len(dico_parametres["nbre_graphes"]) ))
            
    graphes_GR_LG = list();
    
    for (type_operat, mode, p_correl, k_erreur, nbre_graphe) in it.product(
                    dico_parametres["type_operations"],
                    dico_parametres["mode_select_noeuds_1s"],
                    dico_parametres["p_correls"],
                    dico_parametres["k_erreurs"],
                    dico_parametres["nbre_graphes"]
                    ):
        rep_base = dico_parametres["rep"] \
                    + "/" \
                    + type_operat \
                    + "/" \
                    + mode + "_sommets_GR_" + str(dico_parametres["nbre_sommets_GR"]) \
                    + "/" \
                    + "data_p_" + str(p_correl) \
                    + "/" \
                    + "G_" + str(nbre_graphe) + "_" + str(k_erreur) ;
        #print("rep_base={}".format(rep_base)) 
        
        chemin_matrices = rep_base + "/" + "matrices" + "/";
        chemin_datasets = rep_base + "/" + "datasets" + "/";
        path = Path(chemin_matrices); path.mkdir(parents=True, exist_ok=True);
        path = Path(chemin_datasets); path.mkdir(parents=True, exist_ok=True);
        
        # generer matE
        if not dbg_mesures:
            matE_LG, matA_GR, dico_sommet_arete = \
                        fct_aux.creer_GR_LG(
                                dim_mat = dico_parametres['nbre_sommets_GR'], 
                                nbre_lien = dico_parametres['nbre_lien'], 
                                chemin_matrices = chemin_matrices, 
                                chemin_datasets = chemin_datasets,
                                nbre_ts = nbre_ts,
                                effet_joule = effet_joule,
                                without_mesures = False
                                            )
        else:
            matE_LG, matA_GR, dico_sommet_arete = \
                        fct_aux.creer_GR_LG(
                                dim_mat = dico_parametres['nbre_sommets_GR'], 
                                nbre_lien = dico_parametres['nbre_lien'], 
                                chemin_matrices = chemin_matrices, 
                                chemin_datasets = chemin_datasets,
                                nbre_ts = nbre_ts,
                                effet_joule = effet_joule,
                                without_mesures = True
                                            )
            
                                
        num_graph = "G_" + str(nbre_graphe) + "_" + str(k_erreur) + "_p_" +\
                    "".join(str(p_correl).split('.'))
#        num_graph = "G_" + str(nbre_graphe) + "_" + str(k_erreur)
                    
        dico_parametres_new = dico_parametres.copy();
        dico_parametres_new["coef_fct_cout"] = (dico_parametres_new['exposant'],
                           dico_parametres_new['facteur_multiplicatif'],
                           type_operat);
        dico_parametres_new['p_correl'] = p_correl;
        dico_parametres_new['mode_select_noeuds_1'] = mode;
        graphes_GR_LG.append(
                        (matE_LG, 
                         matA_GR, 
                         dico_sommet_arete,
                         chemin_matrices,
                         chemin_datasets,
                         mode, 
                         p_correl, 
                         k_erreur, 
                         num_graph,
                         rep_base,
                         dico_parametres_new
                         )
                    )
    
    return graphes_GR_LG;
###############################################################################
#                   load graphes, modified k edges and 
#                       define parametres ===> fin
###############################################################################

###############################################################################
#                   execution parallele avec envoie mail + piece jointe
#                               ===> debut
###############################################################################
def execution_parallele_with_mail(graphes_GR_LG, parametres, DBG_PARALLELE):
    """
    """
    list_returnss = list();
    if DBG_PARALLELE:
        p = Pool(mp.cpu_count()-1) 
        list_returnss = p.starmap(
                            simi_graph_k.simulation_p_correl_k, 
                            graphes_GR_LG)
        p.terminate()
    else:
        for graphe_GR_LG in graphes_GR_LG:
            list_returns = simi_graph_k.simulation_p_correl_k(*graphe_GR_LG)
            list_returnss.append(list_returns);
            
    flat_returns = list(it.chain.from_iterable(list_returnss))
    print('flat_returns ={}'.format(len(flat_returns)))
    
    # mettre list_dicos dans un dataframe
    dico_results = dict()
    for dico in flat_returns:
        dico_results[dico['nom_graphe']] = dico;
        
    df_results = pd.DataFrame.from_dict(data=dico_results, orient='index')
    df_results.to_pickle(parametres["rep"]\
                      +"/"+"resume_execution_k_"\
                      +"_".join(map(str, parametres["k_erreurs"]))\
                      +"_nbre_sommets_GR_"\
                      +str(parametres["nbre_sommets_GR"])+".csv") 
    
    # mail 
    k_erreurs_join = ", ".join(map(str, parametres["k_erreurs"]))
    modes = ", ".join(parametres["mode_select_noeuds_1s"]);
    criteres = ", ".join(parametres["type_operations"])
    message = "caracteristique execution \n " \
        + "nombre de sommets de Graphe racine = {} \n".format(
                parametres["nbre_sommets_GR"])\
        + "nombre de graphes = {} \n".format( len(parametres["nbre_graphes"]) )  \
        + "k_erreurs = [{}] \n".format(k_erreurs_join) \
        + "modes correction = {} \n".format(modes)  \
        + "critere correction = {} \n \n \n".format(criteres)  \
        + " execution Termine !!!! "
    subject = "Execution discovery Topology article correction k = 1."
    
    path_files = parametres["rep"] \
                + "/" \
                + parametres["type_operations"][0] \
                + "/" \
                + parametres["mode_select_noeuds_1s"][0] + "_sommets_GR_" \
                + str(parametres["nbre_sommets_GR"]) \
                + "/" \
                + "data_p_" + str(parametres["p_correls"][0]) + "/" + "distribution";
    envoie_mail_with_pieces_jointes(message, subject, path_files)
    
###############################################################################
#                   execution parallele avec envoie mail + piece jointe
#                               ===> fin
###############################################################################

###############################################################################
#           definir les parametres generales d'execution  ===> debut
###############################################################################
def parametres_generales(nbre_sommets_GR, k_erreurs, alpha_max, NBRE_GRAPHE,
                         p_correl, mode_select_noeuds_1, rep):
    """
    definir les parametres generales d'execution de la simulation 
    de nos algorithmes.
    
    p_correl = 0; # valeur de correlation ([0,1]) a partir duquel 
    # *  correl = 0 est transforme en 1 => AJOUT D'ARETES (p_correl = 1)
    # * et correl = 1 est transforme a 0 => SUPPRESSION D'ARETES (p_correl = 0)
                           
    """
#    rep = "../data_repeat";
    nbre_lien = 5;
    methode_delete_add_edges = 0; # 1:selection aretes a delete or a ajouter par liste, 0: par proba
    correl_seuil = 0.7 # seuil a partir duquel on a des correlations fausses positives (0->1) et fausses negatives (1->0)
                        # [ signification des correlations fausses {positives,negatives} a partir de MATE]
    loi_stats = "uniforme";#"poisson";
    SEUIL_PROBA = 0.8; # je pense cest pour ajouter des probas a chaque arete # A EFFACER
    algoGreedy = False #True; # False correction tous les noeuds a -1, True: algo Greedy   
    biais = False;
    critere_selection_pi1_pi2 = 0; # 0: moins de modif,1: ajout aretes> supp aretes, 2:ajout aretes < supp aretes,
    number_permutations_nodes_1= 10 #100;
    number_items_pi1_pi2 = 1;
    
    ascendant_1 = True;
    simulation = True;
    seuil_U = 0;
    epsilon = 0.75; 
    
    type_operations = ["lineaire_simul50Graphes_priorite_aucune"]               # "lineaire","cloche", "lineaire_simul50Graphes_priorite_supp", "lineaire_simul50Graphes_priorite_ajout", "lineaire_simul50Graphes_priorite_aucune"
    p_correls = [p_correl];
    correl_seuils = [0.6];                                                     # correl_seuils = [0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9]; 
    mode_select_noeuds_1s = [mode_select_noeuds_1];                                     # ["degreMin","coutMin","aleatoire"]; 
    facteur_multiplicatif = 1; exposant = 0;                                   # 0: fct_unitaire, 1:fct_normal, 2: fct_quadratique, 4:fct_quadruple, 5: fct_quintuple
    type_fct_cout = "lineaire"                                                 # "cloche" ou "lineaire"
    coef_fct_cout = (exposant, facteur_multiplicatif, type_fct_cout)
    DBG = True;
    
    nbre_graphes = range(1, NBRE_GRAPHE+1, 1)
    
    parametres = {
            "rep": rep, 
            "biais": biais, 
            "nbre_lien": nbre_lien, 
            "ascendant_1": ascendant_1,
            "simulation": simulation,
            "seuil_U": seuil_U,
            "epsilon": epsilon,
            "number_items_pi1_pi2": number_items_pi1_pi2, 
            "number_permutations_nodes_1": number_permutations_nodes_1, 
            "methode_delete_add_edges": methode_delete_add_edges, 
            "loi_stats": loi_stats, 
            "p_correls": p_correls, 
            "correl_seuil": correl_seuil,
            "correl_seuils": correl_seuils,
            "algoGreedy": algoGreedy, 
            "exposant": exposant,
            "facteur_multiplicatif": facteur_multiplicatif,
#            "coef_fct_cout": coef_fct_cout, 
            "type_operations": type_operations, 
            "mode_select_noeuds_1s": mode_select_noeuds_1s,
            "critere_selection_pi1_pi2": critere_selection_pi1_pi2,
            "nbre_sommets_GR": nbre_sommets_GR,
            "k_erreurs": k_erreurs,
            "nbre_graphes": nbre_graphes,
            "alpha_max": alpha_max,
            "DBG": DBG
            }
    return parametres;
###############################################################################
#           definir les parametres generales d'execution  ===> fin
###############################################################################
    
###############################################################################
#           analyser le resume d'execution  ===> debut
###############################################################################
def create_dico_by_sommet_cliques(cliques):
    dico_sommet_mat_GR = dict();
    for clique in cliques:
        aretes = list(map(lambda item: set(item.split("_")),clique))
        set_inter = set.intersection(*aretes)
        if len(set_inter) == 1:
            sommet = set_inter.pop()
            if sommet not in dico_sommet_mat_GR:
                dico_sommet_mat_GR[sommet] = [clique]
            else:
                dico_sommet_mat_GR[sommet].append(clique)
    return dico_sommet_mat_GR;

def create_dico_by_sommet_number(cliques):
    dico_sommet_mat_GR = dict();
    for clique in cliques:
        aretes = list(map(lambda item: set(item.split("_")),clique))
        set_inter = set.intersection(*aretes)
        if len(set_inter) == 1:
            sommet = set_inter.pop()
            if sommet not in dico_sommet_mat_GR:
                dico_sommet_mat_GR[sommet] = 1;
            else:
                dico_sommet_mat_GR[sommet] += 1;
    return dico_sommet_mat_GR;

def create_sommet_set_cliques(dico_by_sommet):
    dico_sommet_set_cliques = dict();
    for sommet, cliques_item in dico_by_sommet.items():
        dico_sommet_set_cliques[sommet] = set(
                                        it.chain.from_iterable(cliques_item))
    return dico_sommet_set_cliques;

def analyse_df_res(parametres):
    """
    analyse du resume de l'execution des graphes
    """
    nom_resume_execution = "resume_execution_k_{}_".format(parametres['k_erreurs'][0])\
                            +"nbre_sommets_GR_{}.csv".format(parametres['nbre_sommets_GR'])
#    rep_base = parametres["rep"] #+ "/" + parametres["type_operations"][0] 
    rep_base = parametres["rep"] \
            + "/" \
            + parametres["type_operations"][0] \
            + "/" \
            + parametres["mode_select_noeuds_1s"]+ "_sommets_GR_" + str(parametres["nbre_sommets_GR"]) \
            + "/" \
            + "data_p_" + str(parametres["p_correls"][0]);
    file = rep_base + "/" + nom_resume_execution;
    
    f = lambda row:"_".join([row.split("_")[0],row.split("_")[1]])

    df_res = pd.read_pickle(file);

    df_res['alpha'] = df_res['alpha'].astype(int);
    df_res["dico_by_sommet"] = df_res["cliques"].apply(
                                        create_dico_by_sommet_cliques)
    df_res["dico_by_sommet_number"] = df_res["cliques"].apply(
                                        create_dico_by_sommet_number)
    df_res["dico_sommet_set_cliques"] = df_res["dico_by_sommet"].apply(
                                        create_sommet_set_cliques)
    
    return df_res;
###############################################################################
#           analyser le resume d'execution  ===> fin
###############################################################################

    