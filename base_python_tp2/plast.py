"""
plast.py — Implémentation simplifiée de BLASTN (PLAST) pour le TP2 de bio-informatique.

Ce script cherche une séquence de recherche dans une banque de séquences au format FASTA,
en utilisant une graine (contiguë ou espacée façon PatternHunter), étend les HSP
(High Scoring Pairs) et calcule bitscore et e-value.

Exemple d'utilisation :

    python plast.py -i CGTAGTCGGCTAACGCATACGCTTGATAAGCGTAAGAGCCC \
        -db tRNAs.fasta -E 5 -ss 1e-3 -seed '11111111111'

Options :
    -i    séquence nucléotidique à rechercher (séquence de recherche)
    -db   fichier FASTA contenant la banque de séquences
    -E    seuil de coupure pour l'extension gloutonne des HSP
    -ss   seuil de significativité sur la e-value
    -seed graine utilisée (contiguë '111...1' ou espacée type PatternHunter)
"""

import argparse
import math
import time
from utils import read_fasta_sequences


class HSP:
    """
    Représente un High Scoring Pair (HSP) entre la séquence de recherche
    et une séquence de la banque.

    Attributs:
    - id_sequence_banque (str): Identifiant de la séquence dans la banque.
    - debut_recherche (int): Position de début dans la séquence de recherche.
    - debut_banque (int): Position de début dans la séquence de la banque.
    - longueur (int): Longueur actuelle de l'alignement.
    - score (int): Score brut de l'alignement.
    - bitscore (float | None): Bitscore calculé à partir du score brut.
    - evalue (float | None): E-value associée à l'alignement.
    """

    def __init__(
        self,
        id_sequence_banque: str,
        debut_recherche: int,
        debut_banque: int,
        longueur: int,
    ) -> None:
        self.id_sequence_banque = id_sequence_banque
        self.debut_recherche = debut_recherche
        self.debut_banque = debut_banque
        self.longueur = longueur
        self.score = 0
        self.bitscore = None
        self.evalue = None


MATCH = 5
MISMATCH = -4


def score_pair(base_recherche: str, base_banque: str) -> int:
    """
    Calcule le score entre deux nucléotides.

    Paramètres:
    - base_recherche (str): Nucléotide dans la séquence de recherche.
    - base_banque (str): Nucléotide dans la séquence de la banque.

    Renvoie:
    - int: +5 si les nucléotides sont identiques, -4 sinon.
    """
    if base_recherche == base_banque:
        return MATCH
    return MISMATCH


def score_segment(
    sequence_recherche: str,
    sequence_banque: str,
    position_recherche_debut: int,
    position_banque_debut: int,
    longueur: int,
) -> int:
    """
    Calcule le score total d'un segment aligné de taille donnée.

    Paramètres:
    - sequence_recherche (str): Séquence de recherche complète.
    - sequence_banque (str): Séquence de la banque complète.
    - position_recherche_debut (int): Index de début dans la séquence de recherche.
    - position_banque_debut (int): Index de début dans la séquence de la banque.
    - longueur (int): Longueur du segment aligné.

    Renvoie:
    - int: Score brut du segment aligné.
    """
    score_total = 0
    for offset in range(longueur):
        score_total += score_pair(
            sequence_recherche[position_recherche_debut + offset],
            sequence_banque[position_banque_debut + offset],
        )
    return score_total


def extend_hsp(
    hsp: HSP,
    sequence_recherche: str,
    sequence_banque: str,
    seuil_extension: float,
) -> HSP:
    """
    Étend un HSP dans les deux directions en suivant une heuristique gloutonne.

    Paramètres:
    - hsp (HSP): HSP initial (issu de la graine).
    - sequence_recherche (str): Séquence de recherche.
    - sequence_banque (str): Séquence de la banque correspondante.
    - seuil_extension (float): Seuil d'arrêt (E) basé sur la chute de score.

    Renvoie:
    - HSP: HSP mis à jour avec nouvelles coordonnées, longueur et score.
    """
    # Coordonnées initiales (graine exacte)
    position_recherche_gauche = hsp.debut_recherche
    position_recherche_droite = hsp.debut_recherche + hsp.longueur - 1
    position_banque_gauche = hsp.debut_banque
    position_banque_droite = hsp.debut_banque + hsp.longueur - 1

    # Score initial sur la graine
    score_courant = score_segment(
        sequence_recherche,
        sequence_banque,
        position_recherche_gauche,
        position_banque_gauche,
        hsp.longueur,
    )
    meilleur_score = score_courant
    meilleures_coordonnees = (
        position_recherche_gauche,
        position_recherche_droite,
        position_banque_gauche,
        position_banque_droite,
    )

    # Boucle d'extension gloutonne
    while True:
        meilleur_delta_extension = None
        meilleur_cote_extension = None  # "left" ou "right"

        # Candidat : étendre à gauche
        if position_recherche_gauche > 0 and position_banque_gauche > 0:
            delta_gauche = score_pair(
                sequence_recherche[position_recherche_gauche - 1],
                sequence_banque[position_banque_gauche - 1],
            )
            meilleur_delta_extension = delta_gauche
            meilleur_cote_extension = "left"

        # Candidat : étendre à droite
        if (
            position_recherche_droite + 1 < len(sequence_recherche)
            and position_banque_droite + 1 < len(sequence_banque)
        ):
            delta_droite = score_pair(
                sequence_recherche[position_recherche_droite + 1],
                sequence_banque[position_banque_droite + 1],
            )

            if meilleur_delta_extension is None:
                meilleur_delta_extension = delta_droite
                meilleur_cote_extension = "right"
            elif delta_droite > meilleur_delta_extension:
                meilleur_delta_extension = delta_droite
                meilleur_cote_extension = "right"

        # Si aucune extension possible, on arrête
        if meilleur_cote_extension is None:
            break

        # Appliquer la meilleure extension
        score_courant += meilleur_delta_extension
        if meilleur_cote_extension == "left":
            position_recherche_gauche -= 1
            position_banque_gauche -= 1
        else:
            position_recherche_droite += 1
            position_banque_droite += 1

        # Mise à jour du meilleur score
        if score_courant > meilleur_score:
            meilleur_score = score_courant
            meilleures_coordonnees = (
                position_recherche_gauche,
                position_recherche_droite,
                position_banque_gauche,
                position_banque_droite,
            )

        # Condition d'arrêt basée sur la chute de score
        if meilleur_score - score_courant > seuil_extension:
            break

    (
        position_recherche_gauche,
        position_recherche_droite,
        position_banque_gauche,
        position_banque_droite,
    ) = meilleures_coordonnees

    hsp.debut_recherche = position_recherche_gauche
    hsp.debut_banque = position_banque_gauche
    hsp.longueur = position_recherche_droite - position_recherche_gauche + 1
    hsp.score = meilleur_score

    return hsp


def hsps_overlap(hsp_1: HSP, hsp_2: HSP) -> bool:
    """
    Vérifie si deux HSP se chevauchent à la fois dans la séquence de recherche
    et dans la séquence de la banque.

    Paramètres:
    - hsp_1 (HSP): Premier HSP.
    - hsp_2 (HSP): Deuxième HSP.

    Renvoie:
    - bool: True si les deux HSP se chevauchent dans les deux séquences.
    """
    # Intervalles dans la séquence de recherche
    h1_r_debut = hsp_1.debut_recherche
    h1_r_fin = hsp_1.debut_recherche + hsp_1.longueur - 1
    h2_r_debut = hsp_2.debut_recherche
    h2_r_fin = hsp_2.debut_recherche + hsp_2.longueur - 1

    # Intervalles dans la séquence de la banque
    h1_b_debut = hsp_1.debut_banque
    h1_b_fin = hsp_1.debut_banque + hsp_1.longueur - 1
    h2_b_debut = hsp_2.debut_banque
    h2_b_fin = hsp_2.debut_banque + hsp_2.longueur - 1

    chevauchement_recherche = not (h1_r_fin < h2_r_debut or h2_r_fin < h1_r_debut)
    chevauchement_banque = not (h1_b_fin < h2_b_debut or h2_b_fin < h1_b_debut)

    return chevauchement_recherche and chevauchement_banque


def merge_two_hsps(
    hsp_1: HSP,
    hsp_2: HSP,
    sequence_recherche: str,
    sequence_banque: str,
) -> HSP:
    """
    Fusionne deux HSP chevauchants pour un même sujet.

    Paramètres:
    - hsp_1 (HSP): Premier HSP (sera mis à jour).
    - hsp_2 (HSP): Deuxième HSP.
    - sequence_recherche (str): Séquence de recherche.
    - sequence_banque (str): Séquence de la banque.

    Renvoie:
    - HSP: HSP fusionné (hsp_1 modifié).
    """
    # Union des coordonnées dans la séquence de recherche
    nouvelle_recherche_debut = min(hsp_1.debut_recherche, hsp_2.debut_recherche)
    nouvelle_recherche_fin = max(
        hsp_1.debut_recherche + hsp_1.longueur - 1,
        hsp_2.debut_recherche + hsp_2.longueur - 1,
    )

    # Union des coordonnées dans la séquence de la banque
    nouvelle_banque_debut = min(hsp_1.debut_banque, hsp_2.debut_banque)
    nouvelle_banque_fin = max(
        hsp_1.debut_banque + hsp_1.longueur - 1,
        hsp_2.debut_banque + hsp_2.longueur - 1,
    )

    nouvelle_longueur = nouvelle_recherche_fin - nouvelle_recherche_debut + 1

    # Sécurité : ne pas dépasser les bornes
    max_longueur_recherche = len(sequence_recherche) - nouvelle_recherche_debut
    max_longueur_banque = len(sequence_banque) - nouvelle_banque_debut
    longueur_max_possible = min(max_longueur_recherche, max_longueur_banque)

    if nouvelle_longueur > longueur_max_possible:
        nouvelle_longueur = longueur_max_possible
        nouvelle_recherche_fin = nouvelle_recherche_debut + nouvelle_longueur - 1
        nouvelle_banque_fin = nouvelle_banque_debut + nouvelle_longueur - 1

    nouveau_score = score_segment(
        sequence_recherche,
        sequence_banque,
        nouvelle_recherche_debut,
        nouvelle_banque_debut,
        nouvelle_longueur,
    )

    hsp_1.debut_recherche = nouvelle_recherche_debut
    hsp_1.debut_banque = nouvelle_banque_debut
    hsp_1.longueur = nouvelle_longueur
    hsp_1.score = nouveau_score

    return hsp_1


def same_diagonal(hsp_1: HSP, hsp_2: HSP) -> bool:
    """
    Vérifie si deux HSP sont sur la même diagonale (debut_recherche - debut_banque constant).

    Paramètres:
    - hsp_1 (HSP): Premier HSP.
    - hsp_2 (HSP): Deuxième HSP.

    Renvoie:
    - bool: True si les deux HSP ont la même diagonale.
    """
    return (hsp_1.debut_recherche - hsp_1.debut_banque) == (
        hsp_2.debut_recherche - hsp_2.debut_banque
    )


def merge_overlapping_hsps(
    liste_hsps: list[HSP],
    sequence_recherche: str,
    banque_sequences: dict[str, str],
) -> list[HSP]:
    """
    Fusionne tous les HSP chevauchants qui concernent la même séquence de la banque.

    Paramètres:
    - liste_hsps (list[HSP]): Liste de HSP étendus.
    - sequence_recherche (str): Séquence de recherche.
    - banque_sequences (dict[str, str]): Dictionnaire (id -> séquence) de la banque.

    Renvoie:
    - list[HSP]: Nouvelle liste de HSP après fusion des chevauchements.
    """
    hsps_par_sujet: dict[str, list[HSP]] = {}
    for hsp in liste_hsps:
        hsps_par_sujet.setdefault(hsp.id_sequence_banque, []).append(hsp)

    hsps_fusionnes: list[HSP] = []

    for sujet_id, groupe_hsps in hsps_par_sujet.items():
        if not groupe_hsps:
            continue

        sequence_banque = banque_sequences[sujet_id]

        groupe_trie = sorted(
            groupe_hsps,
            key=lambda h: (h.debut_recherche, h.debut_banque),
        )

        hsp_courant = groupe_trie[0]
        for hsp in groupe_trie[1:]:
            if hsps_overlap(hsp_courant, hsp):
                hsp_courant = merge_two_hsps(
                    hsp_courant,
                    hsp,
                    sequence_recherche,
                    sequence_banque,
                )
            else:
                hsps_fusionnes.append(hsp_courant)
                hsp_courant = hsp

        hsps_fusionnes.append(hsp_courant)

    return hsps_fusionnes


LAMBDA = 0.192
K = 0.176


def compute_bitscore_and_evalue(
    hsp: HSP,
    taille_totale_banque: int,
    taille_sequence_recherche: int,
) -> None:
    """
    Calcule le bitscore et la e-value pour un HSP donné.

    Paramètres:
    - hsp (HSP): HSP dont on veut calculer les statistiques.
    - taille_totale_banque (int): Somme des longueurs des séquences de la banque.
    - taille_sequence_recherche (int): Longueur de la séquence de recherche.

    Renvoie:
    - None: Modifie directement les attributs bitscore et evalue de l'HSP.
    """
    score_brut = hsp.score
    bitscore = round((LAMBDA * score_brut - math.log(K)) / math.log(2))
    hsp.bitscore = bitscore
    hsp.evalue = taille_totale_banque * taille_sequence_recherche * (2 ** (-bitscore))


def select_significant_hsps(
    liste_hsps: list[HSP],
    taille_totale_banque: int,
    taille_sequence_recherche: int,
    seuil_significativite: float,
) -> list[HSP]:
    """
    Sélectionne les HSP significatifs (e-value <= seuil) et garde
    au plus un HSP (le meilleur) par séquence de la banque.

    Paramètres:
    - liste_hsps (list[HSP]): Liste de HSP après fusion.
    - taille_totale_banque (int): Somme des longueurs des séquences de la banque.
    - taille_sequence_recherche (int): Longueur de la séquence de recherche.
    - seuil_significativite (float): Valeur maximale de e-value acceptée.

    Renvoie:
    - list[HSP]: Liste des HSP significatifs (un par séquence de la banque).
    """
    meilleur_hsp_par_sujet: dict[str, HSP] = {}

    for hsp in liste_hsps:
        compute_bitscore_and_evalue(
            hsp,
            taille_totale_banque,
            taille_sequence_recherche,
        )

        if hsp.evalue is None or hsp.evalue > seuil_significativite:
            continue

        hsp_actuel = meilleur_hsp_par_sujet.get(hsp.id_sequence_banque)
        if hsp_actuel is None:
            meilleur_hsp_par_sujet[hsp.id_sequence_banque] = hsp
        else:
            if hsp.bitscore > hsp_actuel.bitscore:
                meilleur_hsp_par_sujet[hsp.id_sequence_banque] = hsp
            elif hsp.bitscore == hsp_actuel.bitscore and hsp.evalue < hsp_actuel.evalue:
                meilleur_hsp_par_sujet[hsp.id_sequence_banque] = hsp

    return list(meilleur_hsp_par_sujet.values())


def fasta_list(chemin_fasta: str) -> list[tuple[str, str]]:
    """
    Charge une banque FASTA et renvoie une liste (id, séquence).

    Paramètres:
    - chemin_fasta (str): Chemin du fichier FASTA.

    Renvoie:
    - list[tuple[str, str]]: Liste de paires (identifiant, séquence).
    """
    banque_sequences = read_fasta_sequences(chemin_fasta)
    return list(banque_sequences.items())


def get_kmers(sequence_recherche: str, graine: str) -> list[tuple[str, int]]:
    """
    Génère tous les motifs issus de la séquence de recherche selon une graine espacée.

    Paramètres:
    - sequence_recherche (str): Séquence de recherche.
    - graine (str): Chaîne de '1' et '0', ex. '111010010100110111'.
      '1' = position utilisée ; '0' = position "don't care".

    Renvoie:
    - list[tuple[str, int]]: Liste de (motif_compressé, position_debut_recherche).
    """
    longueur_fenetre = len(graine)
    positions_ones = [index for index, symbole in enumerate(graine) if symbole == "1"]

    liste_motifs = []
    for position_debut in range(len(sequence_recherche) - longueur_fenetre + 1):
        motif = "".join(
            sequence_recherche[position_debut + offset] for offset in positions_ones
        )
        liste_motifs.append((motif, position_debut))

    return liste_motifs


def find_seed_hits(
    sequence_recherche: str,
    liste_motifs: list[tuple[str, int]],
    banque_liste: list[tuple[str, str]],
    graine: str,
) -> list[HSP]:
    """
    Trouve les HSP initiaux (hits de graine) dans la banque de séquences.

    Paramètres:
    - sequence_recherche (str): Séquence de recherche.
    - liste_motifs (list[tuple[str, int]]): Motifs générés par get_kmers.
    - banque_liste (list[tuple[str, str]]): Liste (id, séquence) de la banque.
    - graine (str): Graine espacée de type PatternHunter ('1'/'0').

    Renvoie:
    - list[HSP]: Liste de HSP initiaux.
    """
    if not liste_motifs:
        return []

    longueur_fenetre = len(graine)
    positions_ones = [index for index, symbole in enumerate(graine) if symbole == "1"]

    index_motifs: dict[str, list[int]] = {}
    for motif, position_recherche in liste_motifs:
        index_motifs.setdefault(motif, []).append(position_recherche)

    liste_hsps: list[HSP] = []

    for id_sequence_banque, sequence_banque in banque_liste:
        longueur_banque = len(sequence_banque)
        for position_banque in range(longueur_banque - longueur_fenetre + 1):
            motif_banque = "".join(
                sequence_banque[position_banque + offset] for offset in positions_ones
            )
            if motif_banque in index_motifs:
                for position_recherche in index_motifs[motif_banque]:
                    liste_hsps.append(
                        HSP(
                            id_sequence_banque,
                            position_recherche,
                            position_banque,
                            longueur_fenetre,
                        )
                    )

    return liste_hsps


def parse_args() -> argparse.Namespace:
    """
    Analyse les arguments de la ligne de commande.

    Renvoie:
    - argparse.Namespace: Objet contenant les attributs i, db, E, ss, seed.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=True,
        help="séquence nucléotidique à rechercher (séquence de recherche)",
    )
    parser.add_argument(
        "-db",
        required=True,
        help="fichier FASTA de la banque",
    )
    parser.add_argument(
        "-E",
        type=float,
        default=4.0,
        help="seuil pour l’extension des HSP",
    )
    parser.add_argument(
        "-ss",
        type=float,
        default=1e-3,
        help="seuil de significativité (e-value)",
    )
    parser.add_argument(
        "-seed",
        default="11111111111",
        help="graine utilisée pour les k-mers (contiguë ou espacée)",
    )

    return parser.parse_args()


def main() -> None:
    """
    Point d'entrée principal du programme.

    Lit les arguments, charge la banque, génère les HSP, les étend,
    les fusionne, calcule les statistiques, filtre par e-value
    puis affiche les meilleurs HSP par séquence de la banque.
    """
    debut = time.perf_counter()
    
    args = parse_args()

    sequence_recherche = args.i

    banque_sequences = read_fasta_sequences(args.db)
    banque_liste = list(banque_sequences.items())

    liste_motifs = get_kmers(sequence_recherche, args.seed)

    hsps_initiaux = find_seed_hits(
        sequence_recherche,
        liste_motifs,
        banque_liste,
        args.seed,
    )

    hsps_etendus: list[HSP] = []
    for hsp in hsps_initiaux:
        sequence_banque = banque_sequences[hsp.id_sequence_banque]
        hsps_etendus.append(
            extend_hsp(
                hsp,
                sequence_recherche,
                sequence_banque,
                args.E,
            )
        )

    hsps_fusionnes = merge_overlapping_hsps(
        hsps_etendus,
        sequence_recherche,
        banque_sequences,
    )

    taille_totale_banque = sum(len(seq) for seq in banque_sequences.values())
    taille_sequence_recherche = len(sequence_recherche)

    hsps_significatifs = select_significant_hsps(
        hsps_fusionnes,
        taille_totale_banque,
        taille_sequence_recherche,
        args.ss,
    )

    if not hsps_significatifs:
        print(f"Aucun HSP significatif (e <= {args.ss}).")
        return

    hsps_significatifs.sort(key=lambda hsp: hsp.evalue)

    nombre_sequences = len(hsps_significatifs)

    for hsp in hsps_significatifs:
        sequence_banque = banque_sequences[hsp.id_sequence_banque]

        debut_recherche = hsp.debut_recherche
        fin_recherche = hsp.debut_recherche + hsp.longueur
        debut_banque = hsp.debut_banque
        fin_banque = hsp.debut_banque + hsp.longueur

        print(f">{hsp.id_sequence_banque}")
        print(
            f"# Best HSP score:{hsp.score:.2f}, "
            f"bitscore:{hsp.bitscore:.2f}, "
            f"evalue: {hsp.evalue:.2e}"
        )
        print(
            f"{debut_recherche} "
            f"{sequence_recherche[debut_recherche:fin_recherche]} "
            f"{fin_recherche}"
        )
        print(
            f"{debut_banque} "
            f"{sequence_banque[debut_banque:fin_banque]} "
            f"{fin_banque}"
        )
        print("----------------------------------------")

    print(f"Total : {nombre_sequences}")

    fin = time.perf_counter()
    print(f"Temps d'exécution : {fin - debut:.4f} secondes")


if __name__ == "__main__":
    main()