from Bio import SeqIO


def read_sequences_from_fastq(fastq_file) -> dict[str, str]:
    """
    Lit un fichier FASTQ et retourne une liste de sequences

    Args:
        fastq_file (str): Le chemin vers le fichier FASTQ.

    Returns:
        list: La liste des sequences lues.
    """
    sequences = {}
    for record in SeqIO.parse(fastq_file, "fastq"):
        sequences[record.id] = str(record.seq)
    return sequences


def read_single_fasta_sequence(file_path: str) -> str:
    """
    Lit une seule sequence d'un fichier FASTA, en gérant les lignes enroulées.
    Args:
        file_path (str): Le chemin vers le fichier FASTA.

    Returns:
        str: La séquence lue.

    Raises:
        ValueError: Si le fichier FASTA est vide ou ne contient aucune séquence.
    """
    sequence = ""
    with open(file_path, "r") as file:
        lines = file.readlines()
        if not lines:
            raise ValueError("Le fichier FASTA est vide.")

        if not lines[0].startswith(">"):
            raise ValueError("Le fichier FASTA ne commence pas par un en-tête (>).")

        for line in lines[1:]:
            sequence += line.strip()

    if not sequence:
        raise ValueError("Le fichier FASTA ne contient aucune séquence.")

    return sequence


def read_fasta_sequences(file_path: str) -> dict[str, str]:
    """
    Lit un fichier FASTA et retourne un dictionnaire avec les identifiants comme clés et les séquences comme valeurs.
    Gère les lignes enroulées.
    Args:
        file_path (str): Le chemin vers le fichier FASTA.
    Returns:
        dict: Un dictionnaire avec les identifiants comme clés et les séquences comme valeurs.
    """
    sequences = {}
    current_id = None
    current_seq = []
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:]  # Enlever le '>'
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            sequences[current_id] = "".join(current_seq)
    return sequences



