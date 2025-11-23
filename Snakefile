# Snakefile - pipeline phylogénétique complet
import os
from Bio import SeqIO

# fichier d'entrée (unique)
INPUT_FASTA = "data/seqs-RNA.fasta.txt"

# ---------------------------
# 1. split par espèce (checkpoint pour découvrir dynamiquement les espèces)
# ---------------------------
checkpoint split_by_species:
    input:
        fasta=INPUT_FASTA
    output:
        dir=directory("results/per_species")
    run:
        os.makedirs(output.dir, exist_ok=True)
        seqs = SeqIO.parse(input.fasta, "fasta")
        species_files = {}
        for record in seqs:
            # extrait le nom scientifique après le premier espace et nettoie le point final
            try:
                species_name = record.description.split(' ', 1)[1].strip().replace('.', '')
            except IndexError:
                species_name = record.id
            # safe filename (remplace espaces par underscore)
            species_safe = ''.join(c if c.isalnum() or c=='_' else '_' for c in species_name).rstrip('_')
            if species_safe not in species_files:
                species_files[species_safe] = open(f"{output.dir}/{species_safe}.fasta", "w")
            SeqIO.write(record, species_files[species_safe], "fasta")
        for fh in species_files.values():
            fh.close()

# Fonction pour récupérer dynamiquement les espèces après le checkpoint
def get_species(wildcards):
    ck = checkpoints.split_by_species.get(**wildcards)
    per_species_dir = ck.output[0]
    return [f.rsplit(".",1)[0] for f in os.listdir(per_species_dir) if f.endswith(".fasta")]

# ---------------------------
# 2. QC (contrôle de qualité simple)
# ---------------------------
rule quality_check:
    input:
        fasta="results/per_species/{species}.fasta"
    output:
        qc="results/qc/{species}.qc.txt"
    run:
        os.makedirs("results/qc", exist_ok=True)
        seqs = list(SeqIO.parse(input.fasta, "fasta"))
        n = len(seqs)
        lengths = [len(s.seq) for s in seqs]
        with open(output.qc, "w") as out:
            out.write(f"Species={wildcards.species}\n")
            out.write(f"Sequences={n}\n")
            out.write(f"Min_len={min(lengths) if lengths else 0}\n")
            out.write(f"Max_len={max(lengths) if lengths else 0}\n")
            out.write(f"Mean_len={sum(lengths)/len(lengths) if lengths else 0:.2f}\n")

# ---------------------------
# 3. Annotation (ex: ajouter tag simple)
# ---------------------------
rule annotate:
    input:
        fasta="results/per_species/{species}.fasta",
        qc="results/qc/{species}.qc.txt"
    output:
        annotated="results/annotated/{species}.fasta"
    run:
        os.makedirs("results/annotated", exist_ok=True)
        seqs = list(SeqIO.parse(input.fasta, "fasta"))
        for s in seqs:
            s.description = s.description + " | annotated"
        SeqIO.write(seqs, output.annotated, "fasta")

# ---------------------------
# 4. Alignement MAFFT (par espèce)
# ---------------------------
rule mafft_align:
    input:
        fasta="results/annotated/{species}.fasta"
    output:
        aligned="results/aligned/{species}.aln.fasta"
    benchmark:
        "benchmarks/mafft_{species}.txt"
    log:
        "logs/mafft_{species}.log"
    threads: 1
    run:
        import subprocess
        os.makedirs("results/aligned", exist_ok=True)
        os.makedirs("logs", exist_ok=True)
        with open(output.aligned, "w") as out, open(log[0], "w") as err:
            subprocess.run([
                r"C:\Users\Micro\Downloads\mafft-7.526-win64-signed\mafft-win\mafft.bat", 
                "--auto", 
                input.fasta
            ], stdout=out, stderr=err, check=True)

# ---------------------------
# 5. Construction d'arbre IQ-TREE (par espèce)
# ---------------------------
rule iqtree_build:
    input:
        aligned="results/aligned/{species}.aln.fasta"
    output:
        tree="results/trees/{species}.treefile",
        iqtree="results/trees/{species}.iqtree",
        log="results/trees/{species}.log"
    benchmark:
        "benchmarks/iqtree_{species}.txt"
    log:
        "logs/iqtree_{species}.log"
    threads: 4
    run:
        import subprocess
        from Bio import SeqIO
        
        # Vérifier que l'alignement contient assez de séquences
        seqs = list(SeqIO.parse(input.aligned, "fasta"))
        if len(seqs) < 4:
            # Créer des fichiers vides pour satisfaire Snakemake
            with open(output.tree, "w") as f:
                f.write(f"# Pas assez de séquences pour {wildcards.species} (n={len(seqs)})\n")
            with open(output.iqtree, "w") as f:
                f.write(f"# Skipped: {len(seqs)} sequences\n")
            with open(output.log, "w") as f:
                f.write(f"# Skipped: {len(seqs)} sequences\n")
            with open(log[0], "w") as f:
                f.write(f"Skipped {wildcards.species}: only {len(seqs)} sequences (minimum 4 required)\n")
        else:
            # Exécuter IQ-TREE
            os.makedirs("results/trees", exist_ok=True)
            try:
                subprocess.run([
                    "iqtree3",  # Utilisez iqtree3 au lieu de iqtree2
                    "-s", input.aligned,
                    "-nt", str(threads),
                    "-m", "MFP",
                    "-bb", "1000",
                    "-pre", f"results/trees/{wildcards.species}"
                ], check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                # Sauvegarder l'erreur dans le log
                with open(log[0], "w") as f:
                    f.write(f"STDERR:\n{e.stderr}\n\nSTDOUT:\n{e.stdout}\n")
                raise

# ---------------------------
# 6. Règle finale : tous les arbres (agrégation)
# ---------------------------
rule all:
    input:
        lambda wildcards: expand("results/trees/{species}.treefile", species=get_species(wildcards))
    default_target: True

# ---------------------------
# 7. Règle pour générer un rapport récapitulatif
# ---------------------------
rule report:
    input:
        qc=lambda wildcards: expand("results/qc/{species}.qc.txt", species=get_species(wildcards)),
        trees=lambda wildcards: expand("results/trees/{species}.treefile", species=get_species(wildcards))
    output:
        report="results/pipeline_report.txt"
    run:
        with open(output.report, "w") as f:
            f.write("=== Pipeline Phylogénétique - Rapport ===\n\n")
            f.write(f"Nombre d'espèces analysées: {len(input.qc)}\n\n")
            for qc_file in input.qc:
                f.write(f"--- {os.path.basename(qc_file)} ---\n")
                with open(qc_file) as qc:
                    f.write(qc.read())
                f.write("\n")