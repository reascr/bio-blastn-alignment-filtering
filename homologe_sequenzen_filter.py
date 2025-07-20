from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import os
import csv

# Lade die FASTA-Dateien (hier ändern, müssen im selben directory liegen)
fasta_file1 = "NSCV1.fasta"

# Extrahiere die Sequenzen aus den FASTA-Dateien
nscv1 = SeqIO.read(fasta_file1, "fasta").seq

# Speichern der Sequenzen in temporären Dateien
with open("nscv1.fasta", "w") as f:
    SeqIO.write([SeqIO.SeqRecord(nscv1, id="nscv1", description="")], f, "fasta")

# Erstelle eine BLAST-Datenbank
os.system("makeblastdb -in nscv1.fasta -dbtype nucl -out nscv1_db")

# Führe BLASTn durch, um seq1 gegen die seq1-Datenbank zu vergleichen. Hier ändern
blastn_cline = NcbiblastnCommandline(query="nscv1.fasta", db="nscv1_db", evalue=0.001, outfmt=5,
                                     out="blast_results.xml")
stdout, stderr = blastn_cline()

# Name hier irrenführend! ändern
def is_overlapping(hsp, selected_hsps):
    '''Prüft ob eine Sequenz innerhalb einer anderen liegt.'''
    start_hsp = hsp[4]
    end_hsp = hsp[5]

    if selected_hsps == []:
        return False

    for selected_hsp in selected_hsps:
        start_selected = selected_hsp[4]
        end_selected = selected_hsp[5]

        # Liegt der Bereich von hsp vollständig innerhalb des Bereichs von selected_hsp?
        if start_selected <= start_hsp and end_hsp <= end_selected:
            return True

    return False



# List to store all filtered but unsorted alignments
all_alignments = []
seen_alignments = set()
# First pass – collect all non-duplicate matches
with open("blast_results.xml") as result_handle:
    blast_records = NCBIXML.parse(result_handle)

    for record in blast_records:
        for alignment in record.alignments:
            hsps = alignment.hsps
            for hsp in hsps:
                # Calculate alignment length and identity percent
                alignment_length = len(hsp.sbjct)
                identity_percent = (hsp.identities / alignment_length) * 100

                # Filter based on minimum length and identity percent
                if alignment_length >= 500 and identity_percent >= 95 and identity_percent < 100.0:
                    alignment_key = tuple([hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end])
                    alignment_key_reversed = tuple([hsp.sbjct_start, hsp.sbjct_end, hsp.query_start, hsp.query_end])
                    if alignment_key and alignment_key_reversed not in seen_alignments: # filter for duplicates (we keep the first seq)
                        seen_alignments.add(alignment_key)

                        alignment_tuple = (
                            alignment.title,  # Alignment title
                            alignment.hit_id,  # Alignment hit ID
                            hsp.query,  # Query sequence
                            hsp.sbjct,  # Subject sequence
                            hsp.query_start,  # Query start
                            hsp.query_end,  # Query end
                            hsp.sbjct_start,  # Subject start
                            hsp.sbjct_end,  # Subject end
                            alignment_length,  # Alignment length
                            identity_percent  # Identity percentage
                        )
                        all_alignments.append(alignment_tuple)

# Second pass – sort all alignments by length in descending order
all_alignments.sort(key=lambda x: x[8], reverse=True)  # Sort by alignment length (index 8)
print(len(all_alignments)) # 314, 262

# List to store the final non-overlapping alignments
results = []

# Third pass – Check for overlaps and greedily select the longest alignments
selected_hsps = [] # the non overlapping longest sequences


for alignment in all_alignments:
    # Check if the alignment does not overlap with any previously selected alignments
    if not is_overlapping(alignment, selected_hsps):
        # If it doesn't overlap, add to the selected list and results
        selected_hsps.append(alignment)

        # Store the alignment in the results list
        results.append({
            "Alignment title": alignment[0],
            "Alignment hit ID": alignment[1],
            "Query Sequence": alignment[2],
            "Subject Sequence": alignment[3],
            "Startposition query": alignment[4],
            "Endposition query": alignment[5],
            "Startposition subject": alignment[6],
            "Endposition subject": alignment[7],
            "Alignment Length": alignment[8],
            "Identity": alignment[9]
        })

# Output the final results
print(f"{len(results)} non-overlapping, filtered alignments found and stored.")

# Speichern der Ergebnisse in eine CSV-Datei. Hier gewünschten Dateinamen ändern.
with open("filtered_greedy_alignments_1.csv", mode="w", newline="") as file:
    writer = csv.DictWriter(file, fieldnames=["Alignment title", "Alignment hit ID", "Query Sequence", "Subject Sequence",
                                              "Startposition query", "Endposition query", "Startposition subject",
                                              "Endposition subject", "Alignment Length", "Identity"])
    writer.writeheader()
    for row in results:
        writer.writerow(row)

print(f"{len(results)} gefilterte Ergebnisse wurden gespeichert.")
