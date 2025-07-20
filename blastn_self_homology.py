import os
import glob
import csv
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

def is_overlapping(hsp, selected_hsps):
    """
    Check if the current HSP overlaps completely within any previously selected HSP.
    hsp: tuple with start and end positions at indices 4 and 5
    selected_hsps: list of previously selected HSP tuples
    """
    start_hsp = hsp[4]
    end_hsp = hsp[5]

    if not selected_hsps:
        return False

    for selected_hsp in selected_hsps:
        start_selected = selected_hsp[4]
        end_selected = selected_hsp[5]
        if start_selected <= start_hsp and end_hsp <= end_selected:
            return True
    return False

def process_fasta(fasta_path):
    """
    Process one FASTA file:
    - Create BLAST db,
    - Run self-BLAST,
    - Parse results,
    - Filter and select non-overlapping alignments,
    - Save output CSV.
    """
    base_name = os.path.splitext(os.path.basename(fasta_path))[0]
    temp_fasta = f"{base_name}_temp.fasta"
    blast_db = f"{base_name}_db"
    blast_xml = f"{base_name}_blast_results.xml"
    output_csv = f"filtered_greedy_alignments_{base_name}.csv"

    # Read sequence from input FASTA
    sequence = SeqIO.read(fasta_path, "fasta").seq

    # Write to temporary FASTA file (clean input)
    with open(temp_fasta, "w") as f:
        SeqIO.write([SeqIO.SeqRecord(sequence, id=base_name, description="")], f, "fasta")

    # Create BLAST nucleotide database
    os.system(f"makeblastdb -in {temp_fasta} -dbtype nucl -out {blast_db}")

    # Run BLASTn self-comparison
    blastn_cline = NcbiblastnCommandline(
        query=temp_fasta,
        db=blast_db,
        evalue=0.001,
        outfmt=5,
        out=blast_xml
    )
    stdout, stderr = blastn_cline()

    # Parse BLAST XML and filter alignments
    all_alignments = []
    seen_alignments = set()

    with open(blast_xml) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    alignment_length = len(hsp.sbjct)
                    identity_percent = (hsp.identities / alignment_length) * 100

                    if alignment_length >= 500 and 95 <= identity_percent < 100.0:
                        alignment_key = (hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end)
                        alignment_key_reversed = (hsp.sbjct_start, hsp.sbjct_end, hsp.query_start, hsp.query_end)

                        if alignment_key not in seen_alignments and alignment_key_reversed not in seen_alignments:
                            seen_alignments.add(alignment_key)
                            all_alignments.append((
                                alignment.title,
                                alignment.hit_id,
                                hsp.query,
                                hsp.sbjct,
                                hsp.query_start,
                                hsp.query_end,
                                hsp.sbjct_start,
                                hsp.sbjct_end,
                                alignment_length,
                                identity_percent
                            ))

    # Sort alignments by descending length
    all_alignments.sort(key=lambda x: x[8], reverse=True)

    # Select non-overlapping alignments greedily
    results = []
    selected_hsps = []

    for alignment in all_alignments:
        if not is_overlapping(alignment, selected_hsps):
            selected_hsps.append(alignment)
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

    print(f"{base_name}: Found {len(results)} non-overlapping filtered alignments.")

    # Save results to CSV
    with open(output_csv, mode="w", newline="") as file:
        fieldnames = [
            "Alignment title", "Alignment hit ID", "Query Sequence", "Subject Sequence",
            "Startposition query", "Endposition query", "Startposition subject",
            "Endposition subject", "Alignment Length", "Identity"
        ]
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"{base_name}: Results saved to {output_csv}.")

    # Clean up temporary files if you want (optional)
    os.remove(temp_fasta)
    os.remove(blast_xml)
    os.remove(blast_db + ".nhr")
    os.remove(blast_db + ".nin")
    os.remove(blast_db + ".nsq")

def main():
    # Find all fasta files in current directory
    fasta_files = glob.glob("*.fasta")

    if not fasta_files:
        print("No FASTA files found in current directory.")
        return

    for fasta_file in fasta_files:
        print(f"Processing {fasta_file}...")
        process_fasta(fasta_file)

if __name__ == "__main__":
    main()
