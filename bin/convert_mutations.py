import csv
import json
import re
import os
import sys
import argparse

def parse_nucleotide_mutation(mutation):
    """Parse a nucleotide mutation string"""
    match = re.match(r'nuc:([ACGT])(\d+)([ACGT])', mutation)
    if match:
        ref, pos, alt = match.groups()
        return {"ref": ref, "pos": int(pos), "alt": alt}
    return None

def parse_aa_mutation(mutation):
    """Parse an amino acid mutation string"""
    aa_match = re.match(r'aa:([^:]+):([A-Z])(\d+)([A-Z])\(nuc:([ACGT])(\d+)([ACGT])\)', mutation)
    
    if aa_match:
        gff_feature, ref_aa, pos_aa, alt_aa, nuc_ref, nuc_pos, nuc_alt = aa_match.groups()
        gff_feature = f"cds-{gff_feature}"
        pos_aa = int(pos_aa)
        nuc_pos = int(nuc_pos)
        
        pos_in_codon = (nuc_pos - 1) % 3
        
        genetic_code = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        
        aa_to_codons = {}
        for codon, aa in genetic_code.items():
            if aa not in aa_to_codons:
                aa_to_codons[aa] = []
            aa_to_codons[aa].append(codon)
        
        possible_ref_codons = aa_to_codons.get(ref_aa, [])
        
        ref_candidates = [codon for codon in possible_ref_codons if codon[pos_in_codon] == nuc_ref]
        
        possible_alt_codons = aa_to_codons.get(alt_aa, [])
        
        alt_candidates = [codon for codon in possible_alt_codons if codon[pos_in_codon] == nuc_alt]
        
        best_ref_codon = None
        best_alt_codon = None
        
        if ref_candidates and alt_candidates:
            for ref_codon in ref_candidates:
                for alt_codon in alt_candidates:
                    if all(ref_codon[i] == alt_codon[i] for i in range(3) if i != pos_in_codon):
                        best_ref_codon = ref_codon
                        best_alt_codon = alt_codon
                        break
                if best_ref_codon:
                    break
        
        if not best_ref_codon and ref_candidates:
            best_ref_codon = ref_candidates[0]
        if not best_alt_codon and alt_candidates:
            best_alt_codon = alt_candidates[0]
        
        if not best_ref_codon:
            ref_codon_bases = ['N', 'N', 'N']
            ref_codon_bases[pos_in_codon] = nuc_ref
            best_ref_codon = ''.join(ref_codon_bases)
        
        if not best_alt_codon:
            alt_codon_bases = list(best_ref_codon) if best_ref_codon else ['N', 'N', 'N']
            alt_codon_bases[pos_in_codon] = nuc_alt
            best_alt_codon = ''.join(alt_codon_bases)
        
        return {
            "mutation_string": mutation,
            "GFF_FEATURE": gff_feature,
            "ref_codon": best_ref_codon,
            "alt_codon": best_alt_codon,
            "ref_aa": ref_aa,
            "alt_aa": alt_aa,
            "pos_aa": pos_aa,
            "nuc_pos": nuc_pos,
            "nuc_ref": nuc_ref,
            "nuc_alt": nuc_alt
        }
    return None

def convert_gofasta_output(input_csv):
    """Convert GoFASTA mutation outputs to JSON format."""
    base_name = os.path.splitext(input_csv)[0]
    output_json = f"{base_name}.json"
    
    results = []
    mutation_tracker = {}
    
    try:
        with open(input_csv, 'r') as csvfile:
            reader = csv.reader(csvfile)
            header = next(reader)
            
            for row in reader:
                if len(row) < 2 or not row[1].strip():
                    continue
                    
                query, mutations_str = row
                
                match = re.match(r'Consensus_([^_]+)_([^_]+)_', query)
                if not match:
                    continue
                    
                sra_id, region = match.groups()
                
                mutations = mutations_str.split('|')
                
                nuc_mutations = {}
                aa_mutations = []
                
                for mutation in mutations:
                    if mutation.startswith('nuc:'):
                        nuc_mutation = parse_nucleotide_mutation(mutation)
                        if nuc_mutation:
                            key = f"{nuc_mutation['pos']}_{nuc_mutation['ref']}_{nuc_mutation['alt']}"
                            nuc_mutations[key] = {
                                "sra": sra_id,
                                "region": region,
                                "pos": nuc_mutation["pos"],
                                "ref": nuc_mutation["ref"],
                                "alt": nuc_mutation["alt"],
                                "aa_mutations": []
                            }
                            
                    elif mutation.startswith('aa:'):
                        aa_mutation = parse_aa_mutation(mutation)
                        if aa_mutation:
                            aa_mutations.append(aa_mutation)
                
                for aa_mutation in aa_mutations:
                    nuc_pos = aa_mutation["nuc_pos"]
                    nuc_ref = aa_mutation["nuc_ref"]
                    nuc_alt = aa_mutation["nuc_alt"]
                    nuc_key = f"{nuc_pos}_{nuc_ref}_{nuc_alt}"
                    
                    if nuc_key not in nuc_mutations:
                        nuc_mutations[nuc_key] = {
                            "sra": sra_id,
                            "region": region,
                            "pos": nuc_pos,
                            "ref": nuc_ref,
                            "alt": nuc_alt,
                            "aa_mutations": []
                        }
                    
                    aa_entry = {
                        "GFF_FEATURE": aa_mutation["GFF_FEATURE"],
                        "ref_codon": aa_mutation["ref_codon"],
                        "alt_codon": aa_mutation["alt_codon"],
                        "ref_aa": aa_mutation["ref_aa"],
                        "alt_aa": aa_mutation["alt_aa"],
                        "pos_aa": aa_mutation["pos_aa"]
                    }
                    nuc_mutations[nuc_key]["aa_mutations"].append(aa_entry)
                
                for nuc_key, mutation in nuc_mutations.items():
                    full_key = f"{sra_id}_{region}_{mutation['pos']}_{mutation['ref']}_{mutation['alt']}"
                    
                    if full_key in mutation_tracker:
                        continue
                        
                    results.append(mutation)
                    mutation_tracker[full_key] = True
        
        with open(output_json, 'w') as jsonfile:
            json.dump(results, jsonfile, indent=2)
        
        print(f"Converted {len(results)} mutations from {input_csv} to {output_json}")
        return results, output_json
    except Exception as e:
        print(f"Error processing file: {str(e)}", file=sys.stderr)
        raise

def display_examples(results, n=3):
    """Display n example results for verification."""
    examples_with_aa = [r for r in results if r["aa_mutations"]]
    
    print(f"Examples with amino acid mutations:")
    for i, result in enumerate(examples_with_aa[:n]):
        print(json.dumps(result, indent=2))
        if i < min(n, len(examples_with_aa)) - 1:
            print("-" * 40)
    
    if not examples_with_aa:
        print(f"No examples with AA mutations found. Showing first {n} regular results:")
        for i, result in enumerate(results[:n]):
            print(json.dumps(result, indent=2))
            if i < min(n, len(results)) - 1:
                print("-" * 40)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert GoFASTA mutation outputs to JSON format.')
    parser.add_argument('--input', required=True, help='Input CSV file')
    parser.add_argument('--output', required=True, help='Output JSON file')
    args = parser.parse_args()

    try:
        results, output_file = convert_gofasta_output(args.input)
        if output_file != args.output:
            os.rename(output_file, args.output)
        display_examples(results)
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)