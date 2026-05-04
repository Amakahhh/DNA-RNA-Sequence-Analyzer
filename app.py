"""
TideLab — DNA & RNA Sequence Analyser
Ocean-themed molecular biology pipeline: raw nucleotide sequence → characterised protein.
Variable and function names follow a tidal/oceanographic naming convention throughout.
"""

from flask import Flask, render_template, request, jsonify
import re
import requests
import time
from collections import Counter

tidal_app = Flask(__name__)

# ---------------------------------------------------------------------------
# Biological data — mRNA codon table
# Format: codon → (full_name, three_letter_code, single_letter_code)
# Stop codons use 'STOP' as the full name; they terminate translation.
# ---------------------------------------------------------------------------
DEEP_CODON_REEF = {
    'UUU': ('Phenylalanine', 'Phe', 'F'), 'UUC': ('Phenylalanine', 'Phe', 'F'),
    'UUA': ('Leucine',       'Leu', 'L'), 'UUG': ('Leucine',       'Leu', 'L'),
    'CUU': ('Leucine',       'Leu', 'L'), 'CUC': ('Leucine',       'Leu', 'L'),
    'CUA': ('Leucine',       'Leu', 'L'), 'CUG': ('Leucine',       'Leu', 'L'),
    'AUU': ('Isoleucine',    'Ile', 'I'), 'AUC': ('Isoleucine',    'Ile', 'I'),
    'AUA': ('Isoleucine',    'Ile', 'I'), 'AUG': ('Methionine',    'Met', 'M'),
    'GUU': ('Valine',        'Val', 'V'), 'GUC': ('Valine',        'Val', 'V'),
    'GUA': ('Valine',        'Val', 'V'), 'GUG': ('Valine',        'Val', 'V'),
    'UCU': ('Serine',        'Ser', 'S'), 'UCC': ('Serine',        'Ser', 'S'),
    'UCA': ('Serine',        'Ser', 'S'), 'UCG': ('Serine',        'Ser', 'S'),
    'CCU': ('Proline',       'Pro', 'P'), 'CCC': ('Proline',       'Pro', 'P'),
    'CCA': ('Proline',       'Pro', 'P'), 'CCG': ('Proline',       'Pro', 'P'),
    'ACU': ('Threonine',     'Thr', 'T'), 'ACC': ('Threonine',     'Thr', 'T'),
    'ACA': ('Threonine',     'Thr', 'T'), 'ACG': ('Threonine',     'Thr', 'T'),
    'GCU': ('Alanine',       'Ala', 'A'), 'GCC': ('Alanine',       'Ala', 'A'),
    'GCA': ('Alanine',       'Ala', 'A'), 'GCG': ('Alanine',       'Ala', 'A'),
    'UAU': ('Tyrosine',      'Tyr', 'Y'), 'UAC': ('Tyrosine',      'Tyr', 'Y'),
    'UAA': ('STOP',          'Stop', '*'), 'UAG': ('STOP',          'Stop', '*'),
    'CAU': ('Histidine',     'His', 'H'), 'CAC': ('Histidine',     'His', 'H'),
    'CAA': ('Glutamine',     'Gln', 'Q'), 'CAG': ('Glutamine',     'Gln', 'Q'),
    'AAU': ('Asparagine',    'Asn', 'N'), 'AAC': ('Asparagine',    'Asn', 'N'),
    'AAA': ('Lysine',        'Lys', 'K'), 'AAG': ('Lysine',        'Lys', 'K'),
    'GAU': ('Aspartate',     'Asp', 'D'), 'GAC': ('Aspartate',     'Asp', 'D'),
    'GAA': ('Glutamate',     'Glu', 'E'), 'GAG': ('Glutamate',     'Glu', 'E'),
    'UGU': ('Cysteine',      'Cys', 'C'), 'UGC': ('Cysteine',      'Cys', 'C'),
    'UGA': ('STOP',          'Stop', '*'), 'UGG': ('Tryptophan',   'Trp', 'W'),
    'CGU': ('Arginine',      'Arg', 'R'), 'CGC': ('Arginine',      'Arg', 'R'),
    'CGA': ('Arginine',      'Arg', 'R'), 'CGG': ('Arginine',      'Arg', 'R'),
    'AGU': ('Serine',        'Ser', 'S'), 'AGC': ('Serine',        'Ser', 'S'),
    'AGA': ('Arginine',      'Arg', 'R'), 'AGG': ('Arginine',      'Arg', 'R'),
    'GGU': ('Glycine',       'Gly', 'G'), 'GGC': ('Glycine',       'Gly', 'G'),
    'GGA': ('Glycine',       'Gly', 'G'), 'GGG': ('Glycine',       'Gly', 'G'),
}

# Residue masses in Daltons — used to estimate molecular weight of the protein
REEF_RESIDUE_MASSES = {
    'A': 89.09,  'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16,
    'E': 147.13, 'Q': 146.15, 'G': 75.03,  'H': 155.16, 'I': 131.17,
    'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
    'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15,
}

HYDROPHOBIC_SHORELINE = frozenset({'A', 'V', 'I', 'L', 'M', 'F', 'W', 'P'})
BASIC_TIDEPOOLS       = frozenset({'K', 'R', 'H'})   # positively charged at neutral pH
ACIDIC_TIDEPOOLS      = frozenset({'D', 'E'})         # negatively charged at neutral pH


# ---------------------------------------------------------------------------
# Input preprocessing
# ---------------------------------------------------------------------------

def purge_fasta_headers(raw_tide: str) -> str:
    """Strip FASTA header lines (starting with '>') so only the nucleotide sequence survives."""
    filtered = (
        row.strip()
        for row in raw_tide.strip().splitlines()
        if not row.lstrip().startswith('>')
    )
    return ''.join(filtered).upper()


# ---------------------------------------------------------------------------
# Phase 1 — Sequence classification (DNA vs RNA vs invalid)
# ---------------------------------------------------------------------------

def classify_tide_domain(tide_stream: str) -> dict:
    """
    Identify whether the input is DNA, RNA, or biologically invalid.
    Returns a structured report with detection explanation for the user.
    """
    scrubbed_flow = re.sub(r'[\s\-]', '', tide_stream).upper()

    if not scrubbed_flow:
        return {'valid': False, 'reason': 'The sequence appears to be empty. Please enter at least one nucleotide base.'}

    dna_alphabet   = frozenset('ATCGN')
    rna_alphabet   = frozenset('AUCGN')
    present_bases  = frozenset(scrubbed_flow)
    foreign_bases  = present_bases - (dna_alphabet | rna_alphabet)

    if foreign_bases:
        bad_chars = ', '.join(f'<strong>{c}</strong>' for c in sorted(foreign_bases))
        return {
            'valid': False,
            'reason': (
                f'The sequence contains characters that do not belong to DNA or RNA: {bad_chars}. '
                f'DNA uses only A, T, C, G (and N for ambiguous bases). '
                f'RNA uses only A, U, C, G (and N). Please remove any spaces, numbers, or other characters.'
            )
        }

    has_thymine = 'T' in present_bases
    has_uracil  = 'U' in present_bases

    # Biological impossibility — T and U cannot coexist
    if has_thymine and has_uracil:
        return {
            'valid': False,
            'reason': (
                'The sequence contains both Thymine (T) and Uracil (U), which is biologically impossible. '
                'DNA uses Thymine (T) and RNA uses Uracil (U) — never both at once. '
                'Please check your sequence for errors.'
            )
        }

    if has_thymine:
        ocean_layer    = 'DNA'
        detection_note = (
            f'The system scanned your <strong>{len(scrubbed_flow)}-character</strong> sequence and found '
            f'<strong>Thymine (T)</strong>. Thymine is exclusive to DNA — RNA uses Uracil (U) in its place. '
            f'The presence of T is the definitive indicator that this is a DNA sequence.'
        )
    elif has_uracil:
        ocean_layer    = 'RNA'
        detection_note = (
            f'The system scanned your <strong>{len(scrubbed_flow)}-character</strong> sequence and found '
            f'<strong>Uracil (U)</strong>. Uracil is exclusive to RNA — DNA uses Thymine (T) instead. '
            f'The presence of U confirms this is an RNA sequence.'
        )
    else:
        # Only A, C, G — ambiguous; default to DNA (biological convention)
        ocean_layer    = 'DNA'
        detection_note = (
            f'Your sequence of <strong>{len(scrubbed_flow)} bases</strong> contains only A, C, and G — '
            f'no Thymine (T) or Uracil (U) was found. Since neither distinguishing base is present, '
            f'the system treats it as <strong>DNA</strong> by default. '
            f'If this is RNA, add at least one U base or select RNA manually.'
        )

    base_tally  = dict(sorted(Counter(scrubbed_flow).items()))
    gc_count    = base_tally.get('G', 0) + base_tally.get('C', 0)
    gc_content  = round((gc_count / len(scrubbed_flow)) * 100, 1)

    explanation = (
        '<strong>How sequence type detection works:</strong> The key difference between DNA and RNA '
        'is just one base — DNA uses <em>Thymine (T)</em> while RNA uses <em>Uracil (U)</em>. '
        'All other bases (A, C, G) are shared. So the system simply scans your sequence and checks '
        'which of those two distinguishing bases is present. If both T and U appear, the sequence '
        'is flagged as invalid. If neither is found, DNA is assumed by convention.'
        f'<br><br>{detection_note}'
    )

    return {
        'valid':       True,
        'ocean_layer': ocean_layer,
        'sequence':    scrubbed_flow,
        'length':      len(scrubbed_flow),
        'base_tally':  base_tally,
        'gc_content':  gc_content,
        'explanation': explanation,
    }


# ---------------------------------------------------------------------------
# Phase 2 — Transcription (DNA → mRNA)
# ---------------------------------------------------------------------------

def perform_transcription(scrubbed_flow: str, strand_polarity: str) -> dict:
    """
    Produce the mRNA sequence from a DNA input.
    strand_polarity: 'coding' (non-template/sense) or 'template' (antisense).
    The two strand types require different processing because they run antiparallel.
    """
    if strand_polarity == 'template':
        # Sequences are always entered 5'→3' by convention. The template strand runs antiparallel
        # to the mRNA, so to produce the mRNA (5'→3') we must reverse the template first (giving
        # us its 3'→5' reading direction) and then complement each base, then swap T→U.
        # Equivalently: mRNA = reverse_complement(template).replace('T','U')
        complement_chart = str.maketrans('ATCGN', 'TAGCN')
        surface_ripple   = scrubbed_flow[::-1].translate(complement_chart).replace('T', 'U')
        polarity_note = (
            '<strong>Template strand selected:</strong> The template strand (also called the antisense '
            'or minus strand) is the one that RNA polymerase physically reads — it runs in the '
            '3′ → 5′ direction relative to the mRNA being built. Since sequences are written 5′→3′ '
            'by convention, the system first reversed your sequence (to get the 3′→5′ reading direction), '
            'then complemented each base (A↔T, G↔C), and finally replaced all Thymine (T) with '
            'Uracil (U) to produce the mRNA. This is called the reverse complement operation.'
        )
    else:
        # Non-template (coding/sense) strand already mirrors the mRNA — just swap T→U
        surface_ripple = scrubbed_flow.replace('T', 'U')
        polarity_note = (
            '<strong>Non-template (coding) strand selected:</strong> The coding strand (also called '
            'the sense or plus strand) already has the same sequence as the resulting mRNA — the '
            'only difference is that DNA uses Thymine (T) where RNA uses Uracil (U). '
            'The system simply replaced every T with U to produce the mRNA.'
        )

    explanation = (
        '<strong>What is transcription?</strong> Transcription is the first step of converting a '
        'gene into a protein. Inside the cell\'s nucleus, an enzyme called <em>RNA polymerase</em> '
        'moves along the DNA, reading the template strand, and synthesises a complementary RNA '
        'copy called <em>messenger RNA (mRNA)</em>. The mRNA then travels out of the nucleus '
        'to the ribosome, where it serves as the instruction manual for building a protein.'
        f'<br><br>{polarity_note}'
    )

    return {
        'surface_ripple': surface_ripple,
        'mrna_length':    len(surface_ripple),
        'explanation':    explanation,
    }


# ---------------------------------------------------------------------------
# Phase 3 — Translation (mRNA → amino acid chain)
# ---------------------------------------------------------------------------

def slice_into_wave_triplets(surface_ripple: str) -> list:
    """Partition mRNA into codons (3-base windows). Pads the last codon with '?' if incomplete."""
    return [
        surface_ripple[offset:offset + 3].ljust(3, '?')
        for offset in range(0, len(surface_ripple), 3)
    ]


def perform_translation(surface_ripple: str) -> dict:
    """
    Read the mRNA codon-by-codon and build the polypeptide chain.
    Starts from the first AUG (start codon) and halts at any stop codon.
    Every codon — including upstream UTR codons — is logged for the display table.
    """
    wave_triplets  = slice_into_wave_triplets(surface_ripple)
    codon_log      = []
    pearl_necklace = []

    # Find the first AUG start codon
    anchor_index = next(
        (pos for pos, triplet in enumerate(wave_triplets) if triplet == 'AUG'),
        None
    )
    has_start_signal = anchor_index is not None
    read_from        = anchor_index if has_start_signal else 0

    # Log any upstream (5' UTR) codons before the start site
    for upstream_pos in range(read_from):
        codon_log.append({
            'codon':         wave_triplets[upstream_pos],
            'position':      upstream_pos + 1,
            'amino_acid':    "5' UTR (not translated)",
            'three_letter':  '—',
            'single_letter': '·',
            'status':        'utr',
        })

    chain_halted = False
    for wave_pos in range(read_from, len(wave_triplets)):
        triplet = wave_triplets[wave_pos]

        if '?' in triplet:
            codon_log.append({
                'codon':         triplet,
                'position':      wave_pos + 1,
                'amino_acid':    'Incomplete codon — sequence ended mid-triplet',
                'three_letter':  '—',
                'single_letter': '?',
                'status':        'incomplete',
            })
            break

        reef_entry = DEEP_CODON_REEF.get(triplet)
        if reef_entry is None:
            codon_log.append({
                'codon':         triplet,
                'position':      wave_pos + 1,
                'amino_acid':    'Unrecognised codon',
                'three_letter':  '???',
                'single_letter': '?',
                'status':        'unknown',
            })
            continue

        full_name, three_letter, single_letter = reef_entry

        if full_name == 'STOP':
            codon_log.append({
                'codon':         triplet,
                'position':      wave_pos + 1,
                'amino_acid':    f'STOP codon ({triplet}) — translation ends here',
                'three_letter':  'Stop',
                'single_letter': '*',
                'status':        'stop',
            })
            chain_halted = True
            break

        codon_status = 'start' if (wave_pos == read_from and triplet == 'AUG') else 'coding'
        codon_log.append({
            'codon':         triplet,
            'position':      wave_pos + 1,
            'amino_acid':    full_name,
            'three_letter':  three_letter,
            'single_letter': single_letter,
            'status':        codon_status,
        })
        pearl_necklace.append({'name': full_name, 'three': three_letter, 'one': single_letter})

    reef_sequence = ''.join(bead['one'] for bead in pearl_necklace)

    start_text = (
        'The ribosome located the <strong>start codon AUG</strong> (which always codes for '
        'Methionine) and began building the amino acid chain from there. '
        if has_start_signal else
        'No AUG start codon was found, so translation began at position 1. In a living cell, '
        'this sequence might not produce a protein. '
    )
    stop_text = (
        'Translation continued until a <strong>stop codon</strong> was reached, which signals '
        'the ribosome to release the finished polypeptide chain.'
        if chain_halted else
        'The sequence ended before a stop codon was encountered.'
    )

    explanation = (
        '<strong>What is translation?</strong> Translation is the process where the cell\'s '
        '<em>ribosome</em> reads the mRNA and builds a chain of amino acids. The ribosome moves '
        'along the mRNA three bases at a time, reading each group of three — called a <em>codon</em> '
        '— as a specific instruction. Each codon corresponds to one amino acid (or a start/stop signal). '
        'There are 64 possible codons and 20 amino acids, so most amino acids are encoded by more '
        'than one codon (this redundancy is called degeneracy).'
        f'<br><br>{start_text}{stop_text}'
        f'<br><br>Amino acids produced: <strong>{len(pearl_necklace)}</strong>'
    )

    return {
        'codon_log':       codon_log,
        'pearl_necklace':  pearl_necklace,
        'reef_sequence':   reef_sequence,
        'chain_halted':    chain_halted,
        'has_start_signal': has_start_signal,
        'explanation':     explanation,
    }


# ---------------------------------------------------------------------------
# Phase 4 — Protein characterisation
# ---------------------------------------------------------------------------

def characterise_reef_protein(reef_sequence: str, pearl_necklace: list) -> dict:
    """
    Compute the biochemical composition and key properties of the translated protein.
    Uses the amino acid sequence alone (no 3D structure prediction).
    """
    clean_sequence = reef_sequence.replace('*', '').replace('?', '').replace('-', '').replace('·', '')
    if not clean_sequence:
        return {
            'valid': False,
            'reason': 'No amino acids were produced — unable to characterise a protein from this sequence.',
        }

    residue_tally  = Counter(clean_sequence)
    total_residues = len(clean_sequence)

    # Molecular weight: sum of residue masses minus one water per peptide bond
    WATER_BRIDGE   = 18.02
    estimated_mass = (
        sum(REEF_RESIDUE_MASSES.get(aa, 110.0) * n for aa, n in residue_tally.items())
        - (total_residues - 1) * WATER_BRIDGE
    )

    positive_shore   = sum(residue_tally.get(aa, 0) for aa in BASIC_TIDEPOOLS)
    negative_shore   = sum(residue_tally.get(aa, 0) for aa in ACIDIC_TIDEPOOLS)
    hydrophobic_depth = sum(residue_tally.get(aa, 0) for aa in HYDROPHOBIC_SHORELINE)
    hydrophobic_pct  = round((hydrophobic_depth / total_residues) * 100, 1)

    if positive_shore > negative_shore:
        charge_nature = 'Basic (positively charged at neutral pH)'
    elif negative_shore > positive_shore:
        charge_nature = 'Acidic (negatively charged at neutral pH)'
    else:
        charge_nature = 'Neutral'

    # Build the composition table from the pearl necklace (preserves first-seen order)
    seen = {}
    for bead in pearl_necklace:
        key = bead['one']
        if key not in seen and key not in ('*', '?', '·', '-'):
            seen[key] = {'name': bead['name'], 'three': bead['three'], 'one': bead['one']}
    composition_table = [
        {**entry, 'count': residue_tally.get(entry['one'], 0),
         'pct': round(residue_tally.get(entry['one'], 0) / total_residues * 100, 1)}
        for entry in seen.values()
    ]
    composition_table.sort(key=lambda e: e['count'], reverse=True)

    explanation = (
        '<strong>What is a polypeptide chain, and how does it become a protein?</strong> '
        'After translation, you have a long chain of amino acids linked together by <em>peptide bonds</em> '
        '— this is the <em>polypeptide chain</em>. Amino acids each have different chemical personalities: '
        'some are attracted to water, some repel it, some carry an electrical charge. Because of these '
        'differences, the chain spontaneously folds into a specific three-dimensional shape in the cell\'s '
        'watery environment. That folded shape is the <em>protein</em>. The exact shape — determined '
        'entirely by the amino acid sequence — dictates what job the protein does: carrying oxygen, '
        'fighting infection, catalysing reactions, and so on.'
        '<br><br>'
        f'Your polypeptide has <strong>{total_residues} amino acid residues</strong>. '
        f'Estimated molecular weight: <strong>{estimated_mass:.1f} Da ({estimated_mass / 1000:.2f} kDa)</strong>. '
        f'Overall charge character: <strong>{charge_nature}</strong> '
        f'({positive_shore} basic residues vs {negative_shore} acidic residues). '
        f'<strong>{hydrophobic_pct}%</strong> of residues are hydrophobic, which influences '
        f'how tightly the protein folds and whether it might anchor into a cell membrane.'
    )

    return {
        'valid':             True,
        'clean_sequence':    clean_sequence,
        'total_residues':    total_residues,
        'estimated_mass':    round(estimated_mass, 1),
        'charge_nature':     charge_nature,
        'positive_shore':    positive_shore,
        'negative_shore':    negative_shore,
        'hydrophobic_pct':   hydrophobic_pct,
        'composition_table': composition_table,
        'explanation':       explanation,
    }


# ---------------------------------------------------------------------------
# Phase 5 — External database lookup (EBI BLAST → UniProtKB/Swiss-Prot)
# ---------------------------------------------------------------------------

def dispatch_blast_tide(reef_sequence: str) -> dict:
    """
    Submit a protein sequence to EBI's NCBI BLAST service against Swiss-Prot.
    Returns the job ID for async polling, or an error if submission fails.
    """
    if len(reef_sequence) < 6:
        return {'error': 'Sequence too short for a meaningful database search (minimum 6 amino acids required).'}

    ebi_run_url = 'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run'
    payload = {
        'email':      'analysis@tidelab.app',
        'program':    'blastp',
        'database':   'uniprotkb_swissprot',
        'sequence':   reef_sequence,
        'stype':      'protein',
        'alignments': '5',
        'scores':     '5',
        'exp':        '10',
    }

    try:
        tide_response = requests.post(ebi_run_url, data=payload, timeout=20)
        if tide_response.status_code == 200:
            return {'job_id': tide_response.text.strip(), 'status': 'submitted'}

        # Parse the XML error body for a clearer user message
        body_text = tide_response.text or ''
        if 'blocked' in body_text.lower() or 'fair-use' in body_text.lower():
            user_msg = (
                'The EBI BLAST service has temporarily rate-limited this IP address due to fair-use '
                'policy. This is normal when many test requests are submitted in quick succession. '
                'Please try again in a few minutes — the rest of the analysis above is complete.'
            )
        else:
            user_msg = (
                f'BLAST submission failed (HTTP {tide_response.status_code}). '
                f'The EBI service may be temporarily unavailable. Try again shortly.'
            )
        return {'error': user_msg}
    except requests.exceptions.Timeout:
        return {'error': 'Connection to the BLAST database timed out. Please try again.'}
    except requests.exceptions.RequestException as reef_fault:
        return {'error': f'Network error reaching database: {reef_fault}'}


def retrieve_blast_tide_results(job_id: str) -> dict:
    """
    Poll EBI BLAST for the status of a submitted job and return parsed hits when ready.
    The frontend calls this endpoint on an interval until status is 'complete'.
    """
    status_url = f'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{job_id}'
    result_url = f'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{job_id}/json'

    try:
        status_wave   = requests.get(status_url, timeout=10)
        current_tide  = status_wave.text.strip()
    except requests.exceptions.RequestException:
        return {'status': 'error', 'message': 'Lost connection while checking BLAST job status.'}

    if current_tide in ('RUNNING', 'PENDING', 'QUEUED', 'WAITING'):
        return {'status': 'pending', 'message': 'BLAST is still processing your sequence — checking again shortly…'}

    if current_tide != 'FINISHED':
        return {'status': 'error', 'message': f'BLAST job ended with unexpected status: {current_tide}'}

    # Fetch and parse the results
    try:
        result_wave = requests.get(result_url, timeout=20)
        blast_data  = result_wave.json()
    except Exception as parse_fault:
        return {'status': 'error', 'message': f'Could not parse BLAST results: {parse_fault}'}

    deep_sea_matches = []
    for hit in blast_data.get('hits', [])[:5]:
        # Parse organism from description (UniProt format: "... OS=Organism OX=... ...")
        raw_desc = hit.get('description', 'Unknown protein')
        organism_match = re.search(r'OS=([^=]+?)(?:\s+OX=|$)', raw_desc)
        organism = organism_match.group(1).strip() if organism_match else 'Unknown organism'
        protein_name = re.sub(r'\s+OS=.*', '', raw_desc).strip()

        deep_sea_matches.append({
            'accession':    hit.get('acc', 'N/A'),
            'protein_name': protein_name,
            'organism':     organism,
            'score':        hit.get('score', 0),
            'identity':     hit.get('identity', hit.get('ident', 0)),
            'evalue':       hit.get('expectation', 'N/A'),
        })

    db_explanation = (
        '<strong>What do these database results mean?</strong> '
        'Your protein sequence was submitted to <strong>UniProtKB/Swiss-Prot</strong> — a manually '
        'curated database of over half a million expertly reviewed protein sequences from thousands '
        'of organisms — using <strong>BLAST</strong> (Basic Local Alignment Search Tool). BLAST '
        'compares your sequence against all known proteins and ranks the closest matches. '
        '<br><br>'
        '<em>Score / Bits</em>: Higher values mean a better match. '
        '<em>Identity</em>: The percentage of amino acids that match exactly — 100% would mean '
        'an identical sequence. '
        '<em>E-value</em>: The expected number of matches this close by random chance alone. '
        'Values near zero (like 1e-50) mean the match is almost certainly biologically meaningful, '
        'not a coincidence. If a known protein matches with high identity and a tiny E-value, '
        'your sequence very likely encodes the same protein.'
    )

    return {
        'status':      'complete',
        'matches':     deep_sea_matches,
        'match_count': len(deep_sea_matches),
        'explanation': db_explanation,
    }


# ---------------------------------------------------------------------------
# Flask routes
# ---------------------------------------------------------------------------

@tidal_app.route('/')
def tidal_shore():
    return render_template('index.html')


@tidal_app.route('/api/analyse', methods=['POST'])
def analyse_tide_stream():
    """
    Master analysis endpoint — runs the sequence through all five biological stages
    and returns a fully structured JSON payload for the frontend to render.
    """
    incoming        = request.json or {}
    raw_input       = incoming.get('sequence', '').strip()
    strand_polarity = incoming.get('strand_type', 'coding')

    # Pre-processing: strip FASTA headers
    tide_stream = purge_fasta_headers(raw_input)

    # Stage 1: Classify the molecular domain
    domain_report = classify_tide_domain(tide_stream)
    if not domain_report['valid']:
        return jsonify({'error': domain_report['reason']}), 400

    scrubbed_flow = domain_report['sequence']
    ocean_layer   = domain_report['ocean_layer']
    payload       = {'detection': domain_report}

    # Stage 2: Transcription
    if ocean_layer == 'DNA':
        transcription_data = perform_transcription(scrubbed_flow, strand_polarity)
        surface_ripple     = transcription_data['surface_ripple']
    else:
        surface_ripple     = scrubbed_flow
        transcription_data = {
            'surface_ripple': surface_ripple,
            'mrna_length':    len(surface_ripple),
            'explanation': (
                'Your input was already identified as <strong>RNA</strong>, so transcription '
                '(the DNA → mRNA conversion step) was skipped. Transcription only applies to DNA. '
                'The RNA sequence you provided is used directly as the mRNA for translation.'
            ),
        }
    payload['transcription'] = transcription_data

    # Stage 3: Translation
    translation_data = perform_translation(surface_ripple)
    payload['translation'] = translation_data

    # Stage 4: Protein characterisation
    reef_sequence  = translation_data['reef_sequence']
    pearl_necklace = translation_data['pearl_necklace']
    protein_data   = characterise_reef_protein(reef_sequence, pearl_necklace)
    payload['protein'] = protein_data

    # Stage 5: Submit BLAST job (async — frontend polls for results)
    blast_input = reef_sequence.replace('*', '').replace('?', '')
    if protein_data.get('valid') and len(blast_input) >= 6:
        payload['blast_job'] = dispatch_blast_tide(blast_input)
    else:
        payload['blast_job'] = {'error': 'Sequence is too short for a database search (minimum 6 amino acids).'}

    return jsonify(payload)


@tidal_app.route('/api/blast_result/<job_id>', methods=['GET'])
def fetch_blast_tide_result(job_id: str):
    """Frontend polls this route until BLAST results become available."""
    return jsonify(retrieve_blast_tide_results(job_id))


if __name__ == '__main__':
    tidal_app.run(debug=True, port=5000)
