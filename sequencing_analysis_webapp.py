from flask import Flask, request, render_template, jsonify, send_file
import warnings
import threading
import time
from pathlib import Path
from queue import Queue
import zipfile
import io
import logging
import json
import numpy as np
from werkzeug.utils import secure_filename
from flask_cors import CORS

# Biopython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner, substitution_matrices
import Bio
import pandas as pd  # for Excel export

warnings.filterwarnings("ignore")

app = Flask(__name__)
CORS(app)

# ============== CONFIGURATION ==============
BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / "outputs"
OUTPUT_DIR.mkdir(exist_ok=True)

ALLOWED_EXTENSIONS = {'ab1', 'fasta', 'fa', 'fas'}
status_queues = {}

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

# ============== BIOLOGICAL CORE LOGIC ==============
def find_best_orf(dna_seq):
    """Find ORF using ATG if present, else longest AA product."""
    dna = dna_seq.upper().replace("-", "")
    best_prot = ""
    best_start = 0

    for frame in range(3):
        seq = dna[frame:]
        codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
        current = []
        started = False
        start_pos = frame

        for i, codon in enumerate(codons):
            aa = Bio.Seq.translate(codon, to_stop=False)

            if codon == "ATG" and not started:
                current = ["M"]
                started = True
                start_pos = frame + i * 3
                continue

            if started:
                if aa == "*":
                    break
                current.append(aa)

        prot = "".join(current)
        if len(prot) > len(best_prot):
            best_prot = prot
            best_start = start_pos

    return best_prot, best_start

def align_proteins(protA, protB):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alns = aligner.align(protA, protB)
    if not alns:
        return None, None
    aln = alns[0]
    return str(aln.query), str(aln.target)

def analyze_alignment(alnA, alnB):
    mutations = []
    frameshifts = []

    for i, (a, b) in enumerate(zip(alnA, alnB), start=1):
        if a == "-" or b == "-":
            frameshifts.append(i)
        elif a != b:
            mutations.append(f"{b}{i}{a}")

    return mutations, frameshifts

# ============== SEQUENCE ANALYSIS ==============
def analyze_single_sequence(record_path, wt_record, rev_identifier, logger):
    try:
        record_path = str(record_path)
        file_type = "abi" if record_path.endswith(".ab1") else "fasta"
        record = SeqIO.read(record_path, file_type)

        if rev_identifier and rev_identifier in record_path:
            record = record.reverse_complement()

        wt_prot, wt_start = find_best_orf(str(wt_record.seq))
        rec_prot, rec_start = find_best_orf(str(record.seq))

        alnA, alnB = align_proteins(rec_prot, wt_prot)
        if alnA is None:
            raise RuntimeError("Protein alignment failed")

        mutations, frameshifts = analyze_alignment(alnA, alnB)

        return {
            "name": Path(record_path).stem,
            "mutations": mutations,
            "mutation_count": len(mutations),
            "frameshifts": frameshifts,
            "frameshift_count": len(frameshifts),
            "sequence_length": len(record.seq),
            "status": "completed",
        }

    except Exception as e:
        logger.error(str(e))
        return {
            "name": Path(record_path).stem,
            "status": "error",
            "error_message": str(e),
            "mutations": [],
            "mutation_count": 0,
            "frameshifts": [],
            "frameshift_count": 0,
            "sequence_length": 0,
        }

# ============== PIPELINE ==============
def run_sequence_analysis(job_dir, template_path, sequence_files, rev_identifier, status_queue, timestamp):
    logger = logging.getLogger(f"job_{timestamp}")
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(job_dir / f"log_{timestamp}.log")
    logger.addHandler(handler)

    wt = SeqIO.read(template_path, "abi" if str(template_path).endswith(".ab1") else "fasta")

    results = {
        "status": "processing",
        "alignmentTracks": [],
        "summary": {"completed": 0, "errors": 0},
    }

    for f in sequence_files:
        status_queue.put(f"Processing {f.name}")
        data = analyze_single_sequence(f, wt, rev_identifier, logger)
        results["alignmentTracks"].append({"sequenceData": data})

        if data["status"] == "completed":
            results["summary"]["completed"] += 1
        else:
            results["summary"]["errors"] += 1

    results["status"] = "completed"

    with open(job_dir / f"results_{timestamp}.json", "w") as fh:
        json.dump(results, fh, indent=2)

    status_queue.put(f"analysis_complete:{timestamp}")

# ============== ROUTES ==============
@app.route('/')
def index():
    """Render the input form."""
    return render_template('index.html')

@app.route('/submit', methods=['POST'])
def submit():
    """Handle form submission and start analysis."""
    timestamp = int(time.time())
    job_dir = OUTPUT_DIR / f"job_{timestamp}"
    
    try:
        job_dir.mkdir(exist_ok=True)
    except Exception as e:
        return jsonify({"error": f"Failed to create job directory: {str(e)}"}), 500

    # Setup logging
    job_log_file = job_dir / f"analysis_job_{timestamp}.log"
    logger = logging.getLogger(f"seq_job_{timestamp}")
    logger.setLevel(logging.DEBUG)
    if not logger.handlers:
        handler = logging.FileHandler(job_log_file)
        handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s'))
        logger.addHandler(handler)

    logger.debug("Received form submission")

    # Get reverse identifier
    rev_identifier = request.form.get('rev_identifier', '').strip()

    # Handle template file
    if 'template_file' not in request.files:
        return jsonify({"error": "Template file is required"}), 400
    
    template_file = request.files['template_file']
    if template_file.filename == '':
        return jsonify({"error": "No template file selected"}), 400
    
    if not allowed_file(template_file.filename):
        return jsonify({"error": "Invalid template file type. Allowed: .ab1, .fasta, .fa, .fas"}), 400

    # Save template file
    template_filename = secure_filename(template_file.filename)
    template_path = job_dir / f"template_{template_filename}"
    template_file.save(template_path)
    logger.info(f"Template saved: {template_path}")

    # Handle sequence files
    if 'sequence_files' not in request.files:
        return jsonify({"error": "At least one sequence file is required"}), 400
    
    sequence_files = request.files.getlist('sequence_files')
    if not sequence_files or all(f.filename == '' for f in sequence_files):
        return jsonify({"error": "No sequence files selected"}), 400

    # Save sequence files
    saved_sequence_files = []
    for seq_file in sequence_files:
        if seq_file.filename and allowed_file(seq_file.filename):
            filename = secure_filename(seq_file.filename)
            file_path = job_dir / filename
            seq_file.save(file_path)
            saved_sequence_files.append(file_path)
            logger.debug(f"Saved sequence file: {file_path}")

    if not saved_sequence_files:
        return jsonify({"error": "No valid sequence files provided"}), 400

    logger.info(f"Saved {len(saved_sequence_files)} sequence files")

    # Create status queue for this job
    status_queues[str(timestamp)] = Queue()

    # Start analysis in background thread
    try:
        threading.Thread(
            target=run_sequence_analysis,
            args=(job_dir, template_path, saved_sequence_files, rev_identifier, 
                  status_queues[str(timestamp)], timestamp),
            daemon=True
        ).start()
    except Exception as e:
        return jsonify({"error": f"Failed to start analysis: {str(e)}"}), 500

    logger.info(f"Analysis started for job {timestamp}")
    return jsonify({
        "status": "Analysis started",
        "timestamp": timestamp,
        "num_files": len(saved_sequence_files)
    })

# --- The rest of the routes (status, results, download_json, download_zip, download_log, download_excel)
# remain unchanged, just remove any references to `basic_auth.required`.

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8010, debug=False)
