# 🧬 Sequencing Analyzer – Sanger Mutation Detection

A web-based tool for analyzing **Sanger nucleotide sequencing files** and detecting mutations relative to a reference/template sequence.

---

## 🛠️ Technology Stack

* **Frontend:** HTML, TailwindCSS, jQuery
* **Backend:** Python (Flask or FastAPI recommended)
* **Bioinformatics:** Biopython (`pairwise2`, `BLOSUM62`)
* **File Handling:** Supports `.ab1`, `.fasta`, `.fa`, `.fas` nucleotide sequence files
* **Data Export:** JSON, CSV, ZIP, and log files

---

## ⚙️ Requirements

* Python ≥ 3.9
* Biopython
* Flask or FastAPI
* jQuery (included via CDN)
* TailwindCSS (included via CDN)

## 🚀 Usage

1. **Start the backend server:**

```bash
python sequencing_analysis_webapp.py
```

2. **Open the web interface** in your browser:

```
http://localhost:8010
```

3. **Upload files:**

   * **Template/Wildtype nucleotide sequence** (required)
   * **Sequence nucleotide files** to analyze (required)
   * Optional: Reverse complement identifier

4. **Run Analysis:** Click **“Run Analysis”**

   * Progress displayed in real-time
   * Results and mutations table appear when finished

5. **Download results:**

   * JSON, CSV, ZIP, and analysis logs available via download buttons

---

## 📦 Output

* **Results Table:** Summary of samples, mutations, status, and sequence length
* **Downloadable Files:**

  * JSON: Structured analysis data
  * CSV: Tabular format
  * ZIP: All result files
  * Log: Processing messages

---

## ⚠️ Notes

* Ensure backend routes for `/download_json/<timestamp>`, `/download_csv/<timestamp>`, `/download_zip/<timestamp>`, and `/download_log/<timestamp>` are implemented.
* Analysis may take a few seconds to minutes depending on the number of nucleotide sequence files.
