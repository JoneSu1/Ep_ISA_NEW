$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath 'F:\phd\Drophila\3Model_motif_discovering\ISA\Ep_ISA'

Write-Host '=== Step 1: git init ==='
if (-not (Test-Path '.git')) {
    git init
} else {
    Write-Host '(repo already exists, skipping init)'
}

Write-Host "`n=== Step 2: git add + commit ==="
git add -A
$status = git status --porcelain
if ($status) {
    git commit -m "Initial release: Ep_ISA v0.1.0 - TF cooperativity via Fi-NeMo + ISA

- adapter/finemo_io.py: Fi-NeMo hits.tsv -> motif_locs.csv converter
- quickstart.py: EpQuickStart orchestrator (TF/Keras, 4-stage ISA)
- modeling/: TF/Keras predict.py (auto-detect input/output format)
- scoring/: ISA pipeline (single_isa, combi_isa, null, aggregate)
- plotting/ + exploring: 17 visualization functions (dm3 adapted)
- notebooks/Ep_ISA_Analysis.ipynb: full Colab analysis pipeline
- DATA_REQUIREMENTS.md: dm3 reference data setup guide

Author: Junhua"
    Write-Host 'Commit created.'
} else {
    Write-Host '(nothing to commit)'
}

Write-Host "`n=== Step 3: git remote + push ==="
git remote add origin https://github.com/JoneSu1/EP_ISA.git
git branch -M main
git push -u origin main

Write-Host "`n=== Done! ==="
Write-Host 'Repo: https://github.com/JoneSu1/EP_ISA'
