# Sophie Huebler — Dissertation Workspace
## Biostatistics PhD | GVHD Microbiome + Bayesian Hierarchical Methods

---

## Dissertation Overview
Three-project dissertation on microbiome analysis in graft-versus-host disease (GVHD),
grounded in the T32 grant proposal (see proposal/). Projects connect: Project 1 builds
the harmonized dataset that Project 3 will apply Project 2's method to.

## Active Projects

### Project 1 — ODSiData (`../ODSi/ODSiData`)
Curation and harmonization of open-source 16S rRNA microbiome datasets from allogeneic
stem cell transplant / GVHD studies. QIIME2-based processing pipeline, R/phyloseq
downstream analysis. GitHub: https://github.com/s-huebler/ODSiData

### Project 2 — bhCRR (`../bhCRR`)
Bayesian Hierarchical Competing Risks Regression — R package implementing a semi-
supervised spike-and-slab lasso method for high-dimensional competing risks data.
GitHub: https://github.com/s-huebler/bhCRR

### Project 3 (planned)
Application of bhCRR to the harmonized ODSiData cohort.

## Related Archives
- `../SSL/CmpRsk_SSL/` — early scratch work that evolved into bhCRR. Archival,
  not actively developed, but useful historical context.

## Key Paths
| Resource              | Path                                      |
|-----------------------|-------------------------------------------|
| Project 1 data        | ~/Documents/ODSi/ODSiData                 |
| Project 2 package     | ~/Documents/bhCRR                         |
| T32 proposal          | ~/Documents/dissertation-hub/proposal/    |
| Shared references     | ~/Documents/dissertation-hub/shared/references/ |

## Manuscript Locations
- Project 1 manuscript: Overleaf [add link] → sync target: ODSiData/outputs/
- Project 2 manuscript: Overleaf [add link] → sync target: bhCRR/outputs/

## Environment Notes
- Machine: MacBook Pro M4 Pro (arm64 / Apple Silicon)
- QIIME2 requires Rosetta x86 emulation — see Project 1 CLAUDE.md for details
- R dependency management: renv (both projects have renv.lock)
- Default conda (arm64): ~/miniconda3
- QIIME conda (x86): ~/miniconda3-x86_64/envs/qiime2-env

## CHPC (University of Utah)
- Connect: ssh [username]@chpc.utah.edu
- Workflow: develop locally → push to GitHub → pull on CHPC → submit SLURM jobs
- Job scripts live in each project's chpc/ folder

## Git Conventions
- Project 1 commits: "[Study] brief description" e.g. "[Liu2017] fix DADA2 params"
- Project 2 commits: "feat/fix/docs/test: brief description"
- Always run devtools::check(vignettes=FALSE) before pushing bhCRR

## Prompt Requests (IMPORTANT)
When Sophie asks for a "prompt", says something like "help me with a prompt",
"write me a prompt", "draft a prompt", or otherwise requests a prompt, this
ALWAYS means a prompt for **Claude Code** running in the **Positron terminal**.

Formatting rules for these prompts:
- Deliver the prompt as plain text in a fenced code block so it can be copied
  and pasted directly into the Claude Code CLI in the Positron terminal.
- Do NOT wrap the prompt in shell commands, quotes, or `claude` invocations —
  just the raw prompt text Sophie will paste at the Claude Code prompt.
- Keep each prompt self-contained: include relevant file paths (absolute or
  repo-relative), the project context (Project 1 ODSiData vs. Project 2 bhCRR),
  and any conventions Claude Code should follow (e.g., commit style,
  devtools::check() before pushing, renv usage).
- Avoid surrounding commentary inside the code block. Any explanation,
  rationale, or usage notes go OUTSIDE the code block.

Complex / multi-step tasks:
- If the request is complex or multi-step, break it into MULTIPLE separate
  prompts, one per code block, in the order Sophie should run them.
- Number them (Prompt 1, Prompt 2, ...) and give each a short heading
  describing its goal.
- Between prompts, note any expected checkpoint (e.g., "review the diff before
  running Prompt 2", "commit after this step", "run tests here") so Sophie
  knows when to pause.
- Prefer small, reviewable steps over one giant prompt, especially for changes
  that touch package internals, tests, or manuscripts.
