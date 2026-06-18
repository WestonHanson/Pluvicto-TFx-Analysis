# fragCounter.R
## Running `fragCounter.py` on a Cluster

### 1. Load R module
```bash
ml R/4.1.2-foss-2020b
```

### 2. Check your R library paths
```bash
Rscript -e '.libPaths()'
```

### 3. Create a user-writable R library (if needed)
```bash
mkdir -p ~/R/library
export R_LIBS_USER=~/R/library
```
Verify it’s being used:
```bash
Rscript -e '.libPaths()'
```

### 4. Install required packages
```bash
Rscript -e 'install.packages("devtools", repos="https://cloud.r-project.org")'
Rscript -e 'devtools::install_github("mskilab/fragCounter")'
```

### 5. Verify installation
```bash
Rscript -e 'library(fragCounter)'
```

### 6. Add `frag` to your PATH
```bash
export PATH=${PATH}:$(Rscript -e 'cat(paste0(installed.packages()["fragCounter", "LibPath"], "/fragCounter/extdata/"))')
```

### 7. Test the installation
```bash
frag -h
```

### Notes
- You may need to re-run the `export R_LIBS_USER` and `export PATH` steps in new sessions.
- To make this permanent, add them to your `~/.bashrc` or `~/.zshrc`.
