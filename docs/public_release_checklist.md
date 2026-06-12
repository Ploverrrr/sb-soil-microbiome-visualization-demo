# Public Release Checklist

- [x] No private raw data included
- [x] No absolute local paths included
- [x] No manuscript-only files included
- [x] No tracked `.DS_Store` files
- [x] No tracked `__MACOSX` folders
- [x] No large sequencing files
- [x] Shared toy data are simulated/desensitized
- [x] Each module can be run from its own folder
- [x] Root README explains data privacy clearly
- [x] Figures and results are generated from toy inputs
- [x] MIT license added
- [x] Git status checked before publishing

Final recommended manual check before pushing: run `Rscript scripts/run_all_demos.R` in the intended public R environment after installing all required CRAN/Bioconductor packages.
