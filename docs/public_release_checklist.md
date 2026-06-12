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
- [x] R dependency installation helper added
- [x] GitHub Actions workflow added for full demo run
- [ ] Fresh-library reproducibility test completed

Final recommended manual check before pushing:

```bash
mkdir -p /tmp/sb_demo_r_lib
R_LIBS_USER=/tmp/sb_demo_r_lib Rscript scripts/install_r_dependencies.R
R_LIBS_USER=/tmp/sb_demo_r_lib Rscript scripts/create_shared_toy_data.R
R_LIBS_USER=/tmp/sb_demo_r_lib Rscript scripts/run_all_demos.R
```

This final check is intentionally left unchecked until it has been run in the intended public R environment.
