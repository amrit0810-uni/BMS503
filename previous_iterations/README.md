# Previous Iterations

Zip archives of earlier pipeline versions, kept here for reference and rollback.

## How to add a new archive

Before making major changes, export the current state:

```bash
cd "/home/vdiuser/Downloads/BMS503 (120526-one env)"
git archive --format=zip --prefix="BMS503-pipeline/" HEAD \
    -o previous_iterations/BMS503-pipeline-vX.X.zip
git add previous_iterations/
git commit -m "Archive vX.X before <description of upcoming changes>"
git push
```

## Contents

| File | Version | Description |
|------|---------|-------------|
| _(add entries here as you upload archives)_ | | |
