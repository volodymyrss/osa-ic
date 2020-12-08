# osa-ic

example usage:

```bash
osa-ic \
  -v dev201201 \
  -b /isdc/arc/rev_3/ \
  -f ic_list_combined.txt \
  /isdc/arc/rev_3/ic/ibis/rsp/isgr_ebds_mod_0001.fits
```

`-v dev201201` - produce version tag

`-b ...` - bare IC, used as basis to build the new version

`-f ...` - file with list of IC-ingestable files

`...` any other arguments are additional IC-ingestable files
