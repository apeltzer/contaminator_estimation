# Images

The metro map is generated with [nf-metro](https://github.com/seqeralabs/nf-metro) from `contaminationqc_metro_map.mmd`. After editing the `.mmd`, regenerate all renders from this directory:

```bash
uvx nf-metro@2.1.0 render contaminationqc_metro_map.mmd -o contaminationqc_metro_map.svg
uvx nf-metro@2.1.0 render contaminationqc_metro_map.mmd --mode light -o contaminationqc_metro_map_light.png
uvx nf-metro@2.1.0 render contaminationqc_metro_map.mmd --mode dark -o contaminationqc_metro_map_dark.png
```

The SVG (`manifest.diagram`) switches between light and dark palettes with the viewer's colour scheme; the README uses the PNGs in a `<picture>` element.
