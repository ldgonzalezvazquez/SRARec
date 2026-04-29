Organized pipeline structure

inputs/: source files
scripts/: scripts numbered by execution order
outputs/: generated results

Recommended order:
  python3 scripts/0_build_reconstructed_gamete_aas.py ...
  python3 scripts/1_export_phase_table_tsv.py ...
  python3 scripts/2_plot_s_recomb_rate_with_regions.py ...
  python3 scripts/3_plot_global_circos_annotated_zoom.py ...
  python3 scripts/4_plot_phase_maps_from_table.py ...

Methodological note:
- A small mismatch may exist between Circos arcs and the recombination-rate line.
- The recombination-rate line comes from a genome-wide global run.
- Circos arcs come from an S-gene cut and do not account for polymorphisms outside Spike.
