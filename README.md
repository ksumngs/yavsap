# Japanese Encephalitis Virus Analysis Pipeline

A [Nextflow] pipeline for studying Japanese Encephalitis Virus (JEV)
populations within a single sample. :dna::computer::chart_with_upwards_trend:
Yeah, we're still looking for a better name. :shrug:

## Version History

### v0.1.0-alpha

The pipeline is still pretty rough at this point, but it gets the job done, and
the changes I want to make to it _shouldn't_ affect the public API in a breaking
way. At this point, the changes should be more modular, and it's time to switch
to a Git flow branching model.

## To-Do List

- [ ] Write the docs
- [ ] Cleanup the configuration (citing nf-core where needed)
- [ ] Make haplotyping actually work (maybe push it into its own repo?)
- [ ] Allow downloading of every output file from the visualizer
  - [ ] Will require outputing these files into more distinct output folders
- [ ] Add Filtlong to MultiQC (not this repo, I know)

[Nextflow]: https://nextflow.io
