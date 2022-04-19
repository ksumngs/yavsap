
function parse_arguments()
    s = ArgParseSettings()

    # Disable Julia formatter as it doesn't understand the nested table syntax of ArgParse
    #! format: off

    @add_arg_table! s begin
        "strainstable"
            help = "TSV of strains"
            required = true
            arg_type = String
            range_tester = x -> isfile(x)
        "alignments"
            help = "SAM of haplotypes"
            required = true
            arg_type = String
            range_tester = x -> isfile(x)
        "reference"
            help = "FASTA of reference genome"
            required = true
            arg_type = String
            range_tester = x -> isfile(x)
        "--multiqc"
            help = "Do include a MultiQC report section"
            action = :store_false
        "--no-multiqc"
            help = "Do not include a MultiQC report section"
            action = :store_true
        "--krona"
            help = "Do include a Krona chart"
            action = :store_false
        "--no-krona"
            help = "Do not include a Krona chart"
            action = :store_true
        "--newick"
            help = "A file that containing the phylogenetic tree"
            arg_type = String
            range_tester = x -> isfile(x)
    end #add_arg_table

    return parse_args(s)
end #function
