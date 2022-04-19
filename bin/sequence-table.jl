#!/usr/bin/env julia
using ArgParse
using BioAlignments
using BioSequences
using CSV
using DataFrames
using EzXML
using FASTX
using JSON3
using Kelpie
using XAM

include("jl-lib/structs.jl")
include("jl-lib/functions.jl")
include("jl-lib/arguments.jl")
include("jl-lib/static-divs.jl")

const BIN_DIR = dirname(abspath(PROGRAM_FILE))
const LIB_DIR = joinpath(BIN_DIR, "jl-lib")
const STYLESHEET = joinpath(LIB_DIR, "yavsap.css")
const PHYLOTREE_UTILS = joinpath(LIB_DIR, "phylotree-utils.js")
const PHYLOTREE_COLOR = joinpath(LIB_DIR, "phylotree-colorizer.js")

args = parse_arguments()

strains_file = args["strainstable"]
aligned_file = args["alignments"]
reference_file = args["reference"]
include_multiqc = !args["no-multiqc"]
include_krona = !args["no-krona"]
include_tree = !isnothing(args["newick"])
tree_contents = nothing

reference_reader = open(FASTA.Reader, reference_file)
reference_records = collect(reference_reader)
close(reference_reader)

if length(reference_records) > 1
    @warn "More than one genome (chromosome) found in $reference_file. Only the first will be used"
end #if

reference_name = FASTA.identifier(first(reference_records))
reference_sequence = FASTA.sequence(first(reference_records))

strains_table = CSV.read(
    strains_file,
    DataFrame;
    header=["sample", "accession", "strain", "haplotype", "frequency"],
)

alignment_reader = open(SAM.Reader, aligned_file)
alignment_records = collect(alignment_reader)
close(alignment_reader)

if include_tree
    tree_contents = first(readlines(args["newick"]))
end #if

report_sections = cat(
    [
        ReportSection("summary", "Summary", "fas fa-list-ul"),
        if include_multiqc
            ReportSection(
                "multiqc", "Read Quality and Trimming", "fas fa-magnifying-glass"
            )
        else
            []
        end,
        if include_krona
            ReportSection("krona", "Read Classification (Kraken)", "fas fa-chart-pie")
        else
            []
        end,
        ReportSection("igv", "Alignments", "fas fa-bars-staggered"),
        if include_tree
            ReportSection("phylogenetics", "Phylogenetics", "fas fa-code-branch")
        else
            []
        end,
        ReportSection("nextflow", "Nextflow Report", "fas fa-shuffle"),
    ]...;
    dims=1,
)

custom_style_sheet = style(read(STYLESHEET, String))
phylotree_utils_script = read(PHYLOTREE_UTILS, String)
phylotree_colorizer_script = read(PHYLOTREE_COLOR, String)

navbar = nav(
    a(
        i(""; class="fab fa-bootstrap");
        href="#top",
        class="d-block p-3 link-dark text-decoration-none text-white text-center border-bottom",
        title="YAVSAP Report",
        data_bs_toggle="tooltip",
        data_bs_placement="right",
    ),
    ul(
        [
            li(
                a(
                    i(""; class=report_section.fontawesome);
                    href="#$(report_section.href)",
                    class="nav-link py-3 border-bottom text-light",
                    title=report_section.title,
                    data_bs_toggle="tooltip",
                    data_bs_placement="right",
                );
                class="nav-item",
            ) for report_section in report_sections
        ]...;
        class="nav nav-pills nav-flush flex-column mb-auto text-center",
    ),
    html_div(
        a(
            span(
                i(""; class="fas fa-circle-question");
                data_bs_toggle="modal",
                data_bs_target="#help-dialog",
            );
            href="#",
            class="d-block p-3 link-light text-decoration-none text-light text-center",
            title="Help",
            data_bs_toggle="tooltip",
            data_bs_placement="right",
        );
        class="dropdown border-top",
    );
    class="d-flex flex-column flex-shrink-0 bg-dark",
    style="width: 4.5rem",
)

genome_table = html_div(
    table(
        thead(
            tr(
                th("Sample"; colspan=3),
                th("Genotype"; colspan=2),
                th("Sequence"; colspan=4),
                [td(string(i); colspan=5) for i in 5:5:length(reference_sequence)]...,
            ),
        ),
        tbody(
            tr(
                th("Reference"; colspan=4),
                td(
                    a(
                        reference_name;
                        href="https://ncbi.nlm.nih.gov/nuccore/$reference_name",
                    ),
                ),
                [td(base; class=base) for base in reference_sequence]...,
            ),
            cat(
                [
                    sample_rows(s, reference_sequence, alignment_records, strains_table) for
                    s in unique(strains_table.sample)
                ]...;
                dims=1,
            )...,
        );
        id="genome-view",
        class="table",
    );
    id="genome-wrapper",
)

igv_options = Dict(
    "reference" => Dict(
        "id" => "reference",
        "fastaURL" => "reference.fasta",
        "indexURL" => "reference.fasta.fai",
    ),
    "tracks" => [
        Dict(
            "type" => "alignment",
            "format" => "bam",
            "url" => "$samplename.bam",
            "indexURL" => "$samplename.bam.bai",
            "name" => samplename,
        ) for samplename in unique(strains_table.sample)
    ],
)

phylogenetic_scripts = []
if include_tree
    phylogenetic_scripts = [
        script(
            "";
            src="https://cdn.jsdelivr.net/npm/d3@5.16.0/dist/d3.min.js",
            integrity="sha256-Xb6SSzhH3wEPC4Vy3W70Lqh9Y3Du/3KxPqI2JHQSpTw=",
            crossorigin="anonymous",
        ),
        script(
            "";
            src="https://cdn.jsdelivr.net/npm/underscore@1.13.2/underscore-umd.min.js",
            integrity="sha384-6URC9+r9R/tql/uTNEHRNRXFG53gIbvGSowjGwSBcHeQJfuL3QdHF+NsSgWlzqsr",
            crossorigin="anonymous",
        ),
        script(
            "";
            src="https://cdn.jsdelivr.net/npm/phylotree@1.0.13/dist/phylotree.min.js",
            integrity="sha256-y6vYUVnYZ4w+6dZofvpSL2ItmXVGRiN56p5bovmu0Bw=",
            crossorigin="anonymous",
        ),
        script("""
        $phylotree_utils_script
        newick = "$tree_contents";
        $phylotree_colorizer_script
        """),
    ]
end #if

middle_sections = []
if include_multiqc
    push!(
        middle_sections,
        section(
            iframe(""; src="multiqc_report.html", class="iframe-content");
            id="multiqc",
            class="min-vh-100",
        ),
    )
end
if include_krona
    push!(
        middle_sections,
        section(
            iframe(""; src="krona.html", class="iframe-content");
            id="krona",
            class="min-vh-100",
        ),
    )
end

EzXML.prettyprint(
    html(
        head(
            meta(; charset="utf8"),
            meta(; name="viewport", content="width=device-width, initial-scale=1"),
            link(;
                href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css",
                rel="stylesheet",
                integrity="sha384-1BmE4kWBq78iYhFldvKuhfTAU6auU8tT94WrHftjDbrCEXSU1oBoqyl2QvZ6jIW3",
                crossorigin="anonymous",
            ),
            link(;
                href="https://cdn.jsdelivr.net/npm/phylotree@1.0.13/dist/phylotree.css",
                rel="stylesheet",
                integrity="sha256-OMGYTWHSt3pK8AhnFhc18bkINhqrWz0+srfAkIi1Jdg=",
                crossorigin="anonymous",
            ),
            custom_style_sheet,
            title("YAVSAP Results"),
        ),
        body(
            help_modal,
            html_div(
                navbar,
                main(
                    header(h1("YAVSAP Results"); id="top", class="container py-3"),
                    section(
                        h2("Summary"),
                        genome_table;
                        id="summary",
                        class="container min-vh-100",
                    ),
                    middle_sections...,
                    section(""; id="igv", class="min-vh-100"),
                    include_tree ? phylogenetic_section : nothing,
                    section(
                        iframe(""; src="nextflow_report.html", class="iframe-content");
                        id="nextflow",
                        class="min-vh-100",
                    );
                    style="overflow-y: scroll",
                );
                class="d-flex vh-100",
            ),
            script(""; src="https://cdn.jsdelivr.net/npm/jquery@3.2.1/dist/jquery.min.js"),
            script(
                "";
                src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js",
            ),
            script(
                "";
                src="https://cdn.jsdelivr.net/npm/jquery-freeze-table@1.3.0/dist/js/freeze-table.min.js",
            ),
            script(""; src="https://kit.fontawesome.com/8f147eccd6.js"),
            script(
                "";
                src="https://cdn.jsdelivr.net/npm/igv@2.11.0/dist/igv.min.js",
                integrity="sha256-sr6GZtbybttUnYJHVKTjxA/aj9zru7lgZnRUOV3o7Gc=",
                crossorigin="anonymous",
            ),
            script("""
            var tooltipTriggerList = [].slice.call(
                document.querySelectorAll('[data-bs-toggle="tooltip"]')
            );
            var tooltipList = tooltipTriggerList.map(function (tooltipTriggerEl) {
                return new bootstrap.Tooltip(tooltipTriggerEl);
            });
            """),
            script("""
            \$("#genome-wrapper").freezeTable({
                columnNum: 2,
                scrollable: true,
                shadow: true,
            });
            """),
            script("""
            options = $(JSON3.write(igv_options));
            igv
                .createBrowser(document.getElementById("igv"), options)
                .then(function (browser) {
                    igv.browser = browser;
                });
            """),
            phylogenetic_scripts...,
        ),
    ),
)
