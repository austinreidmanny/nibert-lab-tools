#!/bin/bash

function usage() {
    #==== FUNCTION ====================================================================================================#
    #        NAME: usage
    # DESCRIPTION: setup a usage statement that will inform the user how to correctly invoke the
    #              program
    #==================================================================================================================#

    echo -e "\nERROR: Missing input files. Make sure to provide path to taxonomy and counts files \n\n" \
            "Usage: $0 -t ./sample.taxonomy.txt -c sample.counts.txt \n\n" \
            "Optional parameters: \n" \
                  "-o (directory for saving the files at the end; [default=current folder]) \n" \
                  "-s (sample name for the resulting merged table; [default=name according to taxonomy filename] \n\n" \
            "Example of a complex run: \n" \
            "$0 -t ./sample.taxonomy.txt -c sample.counts.txt -o analysis/merged_table/files -s sample001 \n\n" \
            "Exiting program. Please retry with correct parameters... \n" >&2; exit 1;
}

#======================================================================================================================#
# Process all user parameters and setup the computational environment
#======================================================================================================================#

    #==================================================================================================================#
    # Make sure the tool is invoked correctly, with taxonomy and counts tables
    #==================================================================================================================#
    while getopts "t:c:o:s:" arg;
        do
        	case ${arg} in

                t ) # Path to taxonomy table
                    taxonomy_table=${OPTARG}
                        ;;

                c ) # Path to table containing counts of mapped-reads
                    counts_table=${OPTARG}
                        ;;

                o ) # Path to output directory for saving final files
                    output_directory=${OPTARG}
                        ;;

                s ) # Optional sample name, for naming the resulting table
                    sample=${OPTARG}
                        ;;

                * ) # Display help
        		    usage
        		     	;;
        	esac
        done; shift $(( OPTIND-1 ))

    #==================================================================================================================#
    # Make sure required input files exist
    #==================================================================================================================#
    if [[ ! -f ${taxonomy_table} ]] || [[ ! -f ${counts_table} ]]; then
        usage
    fi

    #==================================================================================================================#
    # If user provides output directory, ensure it is valid, then create it
    #==================================================================================================================#
    # If output directory is provided, make sure it exists; if not provided, just use current dir
    if [[ ! -z "${output_directory}" ]]; then
    	mkdir -p ${output_directory}
    fi

    if [[ -z "${output_directory}" ]]; then
        output_directory="./"
    fi

    #==================================================================================================================#
    # If user did not provide sample name, make one based off the taxonomy file
    #==================================================================================================================#
    if [[ ! -z "${sample}" ]]; then
        sample=$(basename ${taxonomy_table} | cut -d "." -f 1)
    fi

    #==================================================================================================================#
    # Print all of that to the user
    #==================================================================================================================#
    echo -e "This tool will append a per-sequence mapped-reads count as a new column onto a taxonomy table. \n\n" \
            "Taxonomy table: ${taxonomy_table} \n" \
            "Mapped-reads counts table: ${counts_table} \n" \
            "Sample name: ${sample}"

#======================================================================================================================#
# Create the merged taxonomy+counts table
#======================================================================================================================#
# Make a header on the output file
echo -e "sequence_name	taxon_id	e_value	superkingdom	kingdom	phylum	class	order	family	genus	species	sequence_length	number_mapped_reads" > \
        ${sample}.taxonomy-and-counts.txt

# Keep all of the taxonomy table; for the counts table just retain contig name, length, & mapped-reads counts
join \
    -t $'\t' \
    ${taxonomy_table} \
    <(cat ${counts_table} | cut -f1,2,3) >> \
    ${sample}.taxonomy-and-counts.txt

# Tell user it's complete
echo -e "Merge complete \n" \
        "Final table at ${output_directory}/${sample}.taxonomy-and-counts.txt"
#=======================================================================================================================#

