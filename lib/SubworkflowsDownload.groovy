//
// This file holds several functions that help reformat the downloaded reference
// genomes
//

class SubworkflowsDownload {

    /**
     * Gets the index that a fasta record occupies within a list of identifiers
     *
     * @param fastaText A string of a fasta record that the identifier will be parsed from
     * @param accessionNumbers The list of identifiers that `fastaText` will be indexed against
     * @return An integer index of where `fastaText` is found within `accessionNumbers`. If no match is found, returns -1.
     */
    public static Integer sortSequence(String fastaText, List accessionNumbers) {
        String id = fastaText.split(' ')[0].replace('>', '')
        Integer i = 0
        for (accessionNumber in accessionNumbers) {
            if (accessionNumber == id) {
                return i
            }
            i += 1
        }
        return -1
    }

    /**
     * Transforms a list by concatenating them with pipes
     *
     * @param accessionNumbers A list of accession numbers to concatenate
     * @return The joined accession numbers
     */
    public static String accessionTransform(String[] accessionNumbers) {
        return accessionNumbers.join('|')
    }
}
