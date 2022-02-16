#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/// summary: |
///   Takes a list of reads files, skips every `skip` entry in the list, and returns
///   a new list of file objects
def skipping_read(List files, Integer skip) {
    readFiles = [];
    for (int i = 0; i <= files.size(); i = i+skip) {
        if (files[i]?.trim()) {
            filepath = file(files[i])
            if ( filepath instanceof List ) {
                filepath.each { readFiles.add(it) }
            }
            else {
                readFiles.add(filepath)
            }
        }
    }
    return readFiles
}
