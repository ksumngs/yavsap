const express = require('express');
const app = express();
const path = require('path');
const fs = require('fs');

app.get('/', function (req, res) {
    res.sendFile(path.join(__dirname+'/index.html'));
})

app.get('/favicon.ico', function(req, res) {
    res.sendFile(path.join(__dirname+'/favicon.ico'));
})

app.use('/data', express.static(__dirname + '/data'));

app.use('/js/igv', express.static(__dirname + '/node_modules/igv/dist'));
app.use('/js/jquery', express.static(__dirname + '/node_modules/jquery/dist'));

app.use('/multiqc_data', express.static(__dirname + '/multiqc_data'));

app.get('/samples', function(req, res) {
    files = fs.readdirSync(path.join(__dirname + '/data'));
    contig_files = files.filter(file => file.endsWith('.bam'));
    res.setHeader('Content-Type', 'application/json');
    res.end(JSON.stringify(contig_files));
})

app.get('/sample', function(req, res) {
    res.sendFile(path.join(__dirname+'/sample.html'));
})

app.get('/reference', function(req, res) {
    files = fs.readdirSync(path.join(__dirname + '/data'));
    fasta_files = files.filter(file => file.endsWith('.fasta'));
    res.setHeader('Content-Type', 'application/json');
    res.end(JSON.stringify(fasta_files))
})

app.get('/multiqc', function(req, res) {
    res.sendFile(path.join(__dirname+'/multiqc_report.html'));
})

const port = process.env.PORT || 3000;
app.listen(port, function() {
    console.log('jev-analysis-pipeline results visualizer running, available at http://localhost:' + port)
})
