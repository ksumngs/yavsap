const express = require('express');
const pug = require('pug');
const app = express();
const path = require('path');
const fs = require('fs');
const { randomWeibull } = require('d3');

app.set('views', path.join(__dirname, '/views'));
app.set('view engine', 'pug');

function getSampleList() {
    files = fs.readdirSync(path.join(__dirname + '/data'));
    bam_files = files.filter(file => file.endsWith('.bam'));
    sample_files = [];
    for (var i = 0; i < bam_files.length; i++) {
        b = bam_files[i];
        if (!b.includes('contigs')) {
            sample_name = b.replace('.bam', '');
            sample_files.push(sample_name);
        }
    }
    return sample_files;
}

function getReferenceGenomeName() {
    files = fs.readdirSync(path.join(__dirname + '/data'));
    fasta_files = files.filter(file => file.endsWith('.fasta'));
    return fasta_files[0].replace('.fasta', '');
}

app.get('/', function (req, res) {
    res.render('index', {refname: getReferenceGenomeName(), samples: getSampleList()});
})

app.get('/alignments/:sample', function(req, res) {
    sampleName = req.params.sample;
    if (!getSampleList().includes(sampleName)) {
        res.send(404);
    }
    igvOptions = {
        reference: {
            id: getReferenceGenomeName(),
            fastaURL: '/data/' + getReferenceGenomeName() + '.fasta'
        },
        tracks: [
            {
                type: 'alignment',
                format: 'bam',
                url: '/data/'+ sampleName + '.bam',
                indexURL: '/data/' + sampleName + '.bam.bai',
                name: sampleName
            },
            {
                type: 'alignment',
                format: 'bam',
                url: '/data/' + sampleName + '.contigs.bam',
                indexURL: '/data/' + sampleName + '.contigs.bam.bai',
                name: sampleName + '.contigs'
            }
        ]
    };
    res.render('alignment', { samplename: req.params.sample, options: igvOptions });
})

app.get('/favicon.ico', function(req, res) {
    res.sendFile(path.join(__dirname+'/favicon.ico'));
})

app.use('/data', express.static(__dirname + '/data'));

app.use('/js/igv', express.static(__dirname + '/node_modules/igv/dist'));
app.use('/js/jquery', express.static(__dirname + '/node_modules/jquery/dist'));
app.use('/js/big-integer', express.static(__dirname + '/node_modules/big-integer'));
app.use('/js/underscore', express.static(__dirname + '/node_modules/underscore'));
app.use('/js/spin', express.static(__dirname + '/node_modules/spin.js'));
app.use('/js/d3', express.static(__dirname + '/node_modules/d3'));
app.use('/js/cjson', express.static(__dirname + '/node_modules/circular-json/build'));
app.use('/js/cblob', express.static(__dirname + '/node_modules/canvas-toBlob'));
app.use('/js/filesaver', express.static(__dirname + '/node_modules/file-saver'));
app.use('/js/twbs', express.static(__dirname + '/node_modules/bootstrap/dist/js'));

// CSS Serving, from node_modules and locally
app.use('/css/twbs', express.static(__dirname + '/node_modules/bootstrap/dist/css'));
app.use('/css/fonts', express.static(__dirname + '/node_modules/bootstrap/dist/fonts'));
app.use('/css/local', express.static(__dirname + '/css'))

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
