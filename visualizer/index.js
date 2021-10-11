const express = require('express');
const pug = require('pug');
const app = express();
const path = require('path');
const fs = require('fs');
const { report } = require('process');

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
            sample_files.push({ samplename: sample_name, hastree: hasPhylogeneticTree(sample_name) });
        }
    }
    return sample_files;
}

function getSampleNames() {
    sample_names = [];
    for (var i = 0; i < getSampleList().length; i++) {
        sn = getSampleList()[i].samplename;
        sample_names.push(sn);
    }
    return sample_names;
}

function getReferenceGenomeName() {
    files = fs.readdirSync(path.join(__dirname + '/data'));
    fasta_files = files.filter(file => file.endsWith('.fasta'));
    return fasta_files[0].replace('.fasta', '');
}

function hasPhylogeneticTree(sample) {
    files = fs.readdirSync(path.join(__dirname + '/data'));
    tree_files = files.filter(file => file.endsWith('.tree'));
    return tree_files.includes(sample + '.tree');
}

function hasContigsAlignment(sample) {
    files = fs.readdirSync(path.join(__dirname + '/data'));
    contig_files = files.filter(file => file.endsWith('.contigs.bam'))
    return contig_files.includes(sample + '.contigs.bam');
}

function getTraceDocuments(prefix, suffix) {
    files = fs.readdirSync(path.join(__dirname + '/.trace'));
    trace_files = files.filter(file => file.startsWith(prefix) && file.endsWith(suffix));
    traces = [];
    for (var i = 0; i < trace_files.length; i++) {
        trace_file = trace_files[i];
        title = trace_file.replace(prefix, '').replace(suffix, '');
        traces.push( { file: trace_file, timestamp: title } );
    }
    return traces;
}

function getNextflowReports() {
    return getTraceDocuments('execution_report_', '.html');
}

function serveTraceDocument(prefix, suffix) {
    return function(req, res) {
        timestamp = req.params.timestamp;
        if (timestamp.toLowerCase() == 'latest') {
            trace_files = getTraceDocuments(prefix, suffix).map(doc => doc.file);
            trace_files.sort();
            latest_document = trace_files[trace_files.length - 1];
            res.sendFile(path.join(__dirname+'/.trace/'+latest_document));
        }
        else {
            if (! getTraceDocuments(prefix, suffix).map(doc => doc.timestamp).includes(timestamp)) {
                res.send(404);
            }
            else {
                res.sendFile(path.join(__dirname+'/.trace/'+prefix+timestamp+suffix));
            }
        }
    }
}

app.get('/', function (req, res) {
    res.render('index',
        {
            refname: getReferenceGenomeName(),
            samples: getSampleList(),
            nfreports: getNextflowReports()
        });
})

app.get('/alignments/:sample', function(req, res) {
    sampleName = req.params.sample;
    if (!getSampleNames().includes(sampleName)) {
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
            }
        ]
    };
    if (hasContigsAlignment(sampleName)) {
        igvOptions.tracks.push({
            type: 'alignment',
            format: 'bam',
            url: '/data/' + sampleName + '.contigs.bam',
            indexURL: '/data/' + sampleName + '.contigs.bam.bai',
            name: sampleName + '.contigs'
        })
    }
    res.render('alignment', { samplename: req.params.sample, options: igvOptions, hascontigs: hasContigsAlignment(sampleName) });
})

app.get('/phylogenetics/:sample', function(req, res) {
    sampleName = req.params.sample;
    if (!hasPhylogeneticTree(sampleName)) {
        res.send(404);
    }
    newickData = fs.readFileSync(path.join(__dirname+'/data/'+sampleName+'.tree'), {encoding: 'utf8', flag: 'r'}).trim();
    res.render('tree', {samplename: sampleName, sampletree: newickData});
})

app.get('/nf-report/:timestamp', serveTraceDocument('execution_report_', '.html'))

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

app.get('/multiqc', function(req, res) {
    res.sendFile(path.join(__dirname+'/multiqc_report.html'));
})

const port = process.env.PORT || 3000;
app.listen(port, function() {
    console.log('jev-analysis-pipeline results visualizer running, available at http://localhost:' + port)
})
