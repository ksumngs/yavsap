const express = require('express');
const pug = require('pug');
const app = express();
const path = require('path');
const fs = require('fs');
const { report } = require('process');

app.set('views', path.join(__dirname, '/_views'));
app.set('view engine', 'pug');

// Useful functions for finding sample names and lists
function getSampleList() {
    files = fs.readdirSync(path.join(__dirname + '/alignment'));
    bam_files = files.filter(file => file.endsWith('.bam'));
    sample_files = [];
    for (var i = 0; i < bam_files.length; i++) {
        b = bam_files[i];
        sample_name = b.replace('.bam', '');
        sample_files.push({
            samplename: sample_name,
            hastree: hasPhylogeneticTree(sample_name),
            hascontigs: hasContigsAlignment(sample_name),
            hasvariants: hasVariantCalls(sample_name)
        });
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
    files = fs.readdirSync(path.join(__dirname + '/reference'));
    fasta_files = files.filter(file => file.endsWith('.fasta'));
    return fasta_files[0].replace('.fasta', '');
}

function hasOutputfile(sample, directory, suffix) {
    if (!fs.existsSync(__dirname + directory)) {
        return false;
    }
    files = fs.readdirSync(path.join(__dirname + directory));
    typed_files = files.filter(file => file.endsWith(suffix));
    return typed_files.includes(sample + suffix);
}

function hasPhylogeneticTree(sample) {
    return hasOutputfile(sample, '/phylogenetics', '.nwk');
}

function hasContigsAlignment(sample) {
    return hasOutputfile(sample, '/assembly/sequence', '.contigs.fasta');
}

function hasVariantCalls(sample) {
    return hasOutputfile(sample, '/variants', '.vcf');
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

function getNextflowTimelines() {
    return getTraceDocuments('execution_timeline_', '.html');
}

function getNextflowTraces() {
    return getTraceDocuments('execution_trace_', '.txt');
}

function getNextflowGraphs() {
    return getTraceDocuments('pipeline_dag_', '.svg');
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

function didFiltering() {
    return fs.existsSync(path.join(__dirname + '/classification'));
}

function isOntResults() {
    if (fs.existsSync(__dirname + '/haplotypes')) {
        files = fs.readdirSync(path.join(__dirname + '/haplotypes'));
        typed_files = files.filter(file => file.endsWith('.yaml'));
        yaml_extension = typed_files.length > 0;
    }
    else {
        yaml_extension = true;
    }
    return fs.existsSync(__dirname + '/variants') || yaml_extension
}

// Server page rendering
app.get('/', function (req, res) {
    res.render('index',
        {
            title: getReferenceGenomeName(),
            refname: getReferenceGenomeName(),
            platform: isOntResults() ? 'ont' : 'pe',
            samples: getSampleList(),
            nfreports: getNextflowReports(),
            nftimelines: getNextflowTimelines(),
            nftraces: getNextflowTraces(),
            nfdags: getNextflowGraphs(),
            haskraken: didFiltering()
        });
})

app.get('/alignments/:sample', function(req, res) {
    sampleName = req.params.sample;
    if (!getSampleNames().includes(sampleName)) {
        res.sendStatus(404);
    }
    igvOptions = {
        reference: {
            id: getReferenceGenomeName(),
            fastaURL: '/reference/' + getReferenceGenomeName() + '.fasta'
        },
        tracks: [
            {
                type: 'alignment',
                format: 'bam',
                url: '/alignment/'+ sampleName + '.bam',
                indexURL: '/alignment/' + sampleName + '.bam.bai',
                name: sampleName
            }
        ]
    };
    res.render('alignment',
        {
            title: req.params.sample + ' Alignments',
            samplename: req.params.sample,
            options: igvOptions,
            hascontigs: hasContigsAlignment(sampleName)
        });
})

app.get('/phylogenetics/:sample', function(req, res) {
    sampleName = req.params.sample;
    if (!hasPhylogeneticTree(sampleName)) {
        res.send(404);
    }
    newickData = fs.readFileSync(path.join(__dirname+'/phylogenetics/'+sampleName+'.nwk'), {encoding: 'utf8', flag: 'r'}).trim();
    res.render('tree',
        {
            title: sampleName + ' Phylogenetic Tree',
            samplename: sampleName,
            sampletree: newickData
        });
})

app.get('/nf-report/:timestamp', serveTraceDocument('execution_report_', '.html'))
app.get('/nf-timeline/:timestamp', serveTraceDocument('execution_timeline_', '.html'))
app.get('/nf-trace/:timestamp', serveTraceDocument('execution_trace_', '.txt'))
app.get('/nf-dag/:timestamp', serveTraceDocument('pipeline_dag_', '.svg'))
app.get('/multiqc', function(req, res) {
    res.sendFile(path.join(__dirname+'/multiqc_report.html'));
})

// Static file serving
app.get('/favicon.ico', function(req, res) {
    res.sendFile(path.join(__dirname+'/favicon.ico'));
})


// Static directory serving
app.use('/reference', express.static(__dirname + '/reference'));
app.use('/classification', express.static(__dirname + '/classification'));
app.use('/alignment', express.static(__dirname + '/alignment'));
app.use('/assembly', express.static(__dirname + '/alignment'));
app.use('/assembly/alignment', express.static(__dirname + '/assembly/alignment'));
app.use('/assembly/sequence', express.static(__dirname + '/assembly/sequence'));
app.use('/variants', express.static(__dirname + '/variants'));
app.use('/haplotypes', express.static(__dirname + '/haplotypes'));
app.use('/multi_alignment', express.static(__dirname + '/multi_alignment'));
app.use('/phylogenetics', express.static(__dirname + '/phylogenetics'));
app.use('/multiqc_data', express.static(__dirname + '/multiqc_data'));

// Javascript serving
app.use('/js/igv', express.static(__dirname + '/node_modules/igv/dist'));
app.use('/js/jquery', express.static(__dirname + '/node_modules/jquery/dist'));
app.use('/js/twbs', express.static(__dirname + '/node_modules/bootstrap/dist/js'));
app.use('/js/popper', express.static(__dirname + '/node_modules/@popperjs/core/dist/umd'));
app.use('/js/fa', express.static(__dirname + '/node_modules/@fortawesome/fontawesome-free/js'));
app.use('/js/d3', express.static(__dirname + '/node_modules/d3/dist'));
app.use('/js/underscore', express.static(__dirname + '/node_modules/underscore'));
app.use('/js/lodash', express.static(__dirname + '/node_modules/lodash'));
app.use('/js/phylotree', express.static(__dirname + '/node_modules/phylotree/dist'));

// CSS Serving, from node_modules and locally
app.use('/css/twbs', express.static(__dirname + '/node_modules/bootstrap/dist/css'));
app.use('/css/local', express.static(__dirname + '/_css'))
app.use('/css/fa', express.static(__dirname + '/node_modules/@fortawesome/fontawesome-free/css'));
app.use('/css/webfonts', express.static(__dirname + '/node_modules/@fortawesome/fontawesome-free/webfonts'));
app.use('/css/phylotree', express.static(__dirname + '/node_modules/phylotree/dist'));

// App startup
const port = process.env.PORT || 3000;
app.listen(port, function() {
    console.log('YAVSAP results visualizer running, available at http://localhost:' + port)
})
