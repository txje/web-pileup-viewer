// a not-so-elegant but critical piece of scoping and closure
function attacher() {
    var o = arguments[0];
    var f = arguments[1];
    var params = [];
    for(var i = 2; i < arguments.length; i++)
        params.push(arguments[i]);
    return function() {
        var newparams = [];
        for(var i in arguments)
            newparams.push(arguments[i]);
        return f.apply(o, params.concat(newparams));
    }
}

function Socket(url, callback) {
    var sock = new WebSocket(url);
    
    sock.onopen = function(evt) {
        callback('open', '');
    };
    sock.onclose = function(evt) {
        callback('close', '');
    };
    sock.onmessage = function(evt) {
        callback('msg', evt.data);
    };
    sock.onerror = function(evt) {
        callback('err', evt.data);
    };
    
    this.send = function(msg) {
        sock.send(msg);
    }
    
    this.close = function() {
        sock.close();
    }
}

function Console(eid) {
    var out = $('#'+eid);
    this.write = function(txt, color) {
        out.append($('<br/>'));
        var div = $('<div>');
        div.css('display', 'inline');
        div.append('> ' + txt);
        if(color) {
            div.css('color', color);
        }
        out.append(div);
    }
}

function PileupDisplay(eid, app) {
    var $canvas = $('#'+eid);
    var ctx = $canvas.get(0).getContext('2d');
    
    this.drawUPReads = function(start, end, reads) {
        var w = $canvas.width();
        var h = $canvas.height();
        $canvas.attr('width', w);
        $canvas.attr('height', h);
        var bpp = (end - start) / w;
        var read_height = Math.max(Math.min(10, parseInt(1/bpp * 10)), 1); // minimum height of 1, maximum height of 10
        // assuming reads are in genomic order
        var read_rows = []; // read (specified by index) taking up each row
        ctx.save();
        ctx.strokeStyle = 'rgba(255,255,255,0.5)';
        ctx.lineWidth = 2;
        ctx.lineJoin = 'round';
        for(var i = 0; i < reads.length; i++) {
            if(reads[i][0] > end)
                break;
            if(reads[i][0]+reads[i][1] < start)
                continue;
            
            // figure out which row to place it
            var row = null;
            for(var r = 0; r < read_rows.length; r++) {
                var oar = reads[read_rows[r]]; // other aligned read
                if(reads[i][0] > oar[0] + oar[1]) {
                    row = r+1;
                    break;
                }
            }
            if(row == null)
                row = r+1; // rows start at 1, not 0
            
            read_rows[row-1] = i;
            
            // draw it
            if(read_height == 1) {
                ctx.beginPath();
                ctx.moveTo((reads[i][0]-start)/bpp, h-row*5);
                ctx.lineTo((reads[i][0]+reads[i][1]-start)/bpp, h-row*5);
                ctx.stroke();
            }
            else {
                ctx.strokeRect((reads[i][0]-start)/bpp, h-row*(read_height+4), reads[i][1]/bpp, read_height);
            }
        }
        ctx.restore();
    }
    
    this.drawPaired = function(start, end, pairs) {
        var w = $canvas.width();
        var h = $canvas.height();
        $canvas.attr('width', w); // clears the canvas
        $canvas.attr('height', h);
        var bpp = (end - start) / w;
        var pair_height = Math.max(Math.min(10, parseInt(1/bpp * 10)), 1); // minimum height of 1, maximum height of 10
        // assuming reads are in genomic order
        var pair_rows = []; // read (specified by index) taking up each row
        ctx.save();
        ctx.strokeStyle = 'rgba(255,255,255,0.5)';
        ctx.lineWidth = 2;
        ctx.lineJoin = 'round';
        for(var i = 0; i < pairs.length; i++) {
            var read0 = pairs[i][0];
            var read1 = pairs[i][1];
            if(read0[0] > end)
                break;
            if(read1[0]+read1[1] < start)
                continue;
            
            // figure out which row to place it
            var row = null;
            for(var r = 0; r < pair_rows.length; r++) {
                // second read of other aligned pair (or the only read if unpaired)
                var oar = pairs[pair_rows[r]][1];
                if(!oar)
                    oar = pairs[pair_rows[r]][0];
                if(read0[0] > oar[0] + oar[1]) {
                    row = r+1;
                    break;
                }
            }
            if(row == null)
                row = r+1; // rows start at 1, not 0
            
            pair_rows[row-1] = i;
            
            if(read1) { // paired
                var darkcolor = 'rgba(50,200,50,0.5)';
                var lightcolor = 'rgba(50,200,50,0.25)';
            }
            else { // unpaired
                var darkcolor = 'rgba(200,50,50,0.5)';
                var lightcolor = 'rgba(200,50,50,0.25)';
            }
            
            // draw it
            if(pair_height == 1) {
                ctx.beginPath();
                ctx.strokeStyle = darkcolor;
                ctx.moveTo((read0[0]-start)/bpp, h-row*5);
                ctx.lineTo((read0[0]+read0[1]-start)/bpp, h-row*5);
                ctx.stroke();
                
                ctx.strokeStyle = lightcolor; // different color for the insert
                ctx.beginPath();
                ctx.moveTo((read0[0]+read0[1]-start)/bpp, h-row*5);
                ctx.lineTo((read1[0]-start)/bpp, h-row*5);
                ctx.stroke();
                
                ctx.strokeStyle = darkcolor;
                ctx.beginPath();
                ctx.moveTo((read1[0]-start)/bpp, h-row*5);
                ctx.lineTo((read1[0]+read1[1]-start)/bpp, h-row*5);
                ctx.stroke();
            }
            else {
                ctx.strokeStyle = darkcolor;
                ctx.strokeRect((read0[0]-start)/bpp, h-row*(pair_height+4), read0[1]/bpp, pair_height);
                ctx.beginPath();
                ctx.strokeStyle = lightcolor;
                ctx.moveTo((read0[0]+read0[1]-start)/bpp, h-row*(pair_height+4)+pair_height/2);
                ctx.lineTo((read1[0]-start)/bpp, h-row*(pair_height+4)+pair_height/2);
                ctx.stroke();
                ctx.strokeStyle = darkcolor;
                ctx.strokeRect((read1[0]-start)/bpp, h-row*(pair_height+4), read1[1]/bpp, pair_height);
            }
        }
        ctx.restore();
    }
    
    var dragging = false;
    var drag_start;
    var imageData;
    
    $canvas.mousedown(attacher(this, function(event) {
        dragging = true;
        drag_start = event.pageX;
        // cache current image
        imageData = ctx.getImageData(0, 0, $canvas.width(), $canvas.height());
    }));
    
    $canvas.mouseup(attacher(this, function(event) {
        dragging = false;
        var dist = event.pageX - drag_start;
        var percent = -1 * dist / $canvas.width();
        pan(percent);
    }));
    
    $canvas.mousemove(attacher(this, function(event) {
        if(dragging) {
            var dist = event.pageX - drag_start;
            ctx.clearRect(0, 0, $canvas.width(), $canvas.height());
            ctx.putImageData(imageData, dist, 0);
        }
    }));
}

function App() {
    var c = new Console('console');
    var seqs = [];
    var disp = new PileupDisplay('graphics', this);

    this.connected = false;
    
    var start = null;
    var end = null;

    var callback = function(evt, msg) {
        // reset the timeout
        clearTimeout(timer);
        timer = setTimeout(disconnect, TIMEOUT);
        
        if(evt == 'open') {
            c.write('connected.');
            this.connected = true;
            $('#conn').text(1);
            s.send('refs');
        }
        else if(evt == 'msg') {
            console.log(msg);
            
            var fields = msg.split('||');
            
            // get reference (chromosome) listing
            if(fields[0] == 'refs') {
                eval('seqs = ' + fields[1]);
                for(var i = seqs.length-1; i >= 0; i--) {
                    // remove NTs
                    if(seqs[i][0].substring(0, 2) == 'NT')
                        seqs.splice(i, 1);
                }
                // sort numerical refs by number
                seqs.sort(function(a,b){
                    var ia = parseInt(a[0]);
                    var ib = parseInt(b[0]);
                    if(isNaN(ia)) {
                        if(isNaN(ib)) {
                            return -1;
                        }
                        return 1;
                    }
                    else if (isNaN(ib))
                        return -1;
                    else
                        return ia-ib;
                });
                c.write('reference sequences (excluding NT): ' + seqs.length);
                
                // add these to the refs select box
                var $sel = $('#refs');
                for(var i = 0; i < seqs.length; i++) {
                    var $op = $('<option>');
                    $op.val(seqs[i][0]);
                    var len = seqs[i][1];
                    if(len > 1000000)
                        len = (len/1000000).toFixed(2) + ' Mb';
                    else if(len > 1000)
                        len = (len/1000).toFixed(2) + ' Kb';
                    else
                        len = len + ' b';
                    $op.append(seqs[i][0] + ' ('+len+')');
                    $sel.append($op);
                }
                
                update(); // show something initially
            }
            
            else if(fields[0] == 'reads') {
                eval('var reads = ' + fields[1]);
                c.write(reads.length + ' reads fetched.');
                disp.drawPaired(start, end, reads);
            }
            
            else if(fields[0] == 'conn') {
                $('#conn').text(fields[1]);
            }
            
            else if(fields[0] == 'err') {
                c.write(fields[1], '#990000');
            }
            
            else {
                c.write('unknown response (disconnecting): ' + msg);
                s.close();
            }
            
        }
        else if(evt == 'err') {
            c.write('error (disconnecting): ' + msg);
        }
        else if(evt == 'close') {
            c.write('disconnected.');
            $('#conn').text('[disconnected]');
            this.connected = false;
            // stop the timeout counter
            clearTimeout(timer);
        }
    }
    
    this.close = function() {
        s.close();
    }

    ws_server = 'localhost:8080'

    var s = new Socket("ws://" + ws_server + "/", attacher(this, callback));
    c.write('connecting...');
    
    this.get = function(chrom, spos, epos) {
        if(!this.connected)
            return false;
        start = spos;
        end = epos;
        c.write('fetching reads in chrom (' + chrom + ') ' + start + ' - ' + end);
        s.send('reads ' + chrom + ' ' + start + ' ' + end);
        return true;
    }
    
    this.clog = function(msg) {
        c.write(msg);
    }
}

var app;
var timer;
var TIMEOUT = 600000; // 10 min

function init() {
    app = new App();
    timer = setTimeout(disconnect, TIMEOUT); // timeout and automatically disconnect after 10 minutes of inactivity
}

function update() {
    var chrom = $('#refs').val();
    try {
        var start = $('#startpos').val();
        var end = $('#endpos').val();
    } catch(e) {
        app.clog('Start and end positions must contain only numbers (ex. 3000000, 150001000)');
        return;
    }
    app.get(chrom, start, end);
}

function pan(percent) {
    try {
        var start = parseInt($('#startpos').val());
        var end = parseInt($('#endpos').val());
    } catch(e) {
        app.clog('Start and end positions must contain only numbers (ex. 3000000, 150001000)');
        return;
    }
    var span = end - start;
    var topan = parseInt(span * percent);
    $('#startpos').val(start + topan);
    $('#endpos').val(end + topan);
    update();
}

function zoom(percent) {
    try {
        var start = parseInt($('#startpos').val());
        var end = parseInt($('#endpos').val());
    } catch(e) {
        app.clog('Start and end positions must contain only numbers (ex. 3000000, 150001000)');
        return;
    }
    var span = end - start;
    var middle = start + span/2;
    var newspan = span * percent;
    $('#startpos').val(parseInt(middle - newspan/2));
    $('#endpos').val(parseInt(middle + newspan/2));
    update();
}

function disconnect() {
    app.close();
}

window.addEventListener("load", init, false);
