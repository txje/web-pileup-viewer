import pysam
import ws_server

MAX_REGION = 100000

class BamServ:
    def __init__(self, f):
        self.bam = pysam.Samfile(f, "rb")

    def serv(self, data):
        if data == 'refs':
            return 'refs||' + str([[self.bam.references[r], int(self.bam.lengths[r])] for r in xrange(self.bam.nreferences)])
        else:
            fields = data.split(' ')
            if fields[0] == 'reads':
                chrom = fields[1]
                start = int(fields[2])
                end = int(fields[3])
                if end < start:
                    return 'err||invalid coordinates'
                elif end - start > MAX_REGION:
                    return 'err||region too large'
                #try:
                reads = []
                for read in self.bam.fetch(chrom, start, end):
                    #reads.append([read.pos, read.rlen]) # simple unpaired reads
                    reads.append([read.pos, read.rlen, read.is_proper_pair, read.mpos])
                # pair up reads
                pairs = []
                while len(reads) > 0:
                    if not reads[0][2]: # not properly paired
                        pairs.append([reads.pop(0)[:2], 0])
                        continue
                    r = 1
                    paired = False
                    while r < len(reads):
                        if reads[0][0] == reads[r][3] and reads[0][3] == reads[r][0]: # found the matching pair
                            pairs.append([reads.pop(0)[:2], reads.pop(r-1)[:2]]) # the order here matters
                            paired = True
                            break
                        r += 1
                    if not paired: # didn't find matching pair, but we know it exists
                        pairs.append([reads[0][:2], [reads.pop(0)[3], 0]]) # make up this read with a known position and 0 length
                #return 'reads||'+str(reads)
                return 'reads||'+str(pairs)
                #except:
                #    return "unable to fetch reads"
            elif fields[0] == 'read':
                qname = fields[1]
                chrom = fields[2]
                start = fields[3]
                end = fields[4]
                for read in self.bam.fetch(chrom, start, end):
                    if read.qname == qname:
                        break
                else:
                    return "read not found"
                properties = ["aend", "alen", "bin", "cigar", "compare", "flag", # "fancy_str"
                    "is_duplicate", "is_paired", "is_proper_pair", "is_qcfail", "is_read1", "is_read2", "is_reverse", "is_secondary", "is_unmapped",
                    "isize", "mapq", "mate_is_reverse", "mate_is_unmapped", "mpos", "mrnm", "opt", "pos",
                    "qend", "qlen", "qname", "qqual", "qstart", "qual", "query", "rlen", "rname", "seq", "tags", "tid"]
                json = {}
                for p in properties:
                    json[p] = getattr(read, p)
                json = str(json).replace('True', 'true').replace('False', 'false').replace('None', 'null') # convert JSON to string and adjust from Python to Javascript booleans
                return json
                
            return 'echo(' + str(data) + ')'

def start(bf):

    bam = BamServ(bf)
    
    print 'server starting on port 8080.'

    server = ws_server.MyServer(8080, bam.serv)
    server.start()

if __name__ == "__main__":
    start('/playpen/C57BL.bam')
