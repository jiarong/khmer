import khmer, sys, threading, time
 
###
 
s = set()
def read_names(rparser, tnum):
    print 'started', tnum
 
    n = 0
    for n, read in enumerate(rparser):
        s.add(read.name)
 
        if n % 1000 == 0:
            print 'sleeping', tnum, n
            time.sleep(0.2)
 
    print 'done', tnum, 'got', n
 
###
 
filename = sys.argv[1]
n_threads = int(sys.argv[2])
 
config = khmer.get_config()
bufsz = config.get_reads_input_buffer_size()
config.set_reads_input_buffer_size(n_threads * 64 * 1024)
 
rparser = khmer.ReadParser(filename, n_threads)
 
print 'starting threads'
threads = []
for tnum in xrange(n_threads):
    print 'starting', tnum
    t = threading.Thread(target=read_names, args=(rparser, tnum ))
    threads.append(t)
    t.start()
 
for t in threads:
    t.join()
 
print 'done; loaded %s sequences' % len(s)
