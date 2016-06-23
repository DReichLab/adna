CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function
PROG=adna-trim adna-ldup

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

adna-trim:kthread.o trim.o
		$(CC) $(CFLAGS) $^ -o $@ -lz -lm -lpthread

adna-ldup:ldup.o bgzf.o hts.o sam.o
		$(CC) $(CFLAGS) $^ -o $@ -lz -lm -lpthread

bgzf.o:bgzf.c bgzf.h khash.h
		$(CC) -c $(CFLAGS) $(DFLAGS) -DBGZF_MT $(INCLUDES) bgzf.c -o $@

ldup.o:ldup.c sam.h bgzf.h hts.h kdq.h khash.h
		$(CC) -c $(CFLAGS) $(DFLAGS) -DBGZF_MT $(INCLUDES) ldup.c -o $@

clean:
		rm -fr gmon.out *.o ext/*.o a.out *~ *.a *.dSYM session* $(PROG)

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

bgzf.o: bgzf.h
hts.o: bgzf.h hts.h kseq.h khash.h ksort.h
ldup.o: sam.h bgzf.h hts.h kdq.h khash.h
sam.o: sam.h bgzf.h hts.h khash.h kseq.h kstring.h
trim.o: kvec.h kstring.h khash.h kseq.h
