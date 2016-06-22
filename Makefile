CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function
PROG=adna-trim
OBJS=kthread.o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

adna-trim:$(OBJS) trim.o
		$(CC) $(CFLAGS) $(OBJS) trim.o -o $@ -lz -lm -lpthread

clean:
		rm -fr gmon.out *.o ext/*.o a.out *~ *.a *.dSYM session* $(PROG)

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE
