#
#Implementation of Vantage point tree pthreads
#*Doinakis Michail
#*e-mail: doinakis@ece.auth.gr
#
SHELL	:=	/bin/bash
CC	=	gcc-7
CFLAGS	=	-O3	-fcilkplus	-fopenmp	-Wall -g -Werror
INCLUDES	=
LDFLAGS	=
LIBS	= -lm
SRC	=	./src/vptree
TYPES	=	sequential 	pthreads	openmp	cilk

MAIN	=	./src/main
.PRECIOUS:	%.a

all:	$(addprefix	$(MAIN)_,	$(TYPES))
		mv $(addsuffix .a, $(addprefix $(SRC)_, $(TYPES))) ./lib/


lib:	$(addsuffix	.a,	$(addprefix	$(SRC)_,$(TYPES)))
		mv $(addsuffix .a, $(addprefix $(SRC)_, $(TYPES))) ./lib/

$(MAIN)_%:	$(MAIN).c	$(SRC)_%.a
	$(CC)	$(CFLAGS)	$(INCLUDES)	-o	$@	$^	$(LDFLAGS)	$(LIBS)

.o.a:
	ar	rcs	$@	$<
.c.o:
	$(CC)	$(CFLAGS)	$(INCLUDES)	-o	$@	-c	$<

clean:
	$(RM)	*.o	*~	$(addprefix	$(MAIN)_,	$(TYPES))	vptree_*.a
