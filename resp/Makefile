include ../config.h

PROGRAM = resp 
BINDIR = ../exe

OBJS= amopen.o lapack.o mexit.o resp.o  
SRCS= amopen.F lapack.F mexit.F resp.F

.F.o: 
	$(FC) -c $(FFLAGS) -o $@ $<

all:		$(PROGRAM)

$(PROGRAM):	$(OBJS) $(LIBS)
		$(FC) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

install:	$(PROGRAM)
		mv $(PROGRAM) $(BINDIR)

clean:
		rm -f $(OBJS)
