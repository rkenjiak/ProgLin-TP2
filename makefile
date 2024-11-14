#@author Rodrigo Kenji Asato Kobayashi @RGA 202319040134 
#GLPK=~/opt
TRACE=NDEBUG

LOADLIBS=-L $(GLPK)/lib -lglpk -lm
cflags= -c -D_REENTRANT -g -Wall -I $(GLPK)/include  -D$(TRACE)

compile = gcc

program = caminho

csources = $(program).c

cobjects = $(csources:.c=.o)

$(program): $(cobjects)
	$(compile) -o $(program) $(cobjects) $(LOADLIBS)

#$(cobjects): 
.c.o:
	$(compile) -o $@ $*.c $(cflags)

clean:
	rm -f *.o $(program)
