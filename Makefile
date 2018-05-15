
PRG    = PASTML
OBJ    = main.o likelihood.o logger.o marginal_approximation.o marginal_likelihood.o param_minimization.o parsimony.o runpastml.o scaling.o states.o tree.o

CFLAGS = -mcmodel=medium -w -O2
LFLAGS = -lm -lgsl

CC     =  gcc $(CFLAGS)

$(PRG) : $(OBJ) 
	$(CC) -o $@ $^ $(LFLAGS)

.c.o:
	$(CC) -c $<

clean:
	rm -rf $(PRG) $(OBJ)

main.o : main.c pastml.h runpastml.h
runpastml.o : runpastml.c pastml.h marginal_likelihood.h likelihood.h marginal_approximation.h param_minimization.h states.h scaling.h logger.h tree.h parsimony.h
likelihood.o : likelihood.c pastml.h scaling.h
logger.o : logger.c pastml.h
marginal_approximation.o : marginal_approximation.c pastml.h likelihood.h
marginal_likelihood.o : marginal_likelihood.c pastml.h scaling.h
param_minimization.o : param_minimization.c pastml.h likelihood.h logger.h
parsimony.o : parsimony.c pastml.h
scaling.o : scaling.c pastml.h likelihood.h
states.o : states.c pastml.h
tree.o : tree.c pastml.h logger.h

