PRG    = PASTML
SRC    = main.c runpastml.c tree.c likelihood.c marginal_likelihood.c marginal_approximation.c states.c param_minimization.c scaling.c logger.c parsimony.c
OBJ    = $(SRC:.c=.o)

CFLAGS = -mcmodel=medium -w -Ic:/gsl
LFLAGS = -Lc:/gsl -lm -lgsl -lgslcblas

CC     =  gcc $(CFLAGS)

$(PRG) : $(OBJ)
	$(CC) -o $@ $^ $(LFLAGS)

.c.o: 
	$(CC) -c $<

clean:
	rm -rf $(PRG) $(OBJ)
	
main.o : main.c pastml.h runpastml.h
runpastml.o : runpastml.c pastml.h marginal_likelihood.h likelihood.h marginal_approximation.h param_minimization.h states.h logger.h tree.h parsimony.h
tree.o : tree.c pastml.h logger.h
likelihood.o : likelihood.c pastml.h scaling.h
marginal_likelihood.o : marginal_likelihood.c pastml.h scaling.h
marginal_approximation.o : marginal_approximation.c pastml.h likelihood.h
states.o : states.c pastml.h
param_minimization.o: param_minimization.c likelihood.h logger.h
scaling.o: scaling.c likelihood.h pastml.h
logger.o: logger.c pastml.h
parsimony.o: parsimony.c pastml.h