
PRG    = PASTML
OBJ    = main.o runpastml.o make_tree.o likelihood.o marginal_likelihood.o marginal_approximation.o output_tree.o output_states.o param_minimization.o scaling.o

CFLAGS = -mcmodel=medium -w
LFLAGS = -lm

CC     =  gcc $(CFLAGS)

$(PRG) : $(OBJ) 
	$(CC) -o $@ $^ $(LFLAGS)

.c.o:
	$(CC) -c $<

clean:
	rm -rf $(PRG) $(OBJ)

main.o : main.c pastml.h runpastml.h
runpastml.o : runpastml.c pastml.h marginal_likelihood.h likelihood.h marginal_approximation.h param_minimization.h scaling.h make_tree.h
make_tree.o : make_tree.c pastml.h
lik.o : lik.c pastml.h
marginal_lik.o : marginal_lik.c pastml.h
marginal_approxi.o : marginal_approxi.c pastml.h
output_tree.o : output_tree.c pastml.h
output_states.o : output_states.c pastml.h
