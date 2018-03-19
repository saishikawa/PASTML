
PRG    = PASTML
OBJ    = main.o runpastml.o make_tree.o likelihood.o marginal_likelihood.o joint_likelihood.o marginal_approximation.o output_tree.o output_states.o output_simulation.o param_minimization.o scaling.o logger.o eigen.o models.o

CFLAGS = -mcmodel=medium -w
LFLAGS = -lm -lgsl

CC     =  gcc $(CFLAGS)

$(PRG) : $(OBJ) 
	$(CC) -o $@ $^ $(LFLAGS)

.c.o:
	$(CC) -c $<

clean:
	rm -rf $(PRG) $(OBJ)

main.o : main.c pastml.h runpastml.h
runpastml.o : runpastml.c pastml.h marginal_likelihood.h likelihood.h marginal_approximation.h param_minimization.h scaling.h make_tree.h logger.h joint_likelihood.h output_states.h output_tree.h output_simulation.h
make_tree.o : make_tree.c pastml.h
likelihood.o : likelihood.c pastml.h
marginal_likelihood.o : marginal_likelihood.c pastml.h
joint_likelihood.o : joint_likelihood.c pastml.h
marginal_approxi.o : marginal_approxi.c pastml.h
logger.o : logger.c pastml.h
scaling.o : scaling.c pastml.h
output_tree.o : output_tree.c pastml.h
output_states.o : output_states.c pastml.h
output_simulation.o : output_simulation.c pastml.h
param_minimization.o : param_minimization.c pastml.h
eigen.o : eigen.c pastml.h
models.o : models.c pastml.h
