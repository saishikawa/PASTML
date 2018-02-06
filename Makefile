
PRG    = PASTML
OBJ    = main.o runpastml.o make_tree.o lik.o marginal_lik.o marginal_approxi.o output_tree.o output_states.o fletcher.o nrutil.o golden.o fletcherJC.o

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
runpastml.o : runpastml.c pastml.h marginal_lik.h lik.h marginal_approxi.h golden.h fletcher.h fletcherJC.h make_tree.h
make_tree.o : make_tree.c pastml.h
lik.o : lik.c pastml.h
marginal_lik.o : marginal_lik.c pastml.h
marginal_approxi.o : marginal_approxi.c pastml.h
output_tree.o : output_tree.c pastml.h
output_states.o : output_states.c pastml.h
fletcher.o : fletcher.c pastml.h nrutil.h
nrutil.o : nrutil.c pastml.h nrutil.h
golden.o : golden.c pastml.h
fletcherJC.o : fletcherJC.c pastml.h nrutil.h
