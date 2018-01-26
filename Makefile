
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

main.o : main.c asrml.h runpastml.h
runpastml.o : runpastml.c asrml.h marginal_lik.h lik.h marginal_approxi.h golden.h fletcher.h fletcherJC.h make_tree.h
make_tree.o : make_tree.c asrml.h
lik.o : lik.c asrml.h
marginal_lik.o : marginal_lik.c asrml.h
marginal_approxi.o : marginal_approxi.c asrml.h
output_tree.o : output_tree.c asrml.h
output_states.o : output_states.c asrml.h
fletcher.o : fletcher.c asrml.h nrutil.h
nrutil.o : nrutil.c asrml.h nrutil.h
golden.o : golden.c asrml.h
fletcherJC.o : fletcherJC.c asrml.h nrutil.h
