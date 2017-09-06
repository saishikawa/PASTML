
PRG    = ASTRAL
OBJ    = main.o make_tree_alloc.o lik_poly.o golden.o joint_poly.o marginal_poly.o relax_poly_PP.o samp_poly.o output_poly.o summary.o

CFLAGS = -mcmodel=medium -w
LFLAGS = -lm

CC     =  gcc $(CFLAGS)

$(PRG) : $(OBJ) 
	$(CC) -o $@ $^ $(LFLAGS)

.c.o:
	$(CC) -c $<

clean:
	rm -rf $(PRG) $(OBJ)

main.o : main.c asrml.h
make_tree_alloc.o : make_tree_alloc.c asrml.h
lik_poly.o : lik_poly.c asrml.h
golden.o : golden.c asrml.h
joint_poly.o: joint_poly.c asrml.h
marginal_poly.o : marginal_poly.c asrml.h
relax_poly_PP.o : relax_poly_PP.c asrml.h
samp_poly.o : samp_poly.c asrml.h
output_poly.o : output_poly.c asrml.h
summary.o : summary.c asrml.h
