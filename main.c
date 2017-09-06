#include "asrml.h"
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

  Tree *s_tree;
  Node *root;

void *check_alloc(int nbrelt, int sizelt){
  void *retval;
  if( (retval=calloc(nbrelt,sizelt)) != NULL ) {
    return retval; 
  }
  printf("Not enough memory\n");
  exit(0);
}

void Generic_Exit(const char *file, int line, const char *function, int code){
  fprintf(stderr,"\n== Err. in file '%s' (line %d), function '%s'\n",file,line,function);
  exit(code);
}

unsigned int tell_size_of_one_tree(char* filename) {
	/* the only purpose of this is to know about the size of a treefile (NH format) in order to save memspace in allocating the string later on */
	/* wew open and close this file independently of any other fopen */
	unsigned int mysize = 0;
	char u;
	FILE* myfile = fopen(filename, "r");
	if (myfile) {
		while ( (u = fgetc(myfile))!= ';' ) { /* termination character of the tree */
			if (u == EOF) break; /* shouldn't happen anyway */
			if (isspace(u)) continue; else mysize++;
		}
		fclose(myfile);
	} /* end if(myfile) */
	return (mysize+1);
}	

int copy_nh_stream_into_str(FILE* nh_stream, char* big_string) {
	int index_in_string = 0;
	char u;
	/* rewind(nh_stream); DO NOT go to the beginning of the stream if we want to make this flexible enough to read several trees per file */
	while ( (u = fgetc(nh_stream))!= ';' ) { /* termination character of the tree */
		if (u == EOF) { big_string[index_in_string] = '\0'; return 0; } /* error code telling that no tree has been read properly */
		if (index_in_string == MAX_TREELENGTH - 1) {
		  fprintf(stderr,"Fatal error: tree file seems too big, are you sure it is an NH tree file? Aborting.\n");
		  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
		}
		if (isspace(u)) continue;
		big_string[index_in_string++] = u; 
	}
	big_string[index_in_string++] = ';';
	big_string[index_in_string] = '\0';
	return 1; /* leaves the stream right after the terminal ';' */
} /*end copy_nh_stream_into_str */

void free_edge(Edge* edge) {
	int i;
	if (edge == NULL) return;
	for (i=0; i<2; i++) if(edge->subtype_counts[i]) free(edge->subtype_counts[i]);
	
	free(edge);
}

void free_node(Node* node, int count, int num_anno) {
        int j;

	if (node == NULL) return;
	//printf("free name %d\n", count);
	if (node->name&&count!=0) free(node->name);
	//printf("free comment %d\n", count);
	if (node->comment) free(node->comment);
	//printf("free neigh %d\n", count);	
	free(node->neigh);
	//printf("free br %d\n", count);
	free(node->br);

        for(j=0;j<num_anno;j++) free(node->pij[j]);
        free(node->pij);
        if(count==0){
          for(j=0;j<num_anno;j++) free(node->rootpij[j]);
          free(node->rootpij);
        }
        free(node->prob);
        free(node->best_i);
        free(node->condlike);
        free(node->condlike_mar);
        free(node->sortedlike);
        free(node->sortedstates);
        free(node->marginal);
        free(node->mar_state);
        free(node->tmp_best);
        free(node->mar_prob);
        free(node->up_like);
        free(node->sum_down);
        free(node->calc_flag);

	free(node);
}

void free_tree(Tree* tree, int num_anno) {
	if (tree == NULL) return;
	int i;
        //printf("free 1\n");
	for (i=0; i < tree->nb_nodes; i++) free_node(tree->a_nodes[i],i,num_anno);
        //printf("free 2\n");
	for (i=0; i < tree->nb_edges; i++) free_edge(tree->a_edges[i]);
        //printf("free 3\n");
	for (i=0; i < tree->nb_taxa; i++) free(tree->taxa_names[i]);
        //printf("free 4\n");

	free(tree->taxa_names);
	free(tree->a_nodes);
	free(tree->a_edges);
	free(tree);

}

//MAIN

int main(int argc, char** argv){
  /*input (-a annotation -t tree -c characters -n tips -x fraction -m model -f frequency)*/
  int i,ii,j,ret,line=0,check,check2,sum_freq=0, retcode;
  int *states, *count, *factors, *locked;
  double sum=0.,mu=0.,lnl=0.,scale,maxlnl=0.,sum_marginal,factor,nano,sec,upbound=10.0,border_frac;
  double *frequency,*root_prob;  
  char **annotations, **character, **tips, **tipnames, filefreq[100], filescale[100], *c_tree, str[] = "int", str2[]="ROOT";
  int num_anno=0,num_tips=0;
  char *annotation_name, *tree_name, *model, *scaling, *keep_ID;
  FILE *treefile,*annotationfile,*ffreq,*fscale;
  struct timespec samp_ini, samp_fin;
  int opt, check_freq=0;
  int tmpstate[4];
  double tmpfreq;
  char tmpchar[5];

  opterr = 0;
  while ((opt = getopt(argc, argv, "a:t:x:m:f:s:I:")) != -1) {

        switch (opt) {
            case 'a':
                annotation_name=optarg;
                //printf("-a %s\n",annotation_name);
                break;

            case 't':
                tree_name=optarg;
                //printf("-t %s\n",tree_name);
                break;

            case 'x':
                border_frac=atof(optarg);
                //printf("-x %lf\n",border_frac);
                break;

            case 'm':
                model=optarg;
                //printf("-m %s\n",model);
                break;

            case 's':
                scaling=optarg;
                //printf("-m %s\n",model);
                break;

            case 'I':
                keep_ID=optarg;
                //printf("-m %s\n",model);
                break;

            case 'f':
                if(num_anno==0){
                  printf("Error: define -c option in advance of -f\n");
                  exit(0);
                }
                frequency=check_alloc(num_anno,sizeof(double));
                //printf("optind = %s\n", argv[optind-1]);
                for(i=0;i<num_anno;i++){
                  frequency[i]=atof(argv[optind-1+i]);
                  //printf("freq %d = %lf\n", i, frequency[i]);
                }
                check_freq=1;
                break;

            default: /* '?' */
                printf("Unexpected options...\nUsage: %s [-a name_of_annotationfile] [-t name_of_treefile] [-c number_of_states(integer)] [-n number_of_tips(integer)] [-x denominator_of_fraction(integer)] [-m name_of_model] [-f frequencies(double_array)]\n", argv[0]);
                break;
        }
  }
  printf("\n*** Information of primary annotations ***\n\n");
  if(check_freq==0) frequency=calloc(MAXCHAR,sizeof(double));
  states=calloc(MAXNSP,sizeof(int));
  count=calloc(MAXCHAR,sizeof(int));
  factors=calloc(MAXNSP,sizeof(int));
  annotations=calloc(MAXNSP,sizeof(char*));
  tips=calloc(MAXNSP,sizeof(char*));
  locked=calloc(MAXNSP,sizeof(int));
  for(i=0;i<MAXNSP;i++){
    annotations[i]=calloc(MAXLNAME,sizeof(char));
    tips[i]=calloc(MAXLNAME,sizeof(char));
  }
  character=calloc(MAXCHAR,sizeof(char*));
  for(i=0;i<MAXCHAR;i++){
    character[i]=calloc(MAXLNAME,sizeof(char));
  }
  if(border_frac==0){
    border_frac=0.0;
  } else {
    border_frac=1.0/border_frac;
  }
  //clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&samp_ini);
  srand((unsigned) time(NULL));
  /*Read annotation from file*/
  annotationfile=fopen(annotation_name, "r");
  if(!annotationfile){
      printf("Cannot find annotation file : %s\n", annotation_name);
      exit(0);
  }
  num_anno=0;
  while( ( ret = fscanf( annotationfile, "%[^,],%s\n", tips[line], annotations[line]) ) != EOF ){
    check=0;
    for(i=0;i<=line;i++){
      if(strcmp(annotations[line],annotations[i])==0){
        states[line]=states[i];
         strcpy(character[states[line]],annotations[line]);
         break;
      } else {
         check2=1;
         for(j=0;j<i;j++){
           if(strcmp(annotations[i],annotations[j])==0) check2=0;
         }
         if(check2==1)check++;
         states[line]=check;
         if(check > num_anno) num_anno=check;
         //strcpy(character[states[line]],annotations[line]);
      }
    }
    count[states[line]]++;
    //printf("%s, annotation %s = %d\n", tips[line], annotations[line], states[line]);
    line++;
    if(line==num_tips){
      //printf("Finish to read the annotation file...\n");
      break;
    }
  }
  fclose(annotationfile);
  //printf("number of character is %d\n",num_anno+1);
  num_anno++;

  for(i=0;i<num_anno;i++){
    sum_freq=sum_freq+count[i];
  }
  printf("\n*** Frequency of each character in the MODEL %s ***\n\n",model);
  for(i=0;i<num_anno;i++){
    if(strcmp(model,"JC")==0){
      frequency[i]=(double)1/num_anno;
    } else if (strcmp(model,"F81_E")==0){
      frequency[i]=(double)count[i]/sum_freq;
    }
    printf("%s = %lf\n", character[i], frequency[i]);
  }

  /*Read tree from file*/

  treefile = fopen(tree_name,"r");
  if (treefile == NULL) {
    fprintf(stderr,"File %s not found or impossible to access media. Aborting.\n", tree_name);
    Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
  }

  unsigned int treefilesize = 3 * tell_size_of_one_tree(tree_name);
  if (treefilesize > MAX_TREELENGTH) {
    fprintf(stderr,"Tree filesize for %s bigger than %d bytes: are you sure it's a valid NH tree? Aborting.\n", tree_name, MAX_TREELENGTH/3);
    Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
  }

  c_tree = (char*) check_alloc(treefilesize+1, sizeof(char)); 
  retcode = copy_nh_stream_into_str(treefile, c_tree);
  if (retcode != 1) { 
    fprintf(stderr,"Unexpected EOF while parsing the reference tree! Aborting.\n"); 
    Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
  }
  fclose(treefile);

  /*Make Tree structure*/
  s_tree  = complete_parse_nh(c_tree, num_anno, keep_ID); /* sets taxname_lookup_table en passant */  
  root = s_tree->a_nodes[0];
  num_tips=s_tree->nb_taxa;
  //printf("number of taxa is %d\n",num_tips);

  //Compute initial likelihood
  for(i=0;i<num_anno;i++){
    sum+=frequency[i]*frequency[i];
  }
  mu=1/(1-sum);
  calc_lik(root, tips, states, num_tips, num_anno, mu, 1.0, model, frequency, &lnl);
  printf("\n*** Initial likelihood of the tree ***\n\n %lf\n",lnl);
  //Optimise likelihood and tree scale with the golden section search
  scale=1.0;
  golden(tips, states, num_tips, num_anno, mu, upbound, model, frequency, &scale);

  if(strcmp(scaling,"F")==0){
    scale=1.0;
  }

  calc_lik(root, tips, states, num_tips, num_anno, mu, scale, model, frequency, &maxlnl);
  printf("\n*** Optimised tree scaling factor ***\n\n %lf\n\n*** Optimised likelihood ***\n\n %.10f\n",scale,maxlnl);

  //Joint reconstruction with Pupko's method
  printf("\n*** Calculating Joint Solution...\n");
  joint(root, tips, states, num_tips, num_anno, frequency);
  printf("\n*** Joint solution, Likelihood ***\n\n %lf\n", s_tree->pupko_value);

  //Marginal likelihood calculation
  printf("\n*** Calculating Posterior Probabilities...\n");
  down_like_marginal(root, num_tips, num_anno, mu, scale, frequency);

  printf("\n*** Predicting the best ancestral states...\n\n");
  make_samples(tips,states,num_tips,num_anno,character,mu,scale,frequency,model,border_frac);

  //free all
  printf("Freeing step...\n");
  free(annotations);
  free(tips);
  free(character);
  free(states); free(count); free(factors); free(frequency), free(locked);
  free_tree(s_tree, num_anno);
  
  //clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&samp_fin);
  //sec=(double)(samp_fin.tv_sec - samp_ini.tv_sec);
  //nano=(double)(samp_fin.tv_nsec - samp_ini.tv_nsec)/1000.0/1000.0/1000.0;

  //printf("\n*** Total execution time ***\n\n %3.10lf seconds\n\n", (sec+nano));
  exit(0);
}
