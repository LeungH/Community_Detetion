#include<stdio.h>
#include<stdlib.h>
#include <time.h>


char input_path[1000];
char output_path[1000];
char label_real_or_generated[50];
int popsize=100;
int T=20;
int generation = 50;
int verbose = 2;
double alpha = 1.0;
double pc = 0.7;
double pm = 1.0;


void check_parameters() {
    if (T%2==0) T=T+1;
    if (pm>1) pm = 1.0;
    if (pc>1) pc = 1.0;
    if (T>popsize) {
	fprintf(stderr,"Parameter error: T(%d)>popsize(%d)\n",T,popsize);
	exit(-1);
    }
    if (popsize<=0) {
	fprintf(stderr,"Parameter error: popsize(%d)<=0\n",popsize);
	exit(-1);
    }
    if (popsize>10000) {
	fprintf(stderr,"Parameter error: popsize(%d)>10000\n",popsize);
	exit(-1);
    }
    if (T<=0) {
	fprintf(stderr,"Parameter error: T(%d)<=0\n",T);
	exit(-1);
    }
    if (generation<0) {
	fprintf(stderr,"Parameter error: generation(%d)<0\n",generation);
	exit(-1);
    }
}


void str_chop(char *source,char *target){

    int i = 0;
    while(source[i] != '\0'){
        target[i] = source[i];
        i++;
    }

}

void str_seperate(char line[], char first[],char second[]){

    int num = 0,i = 0,j = 0;
    int flag = 0;
    while(line[num] != '\n'){
        if(line[num] == ' '){
            flag = 1;
            num ++;
        }
        if(flag == 0){
            first[i] = line[num];
            i++; num++;
        }
        if(flag == 1){
            second[j] = line[num];
            j++;num++;
        }

    }
    first[i] = '\0';
    second[j] = '\0';

}

int para_set(char *flags){
    int temp;
    int maxnumofeveryline = 200;
    char whole_line[600],as[100],as_value[500];
    FILE *fp;

    if((fp = fopen(flags,"r")) == NULL){
        printf("open flags.dat failed\n");
        return -1;
    }
    while(fgets(whole_line, maxnumofeveryline,fp) != NULL){
        str_seperate(whole_line,as,as_value);

        if(as[1] == 'i') str_chop(as_value,input_path);
        else if(as[1] == 'r') str_chop(as_value,output_path);
        else if(as[1] == 'p' && as[2] == 'o') popsize = atoi(as_value);
        else if(as[1] == 'p' && as[2] == 'c') pc = atof(as_value);
        else if(as[1] == 'p' && as[2] == 'm') pm = atof(as_value);
        else if(as[1] == 'T') T = atoi(as_value);
        else if(as[1] == 'g') generation = atoi(as_value);
        else if(as[1] == 'v') verbose = atoi(as_value);
        else if(as[1] == 'a') alpha = atof(as_value);
		else if(as[1] == 'q') str_chop(as_value,label_real_or_generated);
        else printf("Please check flags.dat\n");
    }
    return 1;
}

void get_parameters(int argc,char* argv[]) {
    int ai;
    double af;
    char as[1000];
    int i,m = 0;
    char flags[20];

    unsigned int au,seed = time(NULL);

    if(argv[1][1] == 'f'){
        if(sscanf(argv[2],"%s",&flags) == 1){
            printf("get parameters from flags.dat\n");
            m = para_set(flags);
            if(m == -1) printf("please check flags.dat\n");
            }

        else{
            printf("Please check parameter");
            }
    }


    else{
    for (i=1; i!=argc; i=i+2) {
//	if (argv[i][0]=='-') {
	    switch (argv[i][1]) {

        case 'i':
            if(sscanf(argv[i+1],"%s",&input_path) == 1) printf("input_path ok\n");
            else{
                printf("Please check parameters");
            }
            break;

        case 'r':
            if(sscanf(argv[i+1],"%s",&output_path) == 1) printf("output_path ok\n");
            else{
                printf("Please check parameters\n");
            }
            break;

	case 'p':
            if(argv[i][2] == 'o'){
                if(sscanf(argv[i+1],"%d",&ai)) popsize = ai;
		//printf("popsize okkkkkkkk: %d ",ai);
            }
            else if(argv[i][2] == 'c'){
                if(sscanf(argv[i+1],"%lf",&af)) pc = af;
            }
            else if(argv[i][2] == 'm'){
                if(sscanf(argv[i+1],"%lf",&af)) pm = af;
            }
            else{
                printf("Please check parameter\n");
		    }
		    break;

	case 'T':
		    if (sscanf(argv[i+1],"%d",&ai)==1) T = ai;
		    else {
			printf("Please check parameter\n");
		    }
		    break;

        case 'g':
		    if (sscanf(argv[i+1],"%d",&ai)==1) generation = ai;
		    else {
			printf("Please check parameter\n");
		    }
		    break;

        case 'v':
		    if (sscanf(argv[i+1],"%d",&ai)==1) verbose = ai;
		    else {
			printf("Please check parameter\n");
		    }
		    break;


		case 'a':
		    if (sscanf(argv[i+1],"%lf",&af)==1) alpha = af;
		    else {
			printf("Please check parameter\n");
		    }
		    break;

        case 'q':
            if(sscanf(argv[i+1],"%s",&label_real_or_generated) == 1)
				printf("label_real_or_generated ok\n");
            else{
                printf("Please check parameter");
            }
            break;

		default:
		    fprintf(stderr,"Invalid options : %s\n",argv[i]);
		    exit(-1);
	    }
	//}

  //  fprintf(stderr,"Too many input files\n");
  //  exit(-1);
	}
    }

    //test
    printf("input_path: %s\n",input_path);
    printf("output_path: %s\n",output_path);
    printf("popsize = %d\n",popsize);
    printf("T = %d\n",T);
    printf("generation = %d\n",generation);
    printf("verbose = %d\n",verbose);
    printf("alpha = %lf\n",alpha);
    printf("pc = %lf\n",pc);
    printf("pm = %lf\n",pm);
	printf("label_real_or_generated: %s\n",label_real_or_generated);

    check_parameters();

    srand(seed);
    if (verbose==0) return;

    printf("CMOEA_SN Algorithm\n");
    printf("input_path :%s\n",input_path);
    printf("output_path:%s\n",output_path);
    printf("popsize    : %d\n",popsize);
    printf("T          : %d\n",T);
    printf("generation : %d\n",generation);
    printf("alpha      : %lf\n",alpha);
    printf("pc         : %lf\n",pc);
    printf("pm         : %lf\n",pm);
    printf("label_real_or_generated :%s\n",label_real_or_generated);
    printf("-----------------------------------\n");
}
