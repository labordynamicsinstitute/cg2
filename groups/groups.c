/****************************************
 *
 *   Making groups of mobility
 *   
 *   Ported by :        Amine Ouazad
 *   Original authors : Robert Creecy, Lars Vilhuber
 *   
 ****************************************/

#include <stdio.h>
#include <string.h>

/****************************************
 *
 *   Data declaration
 *
 ****************************************/

typedef char BOOL;

#define TRUE  1
#define FALSE 0
#define ASSERT(x) {if (!(x)) {printf("Assert failed\n"); exit_program(1);} }


typedef struct {
  long pupil;
  long school;
} record;

#define STACK_SCHOOL  1
#define STACK_PUPIL   2

typedef struct {
  long id;
  char type; // can be either STACK_PUPIL or STACK_SCHOOL
} stack_element;

long group;    // Current group explored

long npupils;  // Number of different pupils 
long nschools; // Number of different schools 
long nobs;     // Number of observations

record *byp;  // Records by pupil
record *bys;  // Records by school

stack_element *m;     // Stack for exploration
long mpoint;

long *pg; // Group of pupils
long *sg; // Group of school

long *pindex;  // Points to the last record of the pupil
long *sindex;  // Points to the last record of the school

BOOL *ptraced;   // Pupil traced
BOOL *straced;   // School traced
BOOL *ponstack;  // Pupil on stack
BOOL *sonstack;  // School on stack

long ntraced; // Just for keeping track of how many have been traced

void
exit_program(int exitcode) {
  printf("Freeing memory\n");
  free(byp);      free(bys);
  free(m);  
  free(pg);       free(sg);
  free(pindex);   free(sindex);
  free(ptraced);  free(straced);
  free(ponstack); free(sonstack);

  if (exitcode != 0) {
    printf("Non normal exit , exit code %u\n", exitcode);
    exit(exitcode);
  }
}

BOOL
check_data()
{

  long i; // To scan the whole array
  long thispupil, thisschool; // the current pupil or school
  
  thispupil = 1;
  thisschool = 1;

  for (i=0; i<nobs; i++) {
    if (byp[i].pupil != thispupil) {
      if (byp[i].pupil != thispupil + 1) {
	fprintf(stderr,"Error - by pupil file not sorted or missing sequence, record %u\n",i);
	fprintf(stderr,"This pupil : pupilid %u schoolid %u \n Previous pupil : pupilid %u schoolid %u\n",
		byp[i].pupil, byp[i].school, byp[i-1].pupil, byp[i-1].school);
	exit_program(1);
      }
      thispupil++;
    }

    if (bys[i].school != thisschool) {
      if (bys[i].school != thisschool + 1) {
	fprintf(stderr,"Error - by school file not sorted or missing sequence, record %u\n", i);
	fprintf(stderr,"This pupil : pupilid %u schoolid %u \n Previous pupil : pupilid %u schoolid %u\n",
		bys[i].pupil, bys[i].school, bys[i-1].pupil, bys[i-1].school);
	exit_program(1);
      }
      thisschool++;
    }

  }

  printf("Data checked - By pupil and by firm files correctly sorted and sequenced.\n");

  return TRUE;
}


void
read_size() {

  FILE *fp_groupsin;
  char sizestring[100];
  char *token;

  // Reading the size of the data in groups.in

  fp_groupsin = (FILE*)fopen("groups.in","r");

  if (!fp_groupsin) {
    fprintf(stderr,"Unable to open file groups.in\n");
    exit_program(1);
  }
  printf("Opened file groups.in !\n");

  if (!fgets(sizestring, 100, fp_groupsin)) {
    fprintf(stderr,"Unable to read the size of the dataset in groups.in\n");
    exit_program(1);
  }
  printf("Got string in file groups.in\n");

  token = (char*)strtok(sizestring," \n");
  printf("Got token %s\n",token);
  if (!token) {
    fprintf(stderr,"Unable to read groups.in : incorrect format.\n");
    exit_program(1);
  }
  nobs = atol(token);

  token = (char*)strtok(NULL," \n");
  printf("Got token %s\n", token);
  if (!token) {
    fprintf(stderr,"Unable to read groups.in : incorrect format.\n");
    exit_program(1);
  }
  npupils = atol(token);

  token = (char*)strtok(NULL," \n");
  printf("Got token %s\n", token);
  if (!token) {
    fprintf(stderr,"Unable to read groups.in : incorrect format.\n");
    exit_program(1);
  }
  nschools = atol(token);

  fclose(fp_groupsin);
  printf("Size of the problem : %u observations %u pupils %u schools\n", nobs, npupils, nschools);

}

void
read_datasets() {

  /***  Reads data from cellsbyp.txt and cellsbys.txt
	and stores it in byp and bys */
 
  // Initial declarations

  FILE *fp_cellsbyp, *fp_cellsbys;
  char recordstring[100];
  long j;
  char *token;

  // Starting with cellsbyp.txt

  printf("Reading dataset cellsbyp.txt\n");

  fp_cellsbyp = (FILE*)fopen("cellsbyp.txt","r");
  if (!fp_cellsbyp) {
    fprintf(stderr,"Unable to open file cellsbyp.txt\n");
    exit_program(1);
  }
  
  if (!fgets(recordstring,100,fp_cellsbyp)) {
    fprintf(stderr,"Unable to read line of labels cellsbyp.txt\n");
    exit_program(1);
  }
  printf("First line of cellsbyp.txt : %s\n",recordstring);
  
  for (j = 0; j<nobs; j++) {
    
    if (!fgets(recordstring,100,fp_cellsbyp)) {
      fprintf(stderr,"Unable to read %u th record of cellsbyp.txt\n",j+1);
      exit_program(1);
    }
    
    token = (char*)strtok(recordstring,",\n");
    if (!token) {
      fprintf(stderr,"Unable to read first token of %u th record\n",j+1);
      exit_program(1);
    }
    byp[j].pupil = atol(token);

    token = (char*)strtok(NULL,",\n");
    if (!token) {
      fprintf(stderr,"Unable to read second token of %u th record\n",j+1);
      exit_program(1);
    }
    byp[j].school = atol(token);

    //printf("Record %u : Pupil %u School %u \n",j,byp[j].pupil, byp[j].school);
  }
  fclose(fp_cellsbyp);
  printf("Successfully read cellsbyp.txt\n");

  // Starting with cellsbys.txt

  printf("Reading datasets cellsbys.txt\n");

  fp_cellsbys = (FILE*)fopen("cellsbys.txt","r");
  if (!fp_cellsbys) {
    fprintf(stderr,"Unable to open file cellsbys.txt\n");
    exit_program(1);
  }
  
  if (!fgets(recordstring,100,fp_cellsbys)) {
    fprintf(stderr,"Unable to read line of labels cellsbys.txt\n");
    exit_program(1);
  }
  printf("First line of cellsbys.txt : %s\n",recordstring);
  
  for (j = 0; j<nobs; j++) {
    
    if (!fgets(recordstring,100,fp_cellsbys)) {
      fprintf(stderr,"Unable to read %u th record of cellsbys.txt\n",j+1);
      exit_program(1);
    }
    
    token = (char*)strtok(recordstring,",\n");
    if (!token) {
      fprintf(stderr,"Unable to read first token of %u th record\n",j+1);
      exit_program(1);
    }
    bys[j].pupil = atol(token);

    token = (char*)strtok(NULL,",\n");
    if (!token) {
      fprintf(stderr,"Unable to read second token of %u th record\n",j+1);
      exit_program(1);
    }
    bys[j].school = atol(token);

    //printf("Record %u : Pupil %u School %u \n",j,byp[j].pupil, byp[j].school);

  }
  fclose(fp_cellsbys);
  
  printf("Successfully read cellsbys.txt\n");

}

void 
write_datasets() {
  // Writes a dataset with the group of each record

  FILE *fp_groupfile;
  long j;
  char record[100];

  fp_groupfile = (FILE*) fopen("groups.out","w");
  if (!fp_groupfile) {
    fprintf(stderr,"Unable to create groups.out for group output\n");
    exit_program(1);
  }
  
  // Checks consistency of output data & writes group data to the file
  
  fputs("pupilid,schoolid,group\n",fp_groupfile);

  for (j=0; j<nobs; j++) {
    if (   (sg[byp[j].school-1] != pg[byp[j].pupil - 1]) 
	|| (sg[bys[j].school-1] != pg[bys[j].pupil - 1])) {
      fprintf(stderr,"Inconsistent output ...\n");
      exit_program(1);
    }

    sprintf(record, "%u,%u,%u\n", byp[j].pupil, byp[j].school, pg[byp[j].pupil-1]);
    if (fputs(record,fp_groupfile)==EOF) {
      fprintf(stderr,"Error writing to file groups.out\n");
      exit_program(1);
    }
  }

  fclose(fp_groupfile);

  printf("File groups.out written\n");
}

void
explore () 
// Explore schools and pupils in the stack
{
  long thisschool, thispupil, aschool, person ;
  long lower, upper;

  ASSERT(mpoint>0);
  ASSERT(group>0);

  //printf("Beginning exploration, with group %u\n", group);

  if (m[mpoint-1].type == STACK_SCHOOL) {
    // The element at the top of the stack is a school


    // Remove the element and store schoolid in "thisschool"
    thisschool = m[mpoint-1].id;
    ASSERT(thisschool>0);
    m[mpoint-1].type = 0;
    mpoint -- ;
    sonstack[thisschool-1] = FALSE;

    // Add this school to the current group
    sg[thisschool-1] = group;
    if (group > 1) {
      printf("Added school %u to group %u\n",thisschool, group);
    }
    straced[thisschool-1] = TRUE;
    ntraced ++;

    // Add all pupils of the school to the group and ...
    if (thisschool == 1) {
      lower = 0;
    }
    else {
      lower = sindex[thisschool-1 -1] + 1;
    }
    upper = sindex[thisschool - 1];

    ASSERT(upper<nobs);
    ASSERT(lower>=0);

    for (person = lower ; person <= upper ; person ++) {
      thispupil = bys[person].pupil;
      ASSERT((thispupil >= 0) && (thispupil<=npupils));
      pg[thispupil-1] = group;
      ASSERT(group>0);

      // ... to the stack if not already traced nor on the stack
      if (ptraced[thispupil-1] == FALSE && ponstack[thispupil-1] == FALSE) {
	  
	//printf("Adding pupil %u to the stack.\n",thispupil);

	mpoint ++;
	m[mpoint - 1].id = thispupil;
	m[mpoint - 1].type = STACK_PUPIL;
	ponstack[thispupil-1] = TRUE;
	ASSERT(mpoint < (npupils + nschools));
      }
    }
  } 
  else if (m[mpoint-1].type == STACK_PUPIL) {
    // The element at the top of the stack is a pupil

    // Remove the element from the top of the stack
    thispupil = m[mpoint - 1 ].id;
    mpoint --;
    // Add this pupil to the current group
    pg[ thispupil - 1 ] = group;
    ptraced[ thispupil - 1 ] = TRUE;
    ntraced ++ ;
    ponstack[ thispupil - 1 ] = FALSE;

    // Add all schools of the pupil to the group and ...
    
    if (thispupil == 1) {
      lower = 0;
    } else {
      lower = pindex[thispupil - 2] + 1;
    }
    upper = pindex[thispupil - 1];

    for (aschool = lower ; aschool <= upper ; aschool ++) {
      thisschool = byp[aschool].school;
      sg[thisschool-1]= group;
    // ... to the stack if not already traced nor on the stack
      if (straced[thisschool-1] == FALSE && sonstack[thisschool] == FALSE) {
	  
	//printf("Adding school %u to the stack.\n",thisschool);
	mpoint ++;
	m[mpoint-1].id = thisschool;
	m[mpoint-1].type = STACK_SCHOOL;
	sonstack[thisschool -1 ] = TRUE;
      }
    }
  } 
  else {
    printf("Strange element on the stack, with type %u\n",m[mpoint-1].type);
  }
}

int
main()
{

  long j;

  printf("Grouping schools & pupils by mobility.\n");

  read_size();

  /*   Memory Allocation    */
  printf("Memory allocation\n");
  
  byp = (record*)calloc(nobs,sizeof(record));
  bys = (record*)calloc(nobs,sizeof(record));
  m   = (stack_element*)calloc(npupils+nschools,sizeof(stack_element));
  pg  = (long*)calloc(npupils,sizeof(long));
  sg  = (long*)calloc(nschools,sizeof(long));
  pindex = (long*)calloc(npupils,sizeof(long));
  sindex = (long*)calloc(nschools,sizeof(long));
  ptraced = (BOOL*)calloc(npupils,sizeof(BOOL));
  straced = (BOOL*)calloc(nschools,sizeof(BOOL));
  ponstack = (BOOL*)calloc(npupils,sizeof(BOOL));
  sonstack = (BOOL*)calloc(nschools,sizeof(BOOL));
  
  if ( !byp || !bys || !m || !pg || !sg 
       || !pindex || !sindex 
       || !ptraced || !straced || !ponstack) {
    fprintf(stderr, "Memory allocation error.\n");
    exit_program(1);
  }
 
  /* Initializing memory */
  /*  memset(ptraced, FALSE, sizeof(BOOL)*npupils);
  memset(straced, FALSE, sizeof(BOOL)*nschools);
  memset(ponstack, FALSE, sizeof(BOOL)*npupils);
  memset(sonstack, FALSE, sizeof(BOOL)*nschools);
  memset(m, 0, sizeof(stack_element)*(npupils+nschools));  */

  /* Read datasets */
  read_datasets();

  if (check_data()==FALSE) {
    fprintf(stderr,"Data incorrectly sorted.\n");
    exit_program(1);
  }


  for (j = 1; j < nobs; j++) {
    pindex[byp[j].pupil  - 1] = j;
    sindex[bys[j].school - 1] = j;
  }


  long nextschool;
  group = 1;
  ntraced = 0;
  
  nextschool = 1;
  
  // Add school 1 to the stack for exploration
  mpoint = 1;
  m[mpoint-1].id = 1;
  m[mpoint-1].type = STACK_SCHOOL;
  sonstack[1-1] = TRUE;
  
  while (mpoint > 0) {
    
    // Explore schools and pupils on the stack
    explore();
    
    if (mpoint == 0) {
      // New group : stack exhausted
      
      group ++;
      printf("New group created, group number %u\n", group);
      
      // Look through the dataset if there is a school without a group
      while (nextschool < nschools && sg[nextschool-1] != 0) {
	nextschool ++;
      }
      
      if (sg[nextschool-1] == 0) {
	// There is a school without a group
	// Add it to the stack
	
	printf("Adding school %u to the stack.\n",nextschool);
	printf("State : mpoint %u / group %u\n",mpoint, group);
	mpoint = 1;
	m[mpoint-1].id = nextschool;
	m[mpoint-1].type = STACK_SCHOOL;
	sonstack[nextschool-1] = TRUE;
      }
      
    }
    
  }
  
  printf("Traced %u\n",ntraced);
  printf("Groups defined, writing datasets\n");

  /* Write datasets */
  write_datasets();
  
  exit_program(0);
  return 0;
}
