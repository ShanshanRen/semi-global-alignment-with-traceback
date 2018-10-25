#include<iostream>
#include <limits>
#include<cmath>
#include<algorithm> 
#include<list> 
#include<string.h>
#include<stdio.h>
using namespace std;

int total_number;
int total_match;
int total_insertion;
int total_deletion;
int total_soft;
double diff(timespec start, timespec end)
{
	double a=0;
    if((end.tv_nsec-start.tv_sec)<0)
    {
        a=end.tv_sec-start.tv_sec-1;
        a+=(1000000000+end.tv_nsec-start.tv_nsec)/1000000000.0;
    }
    else
    {
        a=end.tv_sec-start.tv_sec+(end.tv_nsec-start.tv_nsec)/1000000000.0;

    }
    return a;
}

struct Paramter
{
	int match;
	int mismatch;
	int open;
	int extend;
};


struct InputData
{
char read_base[1000];
char reference_base[1000];
};

struct Element
{
	char a;
	int length;
};

enum State 
{
    MATCH,
    INSERTION,
    DELETION,
    CLIP
};

int wd(char a,char b)
{
	return a==b?200:-150;
}

void matrix(char * read, char * ref,int col,int row,int sw[][600],int btrack[][600],int mt[][600])
{
  
    Paramter p;
    p.match=200;
    p.mismatch=-150;
    p.open=-260;
    p.extend=-11;
    int MATRIX_MIN_CUTOFF=(int) -1e8;
	int lowInitValue=numeric_limits<int>::min()/2;
    for(int i=0;i<row+1;i++)
    	for(int j=0;j<col+1;j++)
    		sw[i][j]=btrack[i][j]=mt[i][j]=0;

    int best_gap_v[600],gap_size_v[600],best_gap_h[600],gap_size_h[600];
    for(int i=0;i<col+1;i++)
    {
    	best_gap_v[i]=lowInitValue;
    	gap_size_v[i]=0;
    }
   
	
    for(int i=0;i<row+1;i++)
    {
    	best_gap_h[i]=lowInitValue;
    	gap_size_h[i]=0;
    }
    int *match_h_2=mt[0];
    int * curRow=sw[0];
    for(int i=1;i<row+1;i++)
    {
    	char a_base=ref[i-1];
    	int * lastRow=curRow;
	int *match_h_1=match_h_2;
	match_h_2=mt[i];
    	curRow=sw[i];
    	int * curBtrack=btrack[i];

    	for(int j=1;j<col+1;j++)
    	{
    		char b_base=read[j-1];
    		int step_diag=lastRow[j-1]+wd(a_base,b_base);

    		int prev_gap=lastRow[j]+p.open;
    		best_gap_v[j]+=p.extend;
    		if(prev_gap>best_gap_v[j])
    		{
    			best_gap_v[j]=prev_gap;
    			gap_size_v[j]=1;
    		}
    		else
    			gap_size_v[j]++;
    		int step_down=best_gap_v[j];
    		int kd=gap_size_v[j];


    		prev_gap=curRow[j-1]+p.open;
    		best_gap_h[i]+=p.extend;
    		if(prev_gap>best_gap_h[i])
    		{
    			best_gap_h[i]=prev_gap;
    			gap_size_h[i]=1;
    		}
    		else
    			gap_size_h[i]++;
    		int step_right=best_gap_h[i];
    		int ki=gap_size_h[i];
		
		match_h_2[j]=0;
    		bool diag=(step_diag>=step_down)&&(step_diag>=step_right);
    		if(diag)
    		{
    			curRow[j]=max(MATRIX_MIN_CUTOFF, step_diag);
    			curBtrack[j]=0;
			match_h_2[j]=match_h_1[j-1]+1;
			if(j==1||i==1)
			match_h_2[j]=0;
    		}
    		else
    			if(step_right>=step_down)
    			{
    				curRow[j]=max(MATRIX_MIN_CUTOFF, step_right);
    				curBtrack[j]=-ki;// negative = horizontal
    			}
    			else
    			{
    				curRow[j]=max(MATRIX_MIN_CUTOFF, step_down);
    				curBtrack[j]=kd; // positive=vertical
    			}
    	}
    }
/*
for(int j=1;j<col+1;j++)
{
  for(int i=1;i<row+1;i++)
	printf("%d ", mt[i][j]);
printf("\n\n");
}
*/
}

Element makeElement(State state,int length)
{
	char a;
	Element b;
	switch(state)
	{
		case MATCH: b.a='M'; break;
		case INSERTION:b.a='I';break;
		case DELETION:b.a='D'; break;
		case CLIP:b.a='S'; break;
	}
	b.length=length;
	return b;
} 


void calculateCigar(int col,int row,int sw[][600],int btrack[][600],int mt[][600])
{
	int maxscore=numeric_limits<int>::min();
	int segment_length=0;

	int p1=0,p2=0;

	p2=col;
    for(int i=1;i<row+1;i++)  
    {
        int curScore = sw[i][col];
//	printf("%d ",curScore);
        if (curScore >= maxscore ) 
        {       		
            p1 = i;
            maxscore = curScore;
        }
    }
//	printf("\n");
 
//	printf("%d %d\n",p1,p2);

    int * bottomRow=sw[row];
    for ( int j = 1 ; j < col+1; j++) 
    {
	    int curScore=bottomRow[j];
	                    // data_offset is the offset of [n][j]
//	    printf("%d ",curScore);
	    if ( curScore > maxscore ||(curScore == maxscore && abs(row-j) < abs(p1 - p2) ) ) 
	    {
	        p1 = row;
	        p2 = j ;
	        maxscore = curScore;
	        segment_length = col - j ; // end of sequence 2 is overhanging; we will just record it as 'M' segment
        }
     }

    list<Element> lce;
    if ( segment_length > 0) 
    {
        lce.push_back(makeElement(CLIP, segment_length));
        segment_length = 0;
    }

    State state = MATCH;
    do 
    {
        int btr = btrack[p1][p2];
        State new_state;
        int step_length = 1;
        if ( btr > 0 ) 
        {
            new_state = DELETION;
            step_length = btr;
        } 
        else if ( btr < 0 ) 
		        {
		            new_state = INSERTION;
		            step_length = (-btr);
		        } 
		        else 
		        	{
				new_state = MATCH; // and step_length =1, already set above
				step_length=mt[p1][p2];
				if(step_length==0)
					step_length=1;
				}
//	printf("%d %d %d %d \n",state,new_state,p1,p2);
            // move to next best location in the sw matrix:
        switch( new_state ) 
        {
                case MATCH:
			if(step_length>0)   
			{	p1-=step_length; p2-=step_length; }
			else
				{ p1--;p2--;}
			break; // move back along the diag in the sw matrix
                case INSERTION: p2 -= step_length; break; // move left
                case DELETION:  p1 -= step_length; break; // move up
        }
//	printf("%d %d %d \n",new_state,p1, p2);
            // now let's see if the state actually changed:
        if ( new_state == state ) segment_length+=step_length;//this will only happen when misMatch or Match!
        else 
        {
                // state changed, lets emit previous segment, whatever it was (Insertion Deletion, or (Mis)Match).
            lce.push_back(makeElement(state, segment_length));
            segment_length = step_length;
            state = new_state;
        }
        // next condition is equivalent to  while ( sw[p1][p2] != 0 ) (with modified p1 and/or p2:
    } while ( p1 > 0 && p2 > 0 );


	int alignment_offset;
    lce.push_back(makeElement(state, segment_length));
    if ( p2 > 0 ) lce.push_back(makeElement(CLIP, p2));
    alignment_offset = p1;
    lce.reverse();
   

    //print lce
 //  printf("%d %d\n",maxscore, alignment_offset);
  //  printf("[");
    int i=0;
    for (list<Element>::iterator it=lce.begin(); it != lce.end(); ++it)
    {
  //      if(i>0)
    //        printf(", ");
      // printf("%d%c", it->length,it->a); 
	if(it->a=='M')  total_match++;
	if(it->a=='I')  total_insertion++;
	if(it->a=='D')  total_deletion++;
	
	if(it->a=='S')  total_soft++;
      // i++;
    }
 //   printf("]");
//	printf("\n");
 //   printf("%d\n",alignment_offset);
    
}


void smithwaterman(int size,InputData * inputdata)
{
    struct timespec start,finish;   
	for(int i=0;i<size;i++)
	{
		int col=strlen(inputdata[i].read_base);
		int row=strlen(inputdata[i].reference_base);
        int sw[600][600];
        int btrack[600][600];
	int mt[600][600];
        //printf("%d\n",i);
		matrix(inputdata[i].read_base,inputdata[i].reference_base,col,row,sw,btrack,mt);
		calculateCigar(col,row,sw,btrack,mt);
	}
}



int main(int argc,char * args[])
{
	total_match=0;
	total_insertion=0;
	total_soft=0;
	total_deletion=0;
	 timespec start,finish;
	total_number=0;	
	FILE * file;
	file=fopen(args[1],"r");

 	int size;
 	fscanf(file,"%d",&size);
 	//printf("size=%d\n",size);
 	double computation_time=0;
	 while(!feof(file))
	{
 		InputData * inputdata=(InputData* )malloc(size*(sizeof(InputData)));
 		for(int i=0;i<size;i++)
   		 {
    			fscanf(file,"%s ",inputdata[i].reference_base);
    			fscanf(file,"%s ",inputdata[i].read_base);
		  total_number++;
   	  /*	if(total_number==19475)
  		 {
		printf("%s\n", inputdata[i].reference_base);
		printf("%s\n", inputdata[i].read_base);
  			return 0;
		}
	 */	
		
  		}
 		clock_gettime(CLOCK_MONOTONIC_RAW,&start); 	
		smithwaterman(size,inputdata);
 		clock_gettime(CLOCK_MONOTONIC_RAW,&finish);
                computation_time+=diff(start,finish); 

		fscanf(file,"%d",&size);
	}
	printf("%d %d %d %d %d\n",total_match+total_insertion+total_deletion,total_match,total_insertion,total_deletion,total_soft);
//	printf("computation_time=%e\n",computation_time);
	return 0;
}
